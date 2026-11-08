"""Explicitly selected read-only backends; no global backend switching."""
import json
from functools import cached_property

from .config import PROJECT_ROOT
from .repository import ResearchRepository, DataUnavailable, read_rows

BACKENDS = ('ecograph', 'crcsc')
CMS_SUBTYPES = ('CMS1', 'CMS2', 'CMS3', 'CMS4')
CMS_RUN_ID = '20261108_173948'
CMS_SOURCE_ID = '20261108_172638'


def validate_backend(backend):
    if not isinstance(backend, str) or backend not in BACKENDS:
        raise ValueError('Choose ecograph or crcsc')
    return backend


def catalog(backend):
    validate_backend(backend)
    if backend == 'crcsc':
        return {'id': backend, 'label': 'CRCSC · CMS expression', 'subtypes': list(CMS_SUBTYPES),
                'gnn_condition': 'ppi_pathway', 'gnn_label': 'PPI + pathway GNN',
                'profile_type': 'signed expression signature', 'recorded_available': False}
    return {'id': backend, 'label': 'ECoGraph · GI genomic (original)',
            'subtypes': ['GI.CIN', 'GI.HM-indel', 'GI.HM-SNV', 'GI.GS'],
            'gnn_condition': 'target_only', 'gnn_label': 'Target-only GNN',
            'profile_type': 'mutation observation frequencies', 'recorded_available': True}


class CRCSCRepository(ResearchRepository):
    """Shares PrimeKG evidence tools; overrides only subtype/model semantics."""
    def __init__(self, root=PROJECT_ROOT):
        super().__init__(root)
        self.run = self.root / 'reports/crcsc_drug_gnn' / CMS_RUN_ID
        self.source_run = self.root / 'reports/crcsc' / CMS_SOURCE_ID

    @cached_property
    def manifest(self):
        try:
            result = json.loads((self.run / 'manifest.json').read_text())
        except (OSError, ValueError) as error:
            raise DataUnavailable('Reviewed CRCSC manifest is unavailable') from error
        if not isinstance(result, dict) or result.get('input_hashes_unchanged') is not True or result.get('cms_training_labels') != 0:
            raise DataUnavailable('Unexpected CRCSC experiment manifest')
        return result

    def normalize_subtype(self, value):
        value = value.strip()
        return next((s for s in CMS_SUBTYPES if s.casefold() == value.casefold()), None)

    def get_subtypes(self, disease):
        normalized = disease.strip().lower().replace(' ', '_').replace('-', '_')
        if normalized != 'colorectal_cancer':
            return {'status': 'unsupported_disease', 'disease': normalized, 'subtypes': [],
                    'supported_diseases': ['colorectal_cancer'], 'data_backend': 'crcsc'}
        return {'status': 'ok', 'disease': normalized, 'subtypes': list(CMS_SUBTYPES),
                'data_backend': 'crcsc', 'classification': 'expression-defined consensus molecular subtypes',
                'source': self.source(self.run / 'cms_nodes.csv'),
                'warning': 'CMS1–4 are not equivalent to ECoGraph GI genomic subtypes. Cached hypothesis scores only.'}

    def get_model_evaluation(self):
        manifest = self.manifest
        path = self.run / 'evaluation_metrics.csv'
        rows = read_rows(path)
        grouped = {}
        for row in rows:
            grouped.setdefault(row['condition'], []).append(float(row['test_ap']))
        chosen = next(r for r in rows if r['condition'] == 'ppi_pathway' and r['seed'] == '42')
        return {'status': 'ok', 'data_backend': 'crcsc', 'run_id': CMS_RUN_ID, 'source': self.source(path),
                'task': 'held-out general disease indication recovery against sampled unknowns; no CMS therapeutic labels',
                'test_known_indications': manifest['split']['test']['positives'],
                'test_unknown_pairs': manifest['split']['test']['unknowns'],
                'metric': 'average precision; not accuracy, response probability, or efficacy',
                'conditions': {k: {'runs': len(v), 'mean_test_ap': sum(v)/len(v), 'min_test_ap': min(v), 'max_test_ap': max(v)} for k,v in grouped.items()},
                'selected_cached_model': {'condition': 'ppi_pathway', 'seed': 42, 'profile': 'signed', 'test_ap': float(chosen['test_ap'])},
                'interpretation': 'Popularity outperformed every fitted GNN. Molecular context gave a small indication-AP gain over drug-only, but the leading CMS drug group remains shared and unvalidated.',
                'subtype_efficacy_evaluated': False, 'model_uses_phenotype_edges': False,
                'warning': 'Signed CMS expression is transferred through a decoder trained on unsigned disease-gene associations. Drug-action direction and subtype transfer validity are not established.'}

    @cached_property
    def candidates(self):
        return read_rows(self.run / 'candidate_hypotheses.csv.gz')

    def get_drug_predictions(self, subtype, top_k=10):
        subtype = self.normalize_subtype(subtype)
        if subtype is None:
            return {'status': 'unsupported_subtype', 'data_backend': 'crcsc', 'supported_subtypes': list(CMS_SUBTYPES)}
        if type(top_k) is not int or not 1 <= top_k <= 20:
            return {'status': 'invalid_input', 'message': 'top_k must be an integer from 1 to 20'}
        selected = [r for r in self.candidates if r['condition'] == 'ppi_pathway' and r['seed'] == '42'
                    and r['profile'] == 'signed' and r['subtype'] == subtype and int(r['tied_rank_min']) <= top_k]
        if not selected:
            raise DataUnavailable('No reviewed CRCSC cached prediction rows')
        selected.sort(key=lambda r: (int(r['tied_rank_min']), r['drug_id']))
        predictions = [{'drug_id': r['drug_id'], 'drug_name': r['drug_name'],
                        'rank_min': int(r['tied_rank_min']), 'tie_group_size': int(r['tie_group_size']),
                        'raw_logit': float(r['raw_logit']), 'recorded_target_ncbi_ids': r['recorded_target_ncbi_ids'].split('|'),
                        'training_indication_count': int(r['training_indication_count'])} for r in selected]
        return {'status': 'ok', 'data_backend': 'crcsc', 'subtype': subtype, 'backend': 'cached_gnn_experiment',
                'fresh_inference': False, 'run_id': CMS_RUN_ID, 'model_condition': 'ppi_pathway',
                'seed': 42, 'profile': 'signed', 'source': self.source(self.run / 'candidate_hypotheses.csv.gz'),
                'requested_top_k': top_k, 'returned_count': len(predictions), 'ties_included': True,
                'predictions': predictions, 'literature_review_performed': False,
                'warning': 'Full cutoff ties retained; no uniquely best drug within a tie. Raw logits are not efficacy probabilities. Signed expression is not a drug-action direction or dependency. No literature novelty or therapeutic validity established.',
                'evaluation': self.get_model_evaluation()}

    @cached_property
    def expression_edges(self):
        return read_rows(self.run / 'cms_expression_edges.csv')

    def get_subtype_genes(self, subtype, limit=20, gene=''):
        subtype = self.normalize_subtype(subtype)
        if subtype is None:
            return {'status': 'unsupported_subtype', 'data_backend': 'crcsc', 'supported_subtypes': list(CMS_SUBTYPES)}
        if type(limit) is not int or not 1 <= limit <= 50:
            return {'status': 'invalid_input', 'message': 'limit must be an integer from 1 to 50'}
        rows = [r for r in self.expression_edges if r['x_id'] == 'crcsc:GSE39582:' + subtype]
        if not rows:
            raise DataUnavailable('No reviewed CMS expression profile rows')
        total = len(rows); gene = gene.strip()
        if gene:
            rows = [r for r in rows if (r['y_id'] == gene if gene.isdigit() else self.nodes[r['y_index']]['name'] == gene)]
        rows.sort(key=lambda r: (-abs(float(r['signed_weight'])), r['y_id']))
        genes = [{'gene_id': r['y_id'], 'namespace': 'NCBI', 'symbol': self.nodes[r['y_index']]['name'],
                  'gene_node_index': int(r['y_index']), 'expression_relation': r['relation'],
                  'signed_weight': float(r['signed_weight']), 'discovery_adjusted_effect': float(r['discovery_adjusted_effect']),
                  'discovery_q': float(r['discovery_q'])} for r in rows[:limit]]
        cohort = next(r for r in read_rows(self.source_run / 'cohort_counts.csv') if r['CMS'] == subtype)
        return {'status': 'ok', 'data_backend': 'crcsc', 'subtype': subtype, 'genes': genes,
                'total_genes': total, 'gene_filter': gene, 'matching_genes': len(rows), 'returned_count': len(genes),
                'truncated': len(rows) > limit, 'profile': 'signed', 'profile_type': 'expression_signature',
                'cohort_samples': {'discovery': int(cohort['discovery']), 'validation': int(cohort['validation'])},
                'source': self.source(self.run / 'cms_expression_edges.csv'),
                'evidence_type': 'discovery_adjusted_expression_association',
                'warning': 'These are relative higher/lower expression associations, not mutation frequencies, causal dependencies, activation scores, or drug-effect directions. No match means absent from the selected signature, not biological absence.'}


def make_repository(backend='ecograph', root=PROJECT_ROOT):
    validate_backend(backend)
    return CRCSCRepository(root) if backend == 'crcsc' else ResearchRepository(root)
