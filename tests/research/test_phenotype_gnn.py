"""Leakage, relation removal, and ranking checks on small synthetic graphs."""
import sys
import unittest
from pathlib import Path
import numpy as np
import torch
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'analysis'))
from phenotype_gnn import split_labels, sample_unknown, graph_adjacencies, RelationalGNN, evaluate

class PhenotypeGNNTests(unittest.TestCase):
    def test_whole_gene_split_is_reproducible_and_disjoint(self):
        labels={(g,p) for g in range(20) for p in [g%3,(g+1)%3]}
        first=split_labels(labels,17)
        self.assertEqual(first,split_labels(labels,17))
        genes={k:{g for g,p in v} for k,v in first.items()}
        self.assertFalse(genes['train'] & genes['test'])
        self.assertFalse(genes['validation'] & genes['test'])
        self.assertFalse(genes['validation'] & genes['train'])
        self.assertEqual(set().union(*(set(v) for v in first.values())),labels)

    def test_target_reverse_contains_only_training_labels(self):
        ids={'gene/protein':[1,2,3],'effect/phenotype':[4,5],'pathway':[6],'disease':[7]}
        edges={'protein_protein':{(0,1),(1,0)},'pathway_protein':{(0,0),(2,0)},'disease_protein':{(0,0)},'disease_phenotype_positive':{(0,1)}}
        adj=graph_adjacencies(ids,edges,[(0,1)],'molecular')
        self.assertEqual(adj['annotations'][2]._nnz(),1)
        np.testing.assert_equal(adj['annotations_reverse'][2].to_dense().numpy(),[[0,0,0],[1,0,0]])
        model=RelationalGNN({k:len(v) for k,v in ids.items()},adj,width=4)
        with torch.no_grad():
            before=model.encode(removed=['pathways'])
            adj['pathways'][2].values().fill_(999)
            after=model.encode(removed=['pathways'])
        torch.testing.assert_close(before['gene/protein'],after['gene/protein'])

    def test_unknowns_exclude_training_labels(self):
        train=[(0,0),(1,1)]
        unknown=sample_unknown(train,4,np.random.default_rng(4))
        self.assertFalse(set(map(tuple,unknown)) & set(train))

    def test_filtering_and_uniform_ties(self):
        frame=evaluate(np.array([[1.,1.,1.]]),[0],[(0,0)],{(0,0),(0,1)},3)
        self.assertAlmostEqual(frame.mrr.iloc[0],.75)
        self.assertAlmostEqual(frame.hits1.iloc[0],.5)

if __name__=='__main__':unittest.main()
