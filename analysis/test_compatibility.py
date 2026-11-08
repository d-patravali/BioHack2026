import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from compatibility import bh, compare_sets, map_genes, pathway_audit, exploratory_enrichment, SUBTYPES


class CompatibilityTests(unittest.TestCase):
    def test_bh_preserves_original_order(self):
        np.testing.assert_allclose(bh([.04, .001, .03, 1]), [.0533333333, .004, .0533333333, 1])

    def test_symbol_ambiguity_and_same_id_aliases(self):
        prime = pd.DataFrame([
            ["1", "11", "gene/protein", "A", "NCBI"],
            ["2", "11", "gene/protein", "A", "NCBI"],
            ["3", "12", "gene/protein", "B", "NCBI"],
            ["4", "13", "gene/protein", "B", "NCBI"],
        ], columns=["index", "id", "type", "name", "source"])
        eco = pd.DataFrame([["gene:"+s, "Gene", s] for s in ["A", "B", " A ", "a", "patient_id"]],
                           columns=["~id", "~label", "symbol:String"])
        result = map_genes(eco, prime)
        self.assertEqual(result.status.tolist(), ["exact_symbol", "ambiguous", "normalized_whitespace_symbol",
                                                  "unmatched", "excluded_metadata_artifact"])
        self.assertEqual(result.primekg_id.tolist(), ["11", "", "11", "", ""])

    def test_pathway_label_equality_is_not_stable_id_match(self):
        edges = pd.DataFrame({"~to": ["pathway:REACTOME_APOPTOSIS"]})
        prime = pd.DataFrame([["1", "R-HSA-109581", "pathway", "Apoptosis", "REACTOME"]],
                             columns=["index", "id", "type", "name", "source"])
        with tempfile.TemporaryDirectory() as directory:
            audit, summary = pathway_audit(edges, prime, Path(directory))
        self.assertEqual(audit.status.iloc[0], "unverified_label_candidate")
        self.assertEqual(summary.confident_matches.iloc[0], 0)
        self.assertTrue(pd.isna(summary.biological_overlap_pct.iloc[0]))

    def test_overlap_uses_sets(self):
        sets = dict(zip(SUBTYPES, [{"A", "B"}, {"B", "C"}, {"B"}, set()]))
        result = compare_sets(sets, "test")
        self.assertEqual(result.overlap.iloc[0], 1)
        self.assertAlmostEqual(result.jaccard.iloc[0], 1/3)

    def test_tiny_groups_cannot_pass_burden_screen(self):
        x = pd.DataFrame({"A": [1, 1, 0, 0, 0, 0, 0, 0]})
        labels = pd.Series([s for s in SUBTYPES for _ in range(2)])
        mapping = pd.DataFrame({"ecograph_label": ["A"], "primekg_id": ["1"]})
        targets = pd.DataFrame({"gene_id": ["1"]})
        with tempfile.TemporaryDirectory() as directory:
            result, summary = exploratory_enrichment(x, labels, mapping, targets, Path(directory))
        self.assertFalse(result.burden_screen_eligible.any())
        self.assertTrue(result.burden_cmh_one_sided_p.eq(1).all())
        self.assertFalse(result.screen_candidate.any())
        self.assertTrue(summary.joint_exploratory_screen_candidates.isna().all())


if __name__ == "__main__":
    unittest.main()
