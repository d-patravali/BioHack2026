"""Identifier ambiguity, discovery-only normalization, and adjusted contrast checks."""
import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'analysis'))
from crcsc_experiment import parse_gene_id, standardized, adjusted_contrast, bh, top_signature, sample_id

class CRCSCTests(unittest.TestCase):
    def test_mapping_never_picks_one_of_multiple_genes(self):
        self.assertEqual(parse_gene_id('780///100616237')[1],'ambiguous')
        self.assertEqual(parse_gene_id('780///780'),('780','single_ncbi_id'))
        self.assertEqual(parse_gene_id('---')[1],'unmapped')
        self.assertEqual(parse_gene_id('NCBI:780')[1],'invalid_identifier')
        self.assertEqual(sample_id('GSM12345_test.CEL'),'GSM12345')

    def test_scaling_uses_discovery_only_and_excludes_constants(self):
        train=np.array([[1.,8.],[3.,8.]])
        test=np.array([[102.,9.]])
        a,b,valid=standardized(train,test)
        np.testing.assert_equal(valid,[True,False])
        self.assertAlmostEqual(float(a.mean()),0.)
        self.assertAlmostEqual(float(b[0,0]),100/np.sqrt(2))

    def test_adjustment_recovers_effect_with_batch_confounding(self):
        group=np.array([0,0,0,1,0,1,1,1.])
        batch=np.array([0,0,1,0,1,0,1,1.])
        design=np.column_stack([np.ones(8),batch])
        values=(2*group+10*batch)[:,None]
        beta,p,df=adjusted_contrast(values,group,design)
        self.assertAlmostEqual(beta[0],2.)
        self.assertEqual(df,5)
        self.assertLess(p[0],1e-10)

    def test_bh_selection_uses_only_supplied_discovery_effects(self):
        np.testing.assert_allclose(bh([.01,.04,.03]),[.03,.04,.04])
        np.testing.assert_equal(top_signature(np.array([1.,9.,3.]),np.array([.01,.9,.02]),1),[2])

if __name__=='__main__':unittest.main()
