"""Signed-profile semantics and structural model controls on synthetic graphs."""
import sys
import unittest
from pathlib import Path
import numpy as np
import torch
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'analysis'))
from crcsc_drug_gnn import normalize_signed,CMSDrugGNN,score_profiles
from crc_experiments import weighted_adjacency

class CMSDrugTests(unittest.TestCase):
    def data(self):
        adj=lambda a,b,n,m:weighted_adjacency(a,b,np.ones(len(a)),n,m)
        signed=np.array([[1,0,0],[0,-1,0],[0,0,1],[.5,0,.5]],dtype=np.float32)
        return {'genes':range(3),'adj_ppi':adj([0,1,1,2],[1,0,2,1],3,3),'adj_drug':adj([0,1],[0,0],2,3),'adj_disease':adj([0],[1],1,3),'adj_gpath':adj([0,1,2],[0,0,0],3,1),'adj_pathg':adj([0,0,0],[0,1,2],1,3),'target_mask':np.array([True,True]),'profiles':{'signed':signed,'magnitude':np.abs(signed)}}

    def test_signed_normalization_preserves_sign_and_empty_rows(self):
        np.testing.assert_allclose(normalize_signed([[-1,3],[0,0]]),[[-.25,.75],[0,0]])
        with self.assertRaises(ValueError):normalize_signed([[np.nan,1]])

    def test_zero_swap_and_identical_profile_controls(self):
        data=self.data();model=CMSDrugGNN(data,use_pathways=True)
        scores,_=score_profiles(model,data)
        np.testing.assert_allclose(scores['swapped'],scores['signed'][[1,0,3,2]])
        np.testing.assert_allclose(scores['zero'][0],scores['zero'][1])
        np.testing.assert_allclose(scores['identical'][0],scores['identical'][3])

    def test_identical_drug_neighborhoods_remain_tied(self):
        data=self.data();model=CMSDrugGNN(data,use_pathways=True)
        scores,_=score_profiles(model,data)
        np.testing.assert_allclose(scores['signed'][:,0],scores['signed'][:,1])

    def test_drug_only_ignores_signed_subtype_context(self):
        data=self.data();model=CMSDrugGNN(data,use_context=False)
        scores,_=score_profiles(model,data)
        np.testing.assert_allclose(scores['signed'][0],scores['signed'][1])

if __name__=='__main__':unittest.main()
