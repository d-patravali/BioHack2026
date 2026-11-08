"""Small semantic tests for weighted pooling, profile selection and ablation."""

import sys
import unittest
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
if any(importlib.util.find_spec(name) is None for name in ["torch","sklearn"]):
    raise unittest.SkipTest("CRC diagnostics require optional modeling dependencies: torch and scikit-learn")
import torch

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"analysis"))
from crc_experiments import weighted_adjacency,make_profiles,IndicationGNN
from crc_report import tie_top_mask


class CRCExperimentTests(unittest.TestCase):
    def test_weighted_pooling_and_empty_rows(self):
        adjacency=weighted_adjacency([0,0,1],[0,1,1],[1,3,2],3,2)
        np.testing.assert_allclose(adjacency.to_dense().numpy(),[[.25,.75],[0,1],[0,0]])
        with self.assertRaises(ValueError):
            weighted_adjacency([0],[0],[-1],1,1)

    def test_selective_profile_excludes_singletons(self):
        matrix=pd.DataFrame({"A":[1]*40,"B":[1,0,0,0,0,0,0,0,0,0]*4})
        labels=pd.Series([s for s in ["a","b","c","d"] for _ in range(10)])
        mapped=pd.DataFrame({"node_index":[10,11]},index=["A","B"])
        _,profiles,_,_=make_profiles(matrix,labels,pd.Index([10,11]),mapped)
        np.testing.assert_allclose(profiles["selective_frequency"],[[1,0]]*4)
        np.testing.assert_allclose(profiles["frequency"][0],[1/1.1,.1/1.1],rtol=1e-6)

    def test_drug_only_ignores_disease_identity(self):
        ppi=weighted_adjacency([0,1],[1,0],[1,1],2,2)
        drugs=weighted_adjacency([0],[0],[1],1,2)
        disease=weighted_adjacency([0,1],[0,1],[1,1],2,2)
        model=IndicationGNN(2,ppi,drugs,disease,use_context=False)
        scores=model.score(model.encode(),np.array([[0,0],[1,0]])).detach().numpy()
        np.testing.assert_allclose(scores[0],scores[1])

    def test_tied_top_k_includes_entire_group(self):
        self.assertEqual(tie_top_mask([3,2,2,2,1],2).tolist(),[True,True,True,True,False])


if __name__=="__main__":
    unittest.main()
