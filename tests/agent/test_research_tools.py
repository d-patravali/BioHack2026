"""Offline contracts: ties, provenance, graph roles, aliases and missing data."""

import csv
import gzip
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import httpx2

from agent.repository import ResearchRepository, RUN_ID, query
from agent.config import Settings, NEBIUS_BASE_URL
from agent.runtime import create_client, create_agent, run_turn, DebugHooks


def csv_file(root, name, fields, rows):
    path=root/name
    path.parent.mkdir(parents=True,exist_ok=True)
    opener=gzip.open if path.suffix==".gz" else open
    with opener(path,"wt",newline="",encoding="utf-8") as handle:
        writer=csv.DictWriter(handle,fieldnames=fields)
        writer.writeheader();writer.writerows(rows)


class ResearchToolTests(unittest.IsolatedAsyncioTestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name);self.repo=ResearchRepository(self.root)
        nodes=[
            {"index":"1","id":"2908","name":"NR3C1","type":"gene/protein","source":"NCBI"},
            {"index":"2","id":"2908","name":"NR3C1_alias","type":"gene/protein","source":"NCBI"},
            {"index":"3","id":"10","name":"SAME","type":"gene/protein","source":"NCBI"},
            {"index":"4","id":"11","name":"SAME","type":"gene/protein","source":"NCBI"},
            {"index":"5","id":"DB00596","name":"Fixture compound","type":"drug","source":"DrugBank"},
            {"index":"6","id":"123","name":"Fixture phenotype","type":"effect/phenotype","source":"HPO"},
        ]
        csv_file(self.root,"data/processed/primekg/nodes.csv",list(nodes[0]),nodes)
        csv_file(self.root,"reports/compatibility/primekg_direct_targets.csv",["drug_id","gene_id","display_relation"],
            [{"drug_id":"DB00596","gene_id":"2908","display_relation":"target"},
             {"drug_id":"DB00596","gene_id":"10","display_relation":"enzyme"}])
        csv_file(self.root,"data/processed/primekg/phenotype_disease_gene_edges.csv",["x_index","y_index","relation"],
            [{"x_index":"1","y_index":"6","relation":"phenotype_protein"},
             {"x_index":"6","y_index":"2","relation":"phenotype_protein"}])
        self.run=f"reports/crc_experiments/{RUN_ID}/"
        metrics=[{"condition":name,"seed":"42","test_ap":score} for name,score in
                 [("target_only",".8"),("drug_only",".81"),("drug_popularity",".89")]]
        csv_file(self.root,self.run+"evaluation_metrics.csv",list(metrics[0]),metrics)
        prediction={"seed":"42","profile":"frequency","subtype":"subtype:GI.CIN",
            "drug_id":"DB00596","drug_name":"Fixture compound","tied_rank_min":"1","tie_group_size":"2",
            "raw_logit":"2.5","recorded_target_ncbi_ids":"2908","recorded_indication_for_crc_lexical_candidate":"False"}
        csv_file(self.root,self.run+"tie_aware_top_candidates.csv",list(prediction),
            [prediction,{**prediction,"drug_id":"DB00002","drug_name":"Second fixture"}])
        profile={"profile":"frequency","subtype":"subtype:GI.CIN","gene_node_index":"1","symbol":"NR3C1",
            "n_observed":"2","n_profiled":"10","frequency":".2","pooling_weight":".1"}
        csv_file(self.root,self.run+"profile_gene_evidence.csv.gz",list(profile),[profile])

    def test_predictions_are_cached_and_do_not_break_ties(self):
        result=self.repo.get_drug_predictions("GI.CIN",1)
        self.assertFalse(result["fresh_inference"])
        self.assertEqual(result["returned_count"],2)
        self.assertEqual({r["rank_min"] for r in result["predictions"]},{1})
        self.assertEqual(result["evaluation"]["selected_cached_model"]["test_ap"],.8)

    def test_targets_preserve_stable_identity_and_aliases(self):
        result=self.repo.get_drug_targets("DB00596")
        self.assertEqual(result["target_count"],1)
        self.assertEqual(result["targets"][0]["node_indices"],[1,2])
        self.assertEqual(result["relationship"]["display_relation"],"target")

    def test_reciprocal_phenotype_edges_are_deduplicated(self):
        result=self.repo.get_gene_phenotypes("2908")
        self.assertEqual(result["phenotype_count"],1)
        self.assertEqual(result["phenotypes"][0]["namespace"],"HPO")
        self.assertFalse(result["model_uses_phenotype_edges"])

    def test_ambiguous_gene_is_not_silently_merged(self):
        self.assertEqual(self.repo.get_gene_phenotypes("SAME")["status"],"ambiguous_gene")

    def test_gene_filter_queries_full_retained_profile(self):
        result=self.repo.get_subtype_genes("GI.CIN",1,"2908")
        self.assertEqual(result["genes"][0]["observed_samples"],2)
        self.assertEqual(result["genes"][0]["profiled_samples"],10)
        self.assertEqual(self.repo.get_subtype_genes("GI.CIN",1,"missing")["matching_genes"],0)

    def test_invalid_inputs_and_missing_backend_are_distinct(self):
        self.assertEqual(self.repo.get_drug_predictions("unknown")["status"],"unsupported_subtype")
        self.assertEqual(self.repo.get_drug_predictions("GI.CIN",0)["status"],"invalid_input")
        self.assertEqual(self.repo.get_drug_targets("not-a-drug-id")["status"],"invalid_input")
        absent=ResearchRepository(self.root/"absent")
        self.assertEqual(query(absent,"get_drug_predictions",subtype="GI.CIN")["status"],"backend_unavailable")

    async def test_sdk_executes_prediction_and_evidence_tools_offline(self):
        requests=[];logs=[]
        def respond(request):
            requests.append((str(request.url),json.loads(request.content)))
            if len(requests)==1:
                calls=[("get_drug_predictions",{"subtype":"GI.CIN","top_k":1})]
            elif len(requests)==2:
                calls=[("get_drug_targets",{"drug_id":"DB00596"}),
                       ("get_gene_phenotypes",{"gene":"2908","limit":10}),
                       ("get_subtype_genes",{"subtype":"GI.CIN","limit":1,"gene":"2908"})]
            else: calls=[]
            message={"role":"assistant","content":None if calls else "Cached tied hypotheses; graph context is not efficacy evidence."}
            if calls:
                message["tool_calls"]=[{"id":f"call_{len(requests)}_{i}","type":"function",
                    "function":{"name":name,"arguments":json.dumps(arguments)}} for i,(name,arguments) in enumerate(calls)]
            return httpx2.Response(200,json={"id":f"chat_{len(requests)}","object":"chat.completion",
                "created":1,"model":"fixture-model","choices":[{"index":0,"message":message,
                "finish_reason":"tool_calls" if calls else "stop"}],
                "usage":{"prompt_tokens":10,"completion_tokens":10,"total_tokens":20}})
        settings=Settings(api_key="dummy-key",model="fixture-model")
        with patch("agent.tools.repository",return_value=self.repo):
            async with create_client(settings,http_client=httpx2.AsyncClient(transport=httpx2.MockTransport(respond))) as client:
                result=await run_turn(create_agent(settings,client),[],"Investigate the predictions",DebugHooks("dummy-key",logs.append))
        self.assertEqual(len(requests),3)
        self.assertTrue(all(url==NEBIUS_BASE_URL+"chat/completions" for url,_ in requests))
        outputs=[m["content"] for m in requests[-1][1]["messages"] if m["role"]=="tool"]
        self.assertTrue(any("cached_gnn_experiment" in text for text in outputs))
        self.assertTrue(any("recorded_direct_association" in text for text in outputs))
        self.assertTrue(any("cohort_observation_not_validated_enrichment" in text for text in outputs))
        self.assertIn("Cached",result.final_output)
        self.assertEqual(sum("[tool call]" in line for line in logs),4)


if __name__=="__main__":
    unittest.main()
