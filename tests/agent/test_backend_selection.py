"""Backend isolation, CRCSC contracts and SDK routing; no provider calls."""
import ast
import csv
import gzip
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import httpx2
from agents import RunConfig
from agents.tool_context import ToolContext
from starlette.testclient import TestClient

from agent.backends import CRCSCRepository, CMS_RUN_ID, CMS_SOURCE_ID, make_repository
from agent.repository import ResearchRepository, query
from agent.tools import tools_for_backend
from agent.runtime import create_agent, create_client, run_turn, DebugHooks
from agent.config import Settings, NEBIUS_BASE_URL
from ui.server import create_app


def write_csv(root,name,rows):
    path=root/name;path.parent.mkdir(parents=True,exist_ok=True)
    opener=gzip.open if path.suffix=='.gz' else open
    with opener(path,'wt',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


class BackendSelectionTests(unittest.IsolatedAsyncioTestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name);self.repo=CRCSCRepository(self.root)
        self.run=f'reports/crcsc_drug_gnn/{CMS_RUN_ID}/'
        path=self.root/self.run;path.mkdir(parents=True)
        (path/'manifest.json').write_text(json.dumps({'input_hashes_unchanged':True,'cms_training_labels':0,'split':{'test':{'positives':853,'unknowns':853}}}))
        write_csv(self.root,'data/processed/primekg/nodes.csv',[
            {'index':'1','id':'2908','type':'gene/protein','name':'NR3C1','source':'NCBI'},
            {'index':'5','id':'DB00596','type':'drug','name':'Fixture drug','source':'DrugBank'}])
        write_csv(self.root,self.run+'evaluation_metrics.csv',[
            {'condition':c,'seed':'42','test_ap':a} for c,a in [('ppi_pathway','.82'),('ppi','.819'),('drug_only','.818'),('drug_popularity','.894')]])
        candidate={'condition':'ppi_pathway','seed':'42','profile':'signed','subtype':'CMS1','drug_id':'DB00596','drug_name':'Fixture drug','tied_rank_min':'1','tie_group_size':'2','raw_logit':'2.5','training_indication_count':'5','recorded_target_ncbi_ids':'2908'}
        write_csv(self.root,self.run+'candidate_hypotheses.csv.gz',[
            candidate,{**candidate,'drug_id':'DB00002'},
            {**candidate,'profile':'magnitude','raw_logit':'99'},
            {**candidate,'seed':'43','raw_logit':'98'}])
        write_csv(self.root,self.run+'cms_expression_edges.csv',[
            {'x_id':'crcsc:GSE39582:CMS1','y_id':'2908','y_index':'1','relation':'lower_expression_signature','signed_weight':'-.1','discovery_adjusted_effect':'-.7','discovery_q':'.001'}])
        write_csv(self.root,f'reports/crcsc/{CMS_SOURCE_ID}/cohort_counts.csv',[
            {'CMS':'CMS1','discovery':'76','validation':'15'}])

    def test_defaults_remain_original_and_wrong_catalog_rejected(self):
        self.assertIs(type(make_repository(root=self.root)),ResearchRepository)
        self.assertEqual(self.repo.get_drug_predictions('GI.CIN')['status'],'unsupported_subtype')
        self.assertEqual(ResearchRepository(self.root).get_drug_predictions('CMS1')['status'],'unsupported_subtype')
        self.assertEqual(self.repo.get_subtypes('colorectal cancer')['subtypes'],['CMS1','CMS2','CMS3','CMS4'])
        with self.assertRaises(ValueError):make_repository('automatic')

    def test_signed_cached_predictions_keep_ties_and_exact_run(self):
        r=self.repo.get_drug_predictions('cms1',1)
        self.assertEqual(r['data_backend'],'crcsc');self.assertEqual(r['returned_count'],2)
        self.assertEqual({p['raw_logit'] for p in r['predictions']},{2.5})
        self.assertFalse(r['fresh_inference']);self.assertFalse(r['literature_review_performed'])
        self.assertEqual(r['evaluation']['selected_cached_model']['test_ap'],.82)
        self.assertEqual(r['run_id'],CMS_RUN_ID)

    def test_expression_never_masquerades_as_mutation_counts(self):
        r=self.repo.get_subtype_genes('CMS1',1,'2908')
        self.assertEqual(r['profile_type'],'expression_signature')
        self.assertEqual(r['cohort_samples']['discovery'],76)
        gene=r['genes'][0];self.assertEqual(gene['signed_weight'],-.1)
        self.assertEqual(gene['expression_relation'],'lower_expression_signature')
        self.assertNotIn('frequency',gene);self.assertNotIn('observed_samples',gene)
        self.assertEqual(self.repo.get_subtype_genes('CMS1',1,'missing')['matching_genes'],0)

    def test_missing_data_does_not_fall_back_to_ecograph(self):
        absent=CRCSCRepository(self.root/'absent')
        self.assertEqual(query(absent,'get_drug_predictions',subtype='CMS1')['status'],'backend_unavailable')
        self.assertEqual(query(absent,'get_model_evaluation')['status'],'backend_unavailable')

    async def test_tools_bound_per_agent_and_sdk_uses_nebius(self):
        original={t.name:t for t in tools_for_backend()}
        requests=[]
        def respond(request):
            body=json.loads(request.content);requests.append((str(request.url),body))
            calls=[('get_subtypes',{'disease':'colorectal_cancer'}),('get_drug_predictions',{'subtype':'CMS1','top_k':1}),('get_subtype_genes',{'subtype':'CMS1','limit':1,'gene':''})] if len(requests)==1 else []
            message={'role':'assistant','content':None if calls else 'Cached CMS expression hypotheses, not efficacy.'}
            if calls:message['tool_calls']=[{'id':f'c{i}','type':'function','function':{'name':name,'arguments':json.dumps(arguments)}} for i,(name,arguments) in enumerate(calls)]
            return httpx2.Response(200,json={'id':f'chat{len(requests)}','object':'chat.completion','created':1,'model':'fixture','choices':[{'index':0,'message':message,'finish_reason':'tool_calls' if calls else 'stop'}],'usage':{'prompt_tokens':10,'completion_tokens':10,'total_tokens':20}})
        settings=Settings(api_key='dummy',model='fixture')
        with patch('agent.backends.make_repository',return_value=self.repo):
            async with create_client(settings,http_client=httpx2.AsyncClient(transport=httpx2.MockTransport(respond))) as client:
                model=create_agent(settings,client,backend='crcsc')
                result=await run_turn(model,[],'Investigate CMS1',DebugHooks('dummy',lambda _:None))
        self.assertIn('Cached CMS',result.final_output)
        self.assertTrue(all(url==NEBIUS_BASE_URL+'chat/completions' for url,_ in requests))
        outputs=[ast.literal_eval(m['content']) for m in requests[-1][1]['messages'] if m['role']=='tool']
        self.assertTrue(any(o.get('profile_type')=='expression_signature' for o in outputs))
        self.assertTrue(all(o.get('data_backend')=='crcsc' for o in outputs))
        context=ToolContext(context=None,tool_name='get_subtypes',tool_call_id='legacy',tool_arguments='{}',run_config=RunConfig(tracing_disabled=True))
        legacy=await original['get_subtypes'].on_invoke_tool(context,'{"disease":"colorectal_cancer"}')
        self.assertIn('GI.CIN',legacy['subtypes']);self.assertNotIn('CMS1',legacy['subtypes'])


class BackendHTTPTests(unittest.TestCase):
    def test_catalog_bootstrap_and_recording_boundaries(self):
        with TestClient(create_app()) as client:
            old=client.get('/api/bootstrap').json();new=client.get('/api/bootstrap?backend=crcsc').json()
            self.assertEqual(old['backend'],'ecograph');self.assertTrue(old['recorded_available'])
            self.assertEqual(new['backend'],'crcsc');self.assertTrue(new['recorded_available'])
            self.assertEqual(new['selected_backend']['subtypes'],['CMS1','CMS2','CMS3','CMS4'])
            self.assertIn('ppi_pathway',new['evaluation']['conditions'])
            self.assertNotIn('target_only',new['evaluation']['conditions'])
            self.assertEqual(client.get('/api/bootstrap?backend=bad').status_code,400)
            self.assertEqual(client.post('/api/chat',json={'session_id':'one','mode':'recorded','backend':'crcsc','scenario_id':'candidate_audit'}).status_code,200)

    def test_conversation_cannot_change_backend_without_reset(self):
        calls=[]
        async def answer(prompt,history,emit,backend='ecograph'):
            calls.append((backend,list(history)))
            return backend,history+[{'role':'user','content':prompt}]
        with patch('ui.server.configuration',return_value={'live_available':True,'model':'fixture'}),TestClient(create_app(answer)) as client:
            def turn(backend):return client.post('/api/chat',json={'session_id':'same','mode':'live','backend':backend,'prompt':'Hello'})
            self.assertEqual(turn('ecograph').status_code,200)
            self.assertEqual(turn('crcsc').status_code,409)
            client.post('/api/reset',json={'session_id':'same'})
            self.assertEqual(turn('crcsc').status_code,200)
            self.assertEqual(calls,[('ecograph',[]),('crcsc',[])])

if __name__=='__main__':unittest.main()
