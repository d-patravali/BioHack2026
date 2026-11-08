"""Actual saved CRCSC replay is offline, backend-specific, and unaltered."""
import json
import unittest
from unittest.mock import patch
from starlette.testclient import TestClient
from ui.server import create_app, read_recording, RECORDING, CMS_RECORDING
from scripts.rehearse_agent import verify_tool_evidence


def events(response):
    return [json.loads(line[6:]) for line in response.text.splitlines() if line.startswith('data: ')]


class CRCSCRecordingTests(unittest.TestCase):
    def test_live_record_metadata_and_tool_checks(self):
        record=read_recording('crcsc')
        self.assertEqual(record['data_backend'],'crcsc')
        self.assertEqual(record['provider'],'Nebius Token Factory')
        self.assertTrue(record['source_hashes_unchanged'])
        self.assertTrue(record['all_tool_checks_passed'])
        self.assertEqual(len(record['scenarios']),2)
        for scenario in record['scenarios']:
            self.assertTrue(all(verify_tool_evidence(scenario['id'],scenario['events'],'crcsc').values()))

    def test_both_crcsc_replays_need_no_credentials_or_provider(self):
        async def forbidden(*args,**kwargs):
            self.fail('Recorded mode invoked the live provider')
        record=read_recording('crcsc')
        with patch('ui.server.configuration',return_value={'live_available':False,'model':'Not configured'}),TestClient(create_app(forbidden)) as client:
            self.assertTrue(client.get('/api/bootstrap?backend=crcsc').json()['recorded_available'])
            for scenario in record['scenarios']:
                response=client.post('/api/chat',json={'session_id':'offline-cms','backend':'crcsc','mode':'recorded','scenario_id':scenario['id']})
                self.assertEqual(response.status_code,200)
                trace=events(response)
                self.assertEqual(trace[0]['backend'],'crcsc')
                answer=next(e for e in trace if e['type']=='answer')
                self.assertEqual(answer['text'],scenario['answer'])
                self.assertEqual(answer['backend'],'crcsc')
                if scenario['id']=='candidate_audit':self.assertIn('unverified',answer['review_note'])
                self.assertEqual(sum(e['type']=='tool_call' for e in trace),sum(e['event']=='tool_call' for e in scenario['events']))
                self.assertEqual(trace[-1]['type'],'done')

    def test_missing_or_wrong_backend_recording_never_falls_back(self):
        with patch('ui.server.CMS_RECORDING',RECORDING):
            self.assertIsNone(read_recording('crcsc'))
            with TestClient(create_app()) as client:
                self.assertFalse(client.get('/api/bootstrap?backend=crcsc').json()['recorded_available'])
                response=client.post('/api/chat',json={'session_id':'mismatch','backend':'crcsc','mode':'recorded','scenario_id':'candidate_audit'})
                self.assertEqual(response.status_code,400)
        self.assertTrue(read_recording('ecograph')['all_tool_checks_passed'])

if __name__=='__main__':unittest.main()
