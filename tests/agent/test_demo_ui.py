"""Offline HTTP checks for the local demo and its existing recorded evidence."""
import json
import unittest
from unittest.mock import patch
from starlette.testclient import TestClient
from ui.server import create_app


def events(response):
    return [json.loads(line[6:]) for line in response.text.splitlines() if line.startswith('data: ')]


class DemoUITests(unittest.TestCase):
    def test_recording_runs_without_provider(self):
        async def forbidden(*args, **kwargs):
            self.fail('Recorded replay must not invoke a provider')
        with patch('ui.server.configuration', return_value={'live_available': False, 'model': 'Not configured'}), TestClient(create_app(forbidden)) as client:
            self.assertEqual(client.get('/').status_code, 200)
            self.assertEqual(client.get('/static/app.js').status_code, 200)
            self.assertEqual(client.get('/static/.env').status_code, 404)
            info = client.get('/api/bootstrap').json()
            self.assertTrue(info['recorded_available'])
            self.assertFalse(info['live_available'])
            for scenario, calls in [('candidate_audit', 6), ('phenotype_gene_evidence', 2)]:
                response = client.post('/api/chat', json={'session_id': 'offline', 'mode': 'recorded', 'scenario_id': scenario})
                self.assertEqual(response.status_code, 200)
                trace = events(response)
                self.assertEqual(sum(e['type'] == 'tool_call' for e in trace), calls)
                self.assertTrue(any(e['type'] == 'answer' and e['mode'] == 'recorded' for e in trace))
                self.assertEqual(trace[-1]['type'], 'done')
            self.assertEqual(client.post('/api/chat', json={'session_id': 'offline', 'mode': 'live', 'prompt': 'Hello'}).status_code, 503)

    def test_live_session_history_reset_and_safe_errors(self):
        histories = []
        async def answer(prompt, history, emit, backend='ecograph'):
            histories.append(list(history))
            if prompt == 'fail':
                raise RuntimeError('secret-should-never-appear')
            await emit({'type': 'tool_call', 'tool': 'get_subtypes', 'call_id': 'one', 'arguments': {}})
            await emit({'type': 'tool_result', 'tool': 'get_subtypes', 'call_id': 'one', 'result': {'status': 'ok'}})
            return 'Four subtypes.', history + [{'role': 'user', 'content': prompt}]
        with patch('ui.server.configuration', return_value={'live_available': True, 'model': 'test'}), TestClient(create_app(answer)) as client:
            def turn(prompt):
                return client.post('/api/chat', json={'session_id': 'live-session', 'mode': 'live', 'prompt': prompt})
            self.assertEqual(events(turn('first'))[-2]['text'], 'Four subtypes.')
            turn('second')
            self.assertEqual(len(histories[1]), 1)
            failed = turn('fail')
            self.assertNotIn('secret-should-never-appear', failed.text)
            self.assertTrue(any(e['type'] == 'error' for e in events(failed)))
            turn('third')
            self.assertEqual(len(histories[-1]), 2)
            self.assertEqual(client.post('/api/reset', json={'session_id': 'live-session'}).status_code, 200)
            turn('after reset')
            self.assertEqual(histories[-1], [])

    def test_invalid_requests(self):
        with TestClient(create_app()) as client:
            for payload in [[], {}, {'session_id': []}, {'session_id': 'valid', 'mode': []}]:
                self.assertEqual(client.post('/api/chat', json=payload).status_code, 400)
            self.assertEqual(client.post('/api/reset', json={'session_id': []}).status_code, 400)
            self.assertEqual(client.post('/api/chat', json={}, headers={'Origin': 'https://other.example'}).status_code, 403)
            self.assertEqual(client.post('/api/chat', content='plain').status_code, 415)
            self.assertEqual(client.get('/', headers={'Host': 'other.example'}).status_code, 400)


if __name__ == '__main__':
    unittest.main()
