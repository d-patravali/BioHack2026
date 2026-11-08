"""Offline SDK integration tests; no API keys or live service required."""

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import httpx2
from agents import RunConfig
from agents.tool_context import ToolContext

from agent.config import NEBIUS_BASE_URL, Settings, load_settings
from agent.runtime import DebugHooks, create_agent, create_client, run_turn
from agent.tools import SUBTYPES, get_subtypes


class ConfigurationTests(unittest.TestCase):
    def test_missing_configuration_is_clear(self):
        with patch.dict("os.environ", {}, clear=True):
            with self.assertRaisesRegex(ValueError, "NEBIUS_API_KEY, NEBIUS_MODEL"):
                load_settings(Path("/tmp/biohack-nonexistent-env"))

    def test_dotenv_and_environment_precedence(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / ".env"
            path.write_text("NEBIUS_API_KEY=dummy-secret\nNEBIUS_MODEL=file-model\n")
            with patch.dict("os.environ", {"NEBIUS_MODEL": "env-model"}, clear=True):
                settings = load_settings(path)
        self.assertEqual(settings.model, "env-model")
        self.assertEqual(settings.api_key, "dummy-secret")
        self.assertNotIn("dummy-secret", repr(settings))


class ToolAndRunnerTests(unittest.IsolatedAsyncioTestCase):
    async def test_supported_and_unsupported_disease(self):
        context = ToolContext(context=None, tool_name="get_subtypes", tool_call_id="call_unit",
                              tool_arguments="{}", run_config=RunConfig(tracing_disabled=True))
        result = await get_subtypes.on_invoke_tool(context, '{"disease":"Colorectal Cancer"}')
        self.assertEqual(result["subtypes"], list(SUBTYPES))
        unsupported = await get_subtypes.on_invoke_tool(context, '{"disease":"unknown"}')
        self.assertEqual(unsupported["status"], "unsupported_disease")
        self.assertEqual(unsupported["subtypes"], [])

    async def test_malformed_arguments_return_tool_error(self):
        context = ToolContext(context=None, tool_name="get_subtypes", tool_call_id="call_bad",
                              tool_arguments='{"unexpected":1}',
                              run_config=RunConfig(tracing_disabled=True, trace_include_sensitive_data=False))
        result = await get_subtypes.on_invoke_tool(context, '{"unexpected":1}')
        self.assertIsInstance(result, str)
        self.assertIn("error", result.lower())

    async def test_real_sdk_loop_uses_only_nebius_chat_completions(self):
        requests = []
        def respond(request):
            body = json.loads(request.content)
            requests.append((str(request.url), body))
            if len(requests) == 1:
                message = {"role": "assistant", "content": None, "tool_calls": [{
                    "id": "call_test", "type": "function", "function": {
                        "name": "get_subtypes", "arguments": '{"disease":"colorectal_cancer"}'}}]}
                reason = "tool_calls"
            else:
                message = {"role": "assistant", "content": ", ".join(SUBTYPES)}
                reason = "stop"
            return httpx2.Response(200, json={"id": f"chat_test_{len(requests)}",
                "object": "chat.completion", "created": 1, "model": "test-nebius-model",
                "choices": [{"index": 0, "message": message, "finish_reason": reason}],
                "usage": {"prompt_tokens": 20, "completion_tokens": 10, "total_tokens": 30}})

        logs = []
        settings = Settings(api_key="dummy-secret", model="test-nebius-model")
        transport = httpx2.MockTransport(respond)
        async with create_client(settings, http_client=httpx2.AsyncClient(transport=transport)) as client:
            agent = create_agent(settings, client)
            with patch.dict("os.environ", {"OPENAI_API_KEY": "", "OPENAI_BASE_URL": "https://invalid.example"}):
                result = await run_turn(agent, [], "List colorectal cancer subtypes", DebugHooks(settings.api_key, logs.append))
        self.assertEqual(len(requests), 2)
        for url, body in requests:
            self.assertEqual(url, NEBIUS_BASE_URL + "chat/completions")
            self.assertEqual(body["model"], settings.model)
        tool_messages = [m for m in requests[1][1]["messages"] if m["role"] == "tool"]
        self.assertEqual(len(tool_messages), 1)
        self.assertIn("GI.HM-indel", tool_messages[0]["content"])
        self.assertIn("GI.GS", result.final_output)
        self.assertTrue(any("[tool call] get_subtypes" in line for line in logs))
        self.assertTrue(any("[tool result] get_subtypes" in line for line in logs))
        self.assertGreater(len(result.to_input_list()), 1)

    async def test_debug_redacts_secret(self):
        logs = []
        hooks = DebugHooks("dummy-secret", logs.append)
        await hooks.on_tool_end(None, None, get_subtypes, "dummy-secret")
        self.assertNotIn("dummy-secret", "".join(logs))
        self.assertIn("[REDACTED]", "".join(logs))


if __name__ == "__main__":
    unittest.main()
