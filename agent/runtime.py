"""One SDK agent with an explicit Nebius Chat Completions model."""

import json

from agents import Agent, OpenAIChatCompletionsModel, RunConfig, RunHooks, Runner, set_tracing_disabled
from openai import AsyncOpenAI

from .config import NEBIUS_BASE_URL, Settings
from .tools import get_subtypes


class DebugHooks(RunHooks):
    def __init__(self, secret: str = "", emit=print):
        self.secret = secret
        self.emit = emit

    def safe(self, value) -> str:
        text = value if isinstance(value, str) else json.dumps(value, default=str)
        return text.replace(self.secret, "[REDACTED]") if self.secret else text

    async def on_tool_start(self, context, agent, tool):
        self.emit(f"[tool call] {tool.name} {self.safe(getattr(context, 'tool_arguments', ''))}")

    async def on_tool_end(self, context, agent, tool, result):
        self.emit(f"[tool result] {tool.name} {self.safe(result)}")


def create_client(settings: Settings, **kwargs) -> AsyncOpenAI:
    return AsyncOpenAI(api_key=settings.api_key, base_url=NEBIUS_BASE_URL,
                       timeout=60., max_retries=1, **kwargs)


def create_agent(settings: Settings, client: AsyncOpenAI) -> Agent:
    # Disable the SDK's separate OpenAI trace exporter as well as per-run tracing.
    set_tracing_disabled(True)
    if str(client.base_url) != NEBIUS_BASE_URL:
        raise ValueError("This application requires the fixed Nebius Token Factory endpoint.")
    return Agent(
        name="Biomedical research agent",
        instructions=(
            "You are a biomedical knowledge-graph research assistant. "
            "Use get_subtypes to retrieve available disease subtypes; do not invent the catalog. "
            "Quote subtype identifiers exactly as returned. Do not expand abbreviations or add "
            "biological descriptions that are absent from the tool result. "
            "Map natural-language 'colorectal cancer' to colorectal_cancer. "
            "Report unsupported diseases honestly. Only subtype listing is implemented. "
            "GNN predictions, prediction comparisons, and molecular graph queries are not connected. "
            "Do not invent drug rankings, tool results, efficacy claims, or clinical recommendations."
        ),
        model=OpenAIChatCompletionsModel(model=settings.model, openai_client=client),
        tools=[get_subtypes],
    )


async def run_turn(agent: Agent, history: list, prompt: str, hooks: DebugHooks):
    """Return updated history only after a successful, bounded SDK run."""
    return await Runner.run(agent, input=[*history, {"role": "user", "content": prompt}],
                            hooks=hooks, max_turns=6,
                            run_config=RunConfig(tracing_disabled=True))
