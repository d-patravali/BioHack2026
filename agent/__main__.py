"""Interactive or single-prompt CLI: python -m agent."""

import argparse
import asyncio
import json
import sys

from agents import AgentsException
from openai import APIConnectionError, APIStatusError, APITimeoutError, AuthenticationError, RateLimitError

from .config import NEBIUS_BASE_URL, load_settings
from .runtime import DebugHooks, create_agent, create_client, run_turn
from .repository import query
from .backends import BACKENDS, make_repository, catalog


def error_message(error):
    # Provider exception bodies may echo submitted content; don't print them.
    if isinstance(error, AuthenticationError):
        return "Nebius authentication failed. Check NEBIUS_API_KEY."
    if isinstance(error, RateLimitError):
        return "Nebius rate limit reached. Try again later."
    if isinstance(error, APITimeoutError):
        return "Nebius request timed out. Try again."
    if isinstance(error, APIConnectionError):
        return "Could not connect to Nebius. Check network access."
    if isinstance(error, APIStatusError):
        return f"Nebius returned HTTP {error.status_code}. Check the configured model and account access."
    if isinstance(error, AgentsException):
        return "Agent execution failed or reached its turn limit. Try a simpler prompt."
    return "Unexpected agent error. Run the offline tests and check installed dependency versions."


async def main_async(args):
    if args.check_data:
        repository = make_repository(args.backend)
        subtype = catalog(args.backend)["subtypes"][0]
        checks = {
            "predictions": query(repository, "get_drug_predictions", subtype=subtype, top_k=5),
            "targets": query(repository, "get_drug_targets", drug_id="DB00596"),
            "subtype_genes": query(repository, "get_subtype_genes", subtype=subtype, limit=1),
            "phenotype_evidence": query(repository, "get_gene_phenotypes", gene="2908", limit=1),
            "evaluation": query(repository, "get_model_evaluation"),
        }
        print(json.dumps({"mode": "offline_local_data_check", "data_backend": args.backend, "fresh_inference": False,
            "checks": {name: value["status"] for name, value in checks.items()},
            "returned_tied_candidates": checks["predictions"].get("returned_count"),
            "registered_backend": "reviewed cached experiment and local graph evidence"}, indent=2))
        return 0 if all(value["status"] == "ok" for value in checks.values()) else 2
    try:
        settings = load_settings()
    except ValueError as error:
        print(f"Configuration error: {error}", file=sys.stderr)
        return 2
    async with create_client(settings) as client:
        agent = create_agent(settings, client, backend=args.backend)
        hooks = DebugHooks(secret=settings.api_key)
        print(f"Provider: Nebius | Model: {hooks.safe(settings.model)}")
        print(f"API: Chat Completions | Endpoint: {NEBIUS_BASE_URL}chat/completions")
        print("OpenAI tracing: disabled")
        print("Data backend: " + catalog(args.backend)["label"])
        print("Tools: " + ", ".join(tool.name for tool in agent.tools))
        if args.check:
            print("Configuration and agent startup OK (no network request made).")
            return 0
        history = []
        if args.prompt is None:
            print("Enter a prompt. /reset clears history; /quit exits. History stays in memory.")
        while True:
            try:
                prompt = args.prompt if args.prompt is not None else input("\nYou> ").strip()
            except EOFError:
                return 0
            if not prompt:
                if args.prompt is not None:
                    print("Prompt must not be empty.", file=sys.stderr)
                    return 2
                continue
            if prompt in {"/quit", "/exit"}:
                return 0
            if prompt == "/reset" and args.prompt is None:
                history = []
                print("History cleared.")
                continue
            try:
                result = await run_turn(agent, history, prompt, hooks)
                print("\nAgent> " + hooks.safe(result.final_output))
                history = result.to_input_list()
            except Exception as error:
                print(error_message(error), file=sys.stderr)
                if args.prompt is not None:
                    return 1
            if args.prompt is not None:
                return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--backend", choices=BACKENDS, default="ecograph", help="Fixed data/model backend; original ecograph remains default.")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--prompt", help="Run one prompt and exit (contacts Nebius).")
    mode.add_argument("--check", action="store_true", help="Validate configuration and startup offline.")
    mode.add_argument("--check-data", action="store_true", help="Check local graph/result backends offline; no credentials or model calls.")
    try:
        return asyncio.run(main_async(parser.parse_args()))
    except KeyboardInterrupt:
        print("\nGoodbye.")
        return 0


if __name__ == "__main__":
    sys.exit(main())
