"""Record two live agent walkthroughs, or replay their saved transcript offline.

Examples:
  python scripts/rehearse_agent.py --live
  python scripts/rehearse_agent.py --replay docs/demo/recordings/RUN/transcript.json
"""

import argparse
import asyncio
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

SCENARIOS = [
    {
        "id": "candidate_audit",
        "title": "Investigate the tied GI.CIN candidate group",
        "prompt": (
            "Investigate the cached leading candidates for GI.CIN. Retrieve the top five, "
            "inspect a shared drug target, check its recorded phenotype associations and "
            "whether it is in the retained GI.CIN positive-mutation profile, and inspect "
            "the popularity and molecular-context ablation evaluation. Use the tools. "
            "Summarize what we can actually conclude in at most 180 words, with source paths."
        ),
    },
    {
        "id": "phenotype_gene_evidence",
        "title": "Inspect APC phenotype–genotype evidence",
        "prompt": (
            "Now retrieve the three most frequent observed genes for GI.CIN and look up "
            "all recorded direct phenotype associations for APC (NCBI324). Use tools. "
            "Report APC's observed/profiled patient counts, separate cohort observations "
            "from graph associations, and state whether the CRC GNN used these phenotype "
            "edges. Do not imply that these annotations explain the drug ranking. "
            "Summarize in at most 180 words, with source paths."
        ),
    },
]


CMS_SCENARIOS = [
    {'id': 'candidate_audit', 'title': 'Investigate the cached CMS2 drug ranking',
     'short_title': 'Why do these drugs tie?', 'description': 'Inspect CMS2 targets and model controls',
     'prompt': 'Which drugs does the model rank highest for CMS2, and how much confidence should I place in that ranking? Check the recorded targets and model evaluation. Keep the answer brief and cite tool source paths.'},
    {'id': 'phenotype_gene_evidence', 'title': 'Inspect CMS1 expression and phenotype context',
     'short_title': 'Inspect CMS expression', 'description': 'Separate expression from mutation evidence',
     'prompt': 'Show the three highest-weight CMS1 expression signature genes, look up recorded phenotype associations for the first gene, and explain what these connections do and do not tell us. Keep the answer brief and cite tool source paths.'},
]


def artifact_hashes(backend='ecograph'):
    paths = [
        "data/processed/primekg/nodes.csv",
        "data/processed/primekg/phenotype_disease_gene_edges.csv",
        "reports/compatibility/primekg_direct_targets.csv",
        "reports/crc_experiments/20261108_034016/evaluation_metrics.csv",
        "reports/crc_experiments/20261108_034016/tie_aware_top_candidates.csv",
        "reports/crc_experiments/20261108_034016/profile_gene_evidence.csv.gz",
    ]
    if backend == 'crcsc':
        paths = paths[:3] + [
            'reports/crcsc_drug_gnn/20261108_173948/manifest.json',
            'reports/crcsc_drug_gnn/20261108_173948/evaluation_metrics.csv',
            'reports/crcsc_drug_gnn/20261108_173948/candidate_hypotheses.csv.gz',
            'reports/crcsc_drug_gnn/20261108_173948/cms_expression_edges.csv',
            'reports/crcsc/20261108_172638/cohort_counts.csv',
        ]
    hashes = {}
    for relative in paths:
        digest = hashlib.sha256()
        with (ROOT / relative).open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
        hashes[relative] = digest.hexdigest()
    return hashes


def verify_tool_evidence(scenario_id, events, backend='ecograph'):
    """Check tool execution/results, not the truthfulness of arbitrary LLM prose."""
    outputs = [e for e in events if e["event"] == "tool_result"]
    def values(name):
        return [e["result"] for e in outputs if e["tool"] == name and isinstance(e["result"], dict)]
    ok = {"all_returned_tools_ok": bool(outputs) and all(isinstance(e["result"], dict) and e["result"].get("status") == "ok" for e in outputs)}
    if backend == 'crcsc':
        predictions = values('get_drug_predictions')
        evaluations = values('get_model_evaluation') + [p['evaluation'] for p in predictions if 'evaluation' in p]
        if scenario_id == 'candidate_audit':
            targets = values('get_drug_targets')
            ok.update({
                'cms2_reviewed_cache_retrieved': any(p.get('subtype') == 'CMS2' and p.get('data_backend') == 'crcsc' and p.get('profile') == 'signed' for p in predictions),
                'complete_29_drug_tie_returned': any(p.get('returned_count') == 29 and p.get('fresh_inference') is False for p in predictions),
                'nr3c1_target_retrieved': any(any(t['gene_id'] == '2908' for t in p.get('targets', [])) for p in targets),
                'crcsc_baseline_and_ablation_available': any({'ppi_pathway', 'drug_only', 'drug_popularity'} <= p.get('conditions', {}).keys() and p.get('data_backend') == 'crcsc' for p in evaluations),
            })
        else:
            profiles = values('get_subtype_genes')
            genes = [g for p in profiles for g in p.get('genes', [])]
            ok.update({
                'cms1_expression_profile_retrieved': any(p.get('subtype') == 'CMS1' and p.get('data_backend') == 'crcsc' and p.get('profile_type') == 'expression_signature' for p in profiles),
                'three_expression_genes_retrieved': {'10855','10578','6398'} <= {g.get('gene_id') for g in genes},
                'expression_not_mutation_counts': bool(genes) and all('signed_weight' in g and 'discovery_q' in g and 'observed_samples' not in g and 'frequency' not in g for g in genes),
                'hpse_no_recorded_direct_phenotypes': any(p.get('gene_id') == '10855' and p.get('phenotype_count') == 0 and p.get('model_uses_phenotype_edges') is False for p in values('get_gene_phenotypes')),
            })
        return ok
    if scenario_id == "candidate_audit":
        predictions = values("get_drug_predictions")
        targets = values("get_drug_targets")
        phenotypes = values("get_gene_phenotypes")
        profiles = values("get_subtype_genes")
        evaluations = values("get_model_evaluation") + [p["evaluation"] for p in predictions if "evaluation" in p]
        ok.update({
            "cached_prediction_tool_executed": bool(predictions),
            "complete_29_drug_tie_returned": any(p.get("returned_count") == 29 and p.get("fresh_inference") is False for p in predictions),
            "nr3c1_target_retrieved": any(any(t["gene_id"] == "2908" for t in p.get("targets", [])) for p in targets),
            "nr3c1_missing_direct_phenotypes_retrieved": any(p.get("gene_id") == "2908" and p.get("phenotype_count") == 0 for p in phenotypes),
            "nr3c1_retained_profile_checked": any(p.get("gene_filter") in {"2908", "NR3C1"} and p.get("matching_genes") == 0 for p in profiles),
            "baseline_and_ablation_available": any({"drug_popularity", "drug_only", "target_only"} <= p.get("conditions", {}).keys() for p in evaluations),
        })
    else:
        genes = [g for p in values("get_subtype_genes") for g in p.get("genes", [])]
        phenotypes = values("get_gene_phenotypes")
        ok.update({
            "subtype_gene_tool_executed": bool(values("get_subtype_genes")),
            "apc_counts_retrieved": any(g["gene_id"] == "324" and g["observed_samples"] == 185 and g["profiled_samples"] == 224 for g in genes),
            "apc_four_direct_phenotypes_retrieved": any(p.get("gene_id") == "324" and p.get("phenotype_count") == 4 and p.get("returned_count") == 4 for p in phenotypes),
            "phenotype_model_boundary_available": any(p.get("model_uses_phenotype_edges") is False for p in phenotypes),
        })
    return ok


def render(record):
    lines = ["# Recorded live agent rehearsal", "",
        f"Recorded UTC: {record['recorded_utc']}", f"Provider/model: {record['provider']} / {record['model']}",
        f"Data backend: {record.get('data_backend', 'ecograph')}",
        "", "This is a saved live tool-calling session, not a fresh model inference run or a systematic agent evaluation.",
        "The GNN scores are cached. Replay requires no network or credentials.", ""]
    for scenario in record["scenarios"]:
        lines.extend([f"## {scenario['title']}", "", "### User prompt", "", scenario["prompt"], "", "### Tool trace", ""])
        for event in scenario["events"]:
            if event["event"] == "tool_call":
                lines.append(f"- `{event['tool']}`: `{json.dumps(event['arguments'], ensure_ascii=False)}`")
        lines.extend(["", "### Actual model response", "", scenario.get("answer", "No final answer returned."),
                      "", "### Tool-evidence checks", ""])
        lines.extend(f"- {name}: **{'PASS' if passed else 'FAIL'}**" for name, passed in scenario["tool_checks"].items())
        if "error" in scenario:
            lines.extend(["", scenario["error"]])
        lines.append("")
    lines.extend(["## Scope", "", "The checks verify tool results and required calls, not all semantic claims in the response. "
        "Inspect the actual answers before presenting. No clinical efficacy, subtype response accuracy, or literature novelty was established.", ""])
    return "\n".join(lines)


async def record_live(output, backend='ecograph'):
    from agent.config import load_settings, NEBIUS_BASE_URL
    from agent.runtime import DebugHooks, create_agent, create_client, run_turn
    from agent.__main__ import error_message

    settings = load_settings()
    before = artifact_hashes(backend)
    output.mkdir(parents=True, exist_ok=False)
    class RecordingHooks(DebugHooks):
        def __init__(self):
            super().__init__(settings.api_key)
            self.events = []
        def cleaned(self, value):
            try:
                return json.loads(self.safe(value))
            except (json.JSONDecodeError, TypeError):
                return self.safe(value)
        async def on_tool_start(self, context, agent, tool):
            self.events.append({"event": "tool_call", "tool": tool.name,
                "call_id": getattr(context, "tool_call_id", None),
                "arguments": self.cleaned(getattr(context, "tool_arguments", "{}"))})
            print(f"  tool: {tool.name}", flush=True)
        async def on_tool_end(self, context, agent, tool, result):
            self.events.append({"event": "tool_result", "tool": tool.name,
                "call_id": getattr(context, "tool_call_id", None), "result": self.cleaned(result)})

    record = {"recorded_utc": datetime.now(timezone.utc).isoformat(), "provider": "Nebius Token Factory",
        "data_backend": backend, "model": settings.model, "endpoint": NEBIUS_BASE_URL+"chat/completions",
        "gnn_backend": "cached experiment, not fresh GNN inference", "scenarios": [], "source_sha256": before}
    def save():
        # A final pass also redacts a key accidentally echoed by the remote model.
        serialized = json.dumps(record, indent=2, ensure_ascii=False).replace(settings.api_key, "[REDACTED]")
        (output/"transcript.json").write_text(serialized+"\n")
        safe_record = json.loads(serialized)
        (output/"walkthrough.md").write_text(render(safe_record))
    async with create_client(settings) as client:
        agent = create_agent(settings, client, backend=backend)
        history = []
        for scenario in (CMS_SCENARIOS if backend == 'crcsc' else SCENARIOS):
            print(scenario["title"], flush=True)
            hooks = RecordingHooks()
            entry = {**scenario, "events": hooks.events, "tool_checks": {}}
            record["scenarios"].append(entry)
            try:
                result = await run_turn(agent, history, scenario["prompt"], hooks)
                entry["answer"] = hooks.safe(result.final_output)
                history = result.to_input_list()
                print(entry["answer"], flush=True)
            except Exception as error:
                entry["error"] = error_message(error)
                print(entry["error"], file=sys.stderr)
            entry["tool_checks"] = verify_tool_evidence(scenario["id"], hooks.events, backend)
            save()
            if "error" in entry:
                break
    record["source_hashes_unchanged"] = before == artifact_hashes(backend)
    record["all_tool_checks_passed"] = len(record["scenarios"]) == 2 and all(
        all(s["tool_checks"].values()) and "error" not in s for s in record["scenarios"])
    save()
    print(f"Saved: {output/'walkthrough.md'}", flush=True)
    return 0 if record["all_tool_checks_passed"] and record["source_hashes_unchanged"] else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--live", action="store_true", help="Contacts Nebius and may incur API charges.")
    mode.add_argument("--replay", type=Path, help="Read a saved JSON transcript offline; makes no API calls.")
    parser.add_argument("--backend", choices=["ecograph", "crcsc"], default="ecograph", help="Select the backend to record; replay uses its saved metadata.")
    parser.add_argument("--output", type=Path, help="New recording directory, never overwrite an existing run.")
    args = parser.parse_args()
    if args.replay:
        print(render(json.loads(args.replay.read_text())))
        return 0
    output = args.output or ROOT/"docs/demo/recordings"/datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
    if output.exists():
        parser.error("Recording directory must be new")
    return asyncio.run(record_live(output, args.backend))


if __name__ == "__main__":
    sys.exit(main())
