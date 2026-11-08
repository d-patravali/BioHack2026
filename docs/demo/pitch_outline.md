# Pitch outline: an evidence-audited biomedical research agent

## One-sentence description

We built a tool-calling research agent that connects recorded phenotype–gene
evidence and colorectal subtype observations with GNN-generated hypotheses,
then inspects baseline and ablation results before presenting conclusions.

## Suggested 90-second story

**Problem (15 seconds).** Biomedical graph models can produce convincing drug
rankings even when the model barely uses disease-specific information. A
researcher needs to inspect where an output comes from and whether its evidence
supports the proposed interpretation.

**Prototype (20 seconds).** Our framework separates graph integration, GNN
experiments and a Nebius-backed agent. Six local tools expose subtype catalogs,
cached model scores, direct drug targets, cohort gene frequencies, explicit
phenotype–gene annotations and model evaluation. Tool calls and source paths
are visible. The agent does not silently turn a graph connection into a
therapeutic claim.

**Demonstration (35 seconds).** Ask for five GI.CIN candidates. The tool returns
29 tied compounds, rather than pretending an arbitrary five are uniquely best.
The agent checks their shared NR3C1 target, missing recorded direct phenotype
annotations, and the retained subtype observation profile. It also retrieves
the popularity and molecular-context ablation comparisons. A second query
shows recorded APC phenotype associations alongside its 185/224 GI.CIN
observation count, explicitly separating the two kinds of evidence.

**Finding and next step (20 seconds).** On our fixed sampled indication-recovery
test, drug popularity AP is 0.8941 and target-only GNN mean AP is 0.8254. A
drug-only ablation reaches 0.8257. This prototype exposes the limitation rather
than marketing the output as a discovered treatment. The next research step is
to evaluate stronger phenotype–genotype representations and independent
subtype evidence, while retaining this audit trail.

## Three slides

1. **Question and architecture:** source KG + ECoGraph observations -> separate
   graph/model experiments -> read-only tools -> research agent. Label the GNN
   result backend **cached**, and phenotype retrieval as a separate evidence path.
2. **Live investigation:** a tied candidate group -> recorded targets -> gene
   observations/phenotypes -> baseline and ablation -> qualified conclusion.
3. **Measured findings and limits:** performance table, identical target
   neighborhoods, distinction between observed associations and efficacy, and
   the next validation milestone.

## Important boundaries and likely questions

- **Does the GNN currently learn from phenotype edges?** No. It uses a molecular
  subset of the KG; the agent separately retrieves phenotype–gene annotations.
  Do not imply that separate retrieval establishes phenotype-aware GNN learning.
- **Are subtype mutation counts also excluded?** No. Their frequencies are inputs
  to the frequency-weighted subtype profile. Counts being model inputs does not
  establish their causal importance for a particular drug ranking.
- **Is the GNN better than the baseline?** Not in this experiment. AP is a
  ranking metric on sampled known/unknown pairs, not treatment accuracy.
- **Are these new drugs for CRC?** No therapeutic discovery or literature
  novelty is established. They are unvalidated graph/model hypotheses.
- **Does the agent execute ablations live?** No. It retrieves the actual saved
  ablation results. Neither the model nor graph is retrained during the demo.
- **Is the agent's benefit systematically evaluated?** No. A successful
  walkthrough demonstrates orchestration, not reduced hallucination rates or
  superior research decisions. Those need a separate evaluation.
- **What is original?** Present the implementation and evaluated case study as
  a hackathon engineering contribution. Do not claim a new GNN algorithm or a
  first-of-its-kind biomedical agent.

## Demonstration procedure

1. Activate `biohack-agent` and run `python -m agent --check-data`.
2. Launch `python -m agent`; use the two prompts in `agent/DEMO.md`.
3. If the network fails, open the recorded `walkthrough.md` or use offline replay
   and explicitly identify it as a **previously recorded live session**.
4. Keep the source report available:
   `reports/crc_experiments/20261108_034016/crc_findings.md`.

Recorder: `python scripts/rehearse_agent.py --live` (Nebius charges may apply).
Replay: `python scripts/rehearse_agent.py --replay PATH/TO/transcript.json`.
