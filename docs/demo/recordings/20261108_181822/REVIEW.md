# Presenter review: CRCSC live recording

Two scenarios were executed against Nebius GLM-5.3-Flash with the CRCSC backend.
Tool-evidence checks passed for both; all eight input-artifact hashes were unchanged.
The first scenario made five tool calls (including two target lookups); the second
made two. Model answers and tool events are preserved exactly as recorded.

## Verified findings

- CMS2 cached seed42 signed-profile predictions return a 29-drug rank1 tie.
- Recorded target IDs and two explicit target lookups identify NR3C1/NCBI2908.
- The model's sampled indication AP is below popularity, with no CMS therapeutic
  labels and no efficacy/novelty validation.
- CMS1's leading absolute-weight signature genes are HPSE, GNLY and SECTM1, all
  higher-expression associations in the displayed subset. These are expression
  weights/FDR, not mutation frequencies or gene-level patient mutation counts.
- HPSE has no recorded direct phenotype association in this local snapshot;
  this is missing knowledge, not a biological negative.

## Prose qualifications

The first answer calls the entire group glucocorticoids/structurally related
corticosteroids. **Drug class and structural similarity were not supplied by the
tools and have not been independently verified for this recording.** The verified
common relationship is the recorded NR3C1 target. The replay displays this caveat.
Its "CRC2/CMS2" wording is a typo; the tool queried CMS2, not another subtype.

The second answer's "PPI/pathway graph only" wording should mean phenotype edges
were excluded. The drug encoder also uses target and disease–gene relationships.
The replay clarifies that scope; retrieved phenotype context is not a score
attribution or causal explanation.

These are saved live LLM answers over cached results, not fresh GNN inference,
a literature review, or a systematic demonstration of agent truthfulness.
