# BioScouter Manuscript Reproducibility Package

This folder freezes the benchmark evidence used by the v3 concise manuscript.
It is intended to be archived with a DOI before journal submission.

## Files

- `benchmark-queries.csv` - legacy ten-query snapshot retained for audit history; not the primary controlled validation panel.
- `benchmark-freeze.json` - legacy frozen snapshot retained so older manuscript numbers can be traced and removed or replaced.
- `controlled-benchmark-queries.csv` - primary 30-query held-out evaluation panel.
- `corpus-construction-queries.csv` - independent corpus-building queries used before held-out evaluation.
- `controlled_evaluation.py` - leakage-resistant controlled evaluation runner and annotation-pool generator.
- `relevance-rubric.md` - manual relevance-labeling rubric used for P@10, nDCG@10, MRR, agreement, and adjudication.
- `rerun_benchmark.py` - helper script for rerunning the query set against a local or deployed BioScouter API.
- `expanded-benchmark-queries.csv` - 30-query expanded coverage panel spanning transcriptomics, proteomics, metabolomics, epigenomics, genomics, single-cell, metagenomics, and multi-omics.
- `run_expanded_benchmark.py` - helper script that records a dated live snapshot for the expanded panel and creates top-10 labeling sheets.
- `expanded-benchmark-*/expanded-benchmark-freeze.json` - dated live expanded benchmark output.
- `expanded-benchmark-*/annotator_A_top10.csv` and `annotator_B_top10.csv` - independent relevance-labeling sheets.
- `score_relevance_labels.py` - computes P@10, nDCG@10, percent agreement, Cohen's kappa, and an adjudication sheet after annotator labels are filled.
- `independent-labeling-protocol.md` - protocol for the two-annotator relevance check.
- `run_source_reliability.py` - cold/warm live-source reliability runner.
- `run_concept_ablation.py` - paired runner comparing candidate sets with and without deterministic concept expansion.
- `usability-ethics-protocol.md` - optional ethics-ready usability protocol; no usability results are claimed until participant data exist.

## Important Interpretation Note

The manuscript values are a frozen snapshot, not a promise that live public APIs
will return identical counts forever. GEO, PRIDE, ENCODE, ENA, GDC, MassIVE,
CPTAC, and other sources update continuously and sometimes change API behavior.
Live reruns should be treated as reproducibility checks against this snapshot.

The expanded benchmark reports capped candidate retrieval (`max_results=100`)
and should not be interpreted as a complete count of all matching records in
the source repositories.

## Rerun Example

Start BioScouter, then run:

```powershell
python BioScouter_Paper\v3-concise\reproducibility\rerun_benchmark.py --base-url http://127.0.0.1:8001
```

If authentication is required:

```powershell
python BioScouter_Paper\v3-concise\reproducibility\rerun_benchmark.py --base-url http://127.0.0.1:8001 --token "<AUTH_TOKEN>"
```

The script writes a timestamped JSON rerun file beside this README. Compare that
file to `benchmark-freeze.json` before updating manuscript claims.

## Controlled Benchmark Example

Validate the controlled query/corpus split:

```powershell
python BioScouter_Paper\v3-concise\reproducibility\controlled_evaluation.py --validate-only
```

Run the controlled evaluation against a local or deployed API:

```powershell
python BioScouter_Paper\v3-concise\reproducibility\controlled_evaluation.py --base-url http://127.0.0.1:8001 --candidate-depth 100 --top-k 10
```

This writes raw responses, ranked outputs, run manifests, pool mappings, and
two blinded annotation sheets. Do not report final retrieval metrics until the
two independent annotator files are complete and adjudicated.

## Expanded Benchmark Example

```powershell
python BioScouter_Paper\v3-concise\reproducibility\run_expanded_benchmark.py --base-url http://127.0.0.1:8001 --max-results 100 --top-n 10
```

After two annotators complete the generated top-10 CSV files, run:

```powershell
python BioScouter_Paper\v3-concise\reproducibility\score_relevance_labels.py --annotator-a expanded-benchmark-YYYYMMDD-HHMMSS\annotator_A_top10.csv --annotator-b expanded-benchmark-YYYYMMDD-HHMMSS\annotator_B_top10.csv
```

Do not report independent-validation metrics until both annotator files are
completed and disagreements are adjudicated.

## Concept Expansion Ablation

The deterministic concept-normalization feature is evaluated as a paired
candidate-set ablation, not as a standalone relevance metric:

```powershell
python BioScouter_Paper\v3-concise\reproducibility\run_concept_ablation.py --base-url http://127.0.0.1:8001 --max-results 100
```

The output records result counts and candidate overlap for each query. Any
claim about improved relevance still requires the independent top-10 labels.
