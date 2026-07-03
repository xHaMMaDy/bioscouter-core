# BioScouter Manuscript Reproducibility Package

This folder freezes the benchmark evidence used by the v3 concise manuscript.
The public scientific core and reproducibility snapshot are available at
https://github.com/xHaMMaDy/bioscouter-core and archived on Zenodo at
https://doi.org/10.5281/zenodo.21143417.

## Files

- `benchmark-queries.csv` - legacy ten-query snapshot retained for audit history; not the primary controlled validation panel.
- `benchmark-freeze.json` - legacy frozen snapshot retained so older manuscript numbers can be traced and removed or replaced.
- `controlled-benchmark-queries.csv` - earlier 30-query held-out panel retained for audit history.
- `benchmark-queries-80.csv` - broader 80-query controlled baseline and ranking-ablation panel across eight omics categories.
- `corpus-construction-queries.csv` - independent corpus-building queries used before held-out evaluation.
- `controlled_evaluation.py` - leakage-resistant controlled evaluation runner and annotation-pool generator.
- `controlled-evaluation-80-*/system-summary.csv` - per-system completion, returned top-10 count, unique-record count, observed-source count, and latency summary.
- `controlled-evaluation-80-*/benchmark_run_config.json` and `run_manifest.json` - exact run parameters, query/corpus hashes, and artifact checksums.
- `controlled-evaluation-80-*/baseline-comparison-summary.csv` - BioScouter hybrid, designated native-source, and OmicsDI comparison rows.
- `controlled-evaluation-80-*/ablation-summary.csv` - keyword-only, frozen-corpus embedding-only, and hybrid ranking conditions.
- `controlled-evaluation-80-*/normalized-top10-results.csv` - normalized ranked records used to audit system outputs.
- `relevance-rubric.md` - manual relevance-labeling rubric used for P@10, nDCG@10, MRR, agreement, and adjudication.
- `rerun_benchmark.py` - helper script for rerunning the query set against a local or deployed BioScouter API.
- `expanded-benchmark-queries.csv` - 30-query expanded coverage panel spanning transcriptomics, proteomics, metabolomics, epigenomics, genomics, single-cell, metagenomics, and multi-omics.
- `run_expanded_benchmark.py` - helper script that records a dated live snapshot for the expanded panel and creates top-10 labeling sheets.
- `expanded-benchmark-*/expanded-benchmark-freeze.json` - dated live expanded benchmark output.
- `expanded-benchmark-*/annotator_A_top10.csv` and `annotator_B_top10.csv` - independent relevance-labeling sheets.
- `score_relevance_labels.py` - legacy two-annotator scorer retained for audit history.
- `score_three_annotator_labels.py` - computes three-annotator validation metrics for the 80-query pooled record set, including final consensus labels, P@10, strict P@10, nDCG@10, MRR, agreement, Fleiss' kappa, bootstrap confidence intervals, paired tests, and adjudicated no-majority rows.
- `prepare_label_reuse.py` - creates the 80-query pooled labeling file and reuses a prior consensus label only when normalized query, source, and accession all match exactly.
- `analyze_benchmark_upgrade.py` - applies the predeclared completion gates and writes a benchmark analysis without inferring relevance from counts.
- `independent-labeling-protocol.md` - legacy protocol for the two-annotator relevance check; the current submission uses the 80-query three-annotator evaluator packet under `independent_relevance_evaluation_80q_20260703/`.
- `run_source_reliability.py` - cold/warm live-source reliability runner.
- `run_concept_ablation.py` - paired runner comparing candidate sets with and without deterministic concept expansion.
- `concept-ablation-*/concept-ablation-summary.csv` - per-query off/on counts, overlap, and elapsed time.
- `concept-ablation-*/concept-ablation-aggregate.csv` - aggregate paired candidate-set summary.
- `concept-ablation-*/concept-ablation-raw.json` and `run_manifest.json` - raw responses and file checksums.
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
python BioScouter_Paper\v3-concise\reproducibility\rerun_benchmark.py --base-url http://127.0.0.1:8001 --token "<JWT>"
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
python BioScouter_Paper\v3-concise\reproducibility\controlled_evaluation.py --base-url http://127.0.0.1:8001 --queries BioScouter_Paper\v3-concise\reproducibility\benchmark-queries-80.csv --expected-query-count 80 --candidate-depth 100 --top-k 10 --semantic-depth 50
```

This writes raw responses, ranked outputs, run manifests, pool mappings, summary
CSVs, and blinded annotation sheets. The current manuscript uses the completed
80-query three-annotator packet in
`independent_relevance_evaluation_80q_20260703/` to report independent
P@10, strict P@10, nDCG@10, MRR, agreement, and adjudicated consensus metrics.

Runs created before per-mode timing scopes were added record BioScouter source
search request time only and exclude local semantic reranking. Their elapsed
fields are retained for audit but must not be reported as end-to-end hybrid
latency. The manuscript uses the separate frozen ten-query SSE benchmark for
time-to-first-result and total-search latency claims.

## Expanded Benchmark Example

```powershell
python BioScouter_Paper\v3-concise\reproducibility\run_expanded_benchmark.py --base-url http://127.0.0.1:8001 --max-results 100 --top-n 10
```

For the legacy two-annotator expanded panel, after both annotators complete the
generated top-10 CSV files, run:

```powershell
python BioScouter_Paper\v3-concise\reproducibility\score_relevance_labels.py --annotator-a expanded-benchmark-YYYYMMDD-HHMMSS\annotator_A_top10.csv --annotator-b expanded-benchmark-YYYYMMDD-HHMMSS\annotator_B_top10.csv
```

For the current 80-query three-annotator packet, run:

```powershell
python BioScouter_Paper\v3-concise\reproducibility\score_three_annotator_labels.py
```

Do not update manuscript independent-validation metrics unless annotator sheets,
adjudication rows, and scoring outputs are all complete.

## Concept Expansion Ablation

The deterministic concept-normalization feature is evaluated as a paired
candidate-set ablation, not as a standalone relevance metric:

```powershell
python BioScouter_Paper\v3-concise\reproducibility\run_concept_ablation.py --base-url http://127.0.0.1:8001 --queries BioScouter_Paper\v3-concise\reproducibility\benchmark-queries-80.csv --max-results 10
```

The output records result counts and candidate overlap for each query. Any
claim about improved relevance still requires the independent top-10 labels.
