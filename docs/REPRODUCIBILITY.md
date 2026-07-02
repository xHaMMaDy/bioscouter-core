# Reproducibility

This repository includes two reproducibility layers.

## Frozen Benchmarks

The `reproducibility/paper` folder contains benchmark query sets, frozen raw outputs, and scripts used to regenerate summary metrics where supported by public APIs.

Important files:

- `benchmark-queries.csv`
- `expanded-benchmark-queries.csv`
- `expanded-benchmark-20260702-153134/expanded-benchmark-freeze.json`
- `expanded-benchmark-20260702-153134/expanded-benchmark-summary.csv`
- `score_relevance_labels.py`
- `run_expanded_benchmark.py`
- `run_concept_ablation.py`
- `run_source_reliability.py`

## Independent Relevance Evaluation

The `reproducibility/independent_relevance_evaluation` folder contains the evaluator packet, returned labels, scoring output, adjudication file, and summary.

Important files:

- `INDEPENDENT_EVALUATION_RESULTS_SUMMARY.md`
- `returned_labels/BioScouter_Evaluator_A_Labeling_Sheet.csv`
- `returned_labels/BioScouter_Evaluator_B_Labeling_Sheet.csv`
- `scoring_output/adjudicated_system_summary.csv`
- `scoring_output/adjudicated_per_query_metrics.csv`
- `score_returned_top10_labels.py`

## Expected Variability

Live API searches may change as external repositories update their records and ranking endpoints. The frozen JSON and CSV outputs are therefore the traceable evidence for manuscript numbers; live reruns are a robustness check, not a guarantee of byte-identical results.

