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
- `score_three_annotator_labels.py`
- `run_expanded_benchmark.py`
- `run_concept_ablation.py`
- `run_source_reliability.py`

## Independent Relevance Evaluation

The `reproducibility/paper/independent_relevance_evaluation_80q_20260703`
folder contains the current three-annotator evaluator packet, returned labels,
scoring output, adjudication file, and summary.

Important files:

- `README.md`
- `RETURNED_LABELS/BioScouter_80Q_Annotator_A_Labeling_Sheet_COMPLETED.csv`
- `RETURNED_LABELS/BioScouter_80Q_Annotator_B_Labeling_Sheet_COMPLETED.csv`
- `RETURNED_LABELS/BioScouter_80Q_Annotator_C_Labeling_Sheet_COMPLETED.csv`
- `COORDINATOR_ONLY/ADJUDICATION_80Q_COMPLETED_Amr.csv`
- `COORDINATOR_ONLY/three_annotator_scoring_output/final_consensus_system_summary.csv`
- `COORDINATOR_ONLY/three_annotator_scoring_output/three_annotator_independent_validation_report.json`

The headline BioScouter hybrid values are P@10 75.8%, strict P@10 23.6%,
and nDCG@10 62.2% over 80 queries and 2,348 unique pooled records after
adjudication of 148 no-majority rows.

## Expected Variability

Live API searches may change as external repositories update their records and ranking endpoints. The frozen JSON and CSV outputs are therefore the traceable evidence for manuscript numbers; live reruns are a robustness check, not a guarantee of byte-identical results.
