# Independent Evaluation Results Summary

Status: valid and complete at the independent-label and draft adjudication stages.

## Input Files Checked

- `returned_labels/BioScouter_Evaluator_A_Labeling_Sheet.csv`
- `returned_labels/BioScouter_Evaluator_B_Labeling_Sheet.csv`

Both files contain:

- 30 benchmark queries
- 300 top-10 result rows
- complete labels using only `0`, `1`, or `2`

## Label Distribution

| Evaluator | Label 2 | Label 1 | Label 0 |
|---|---:|---:|---:|
| A | 149 | 102 | 49 |
| B | 130 | 114 | 56 |

## Pre-Adjudication Metrics

| Evaluator | Mean P@10 | Mean strict P@10 | Mean nDCG@10 |
|---|---:|---:|---:|
| A | 83.7% | 49.7% | 89.7% |
| B | 81.3% | 43.3% | 92.0% |

## Adjudicated Consensus Metrics

| Metric | Value |
|---|---:|
| Consensus P@10 | 80.0% |
| Consensus strict P@10 | 46.0% |
| Consensus nDCG@10 | 90.8% |

## Agreement

- Paired labels: 300
- Exact agreement: 75.7%
- Within-one-category agreement: 98.3%
- Cohen's kappa, unweighted: 0.61
- Cohen's kappa, linear-weighted: 0.67

## Disagreements

There were 73 discordant rows.

Adjudication file:

`ADJUDICATION_TO_COMPLETE.csv`

The `adjudicated_label` column has been filled with draft adjudication labels and short rationale notes.

Original backup before adjudication:

`ADJUDICATION_TO_COMPLETE.before_draft_adjudication.csv`

## Manuscript Use

These results are suitable to report as independent relevance metrics after the corresponding author reviews and accepts the adjudication file. Before submission, inspect `ADJUDICATION_TO_COMPLETE.csv` and confirm the adjudicated labels are acceptable.
