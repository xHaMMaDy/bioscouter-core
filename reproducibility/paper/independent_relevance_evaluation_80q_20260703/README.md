# BioScouter 80-query three-annotator relevance evaluation

This folder contains the public relevance-evaluation materials for the
80-query controlled panel reported in the BioScouter manuscript.

## Contents

- `RETURNED_LABELS/` - completed blinded labels from annotators A, B, and C.
- `COORDINATOR_ONLY/pool_mapping_DO_NOT_SEND_TO_ANNOTATORS.csv` - mapping from
  pooled record IDs back to query/system/rank positions.
- `COORDINATOR_ONLY/ADJUDICATION_80Q_TO_COMPLETE.csv` - the 148 rows without a
  three-annotator majority.
- `COORDINATOR_ONLY/ADJUDICATION_80Q_COMPLETED_Amr.csv` - completed adjudication
  for the 148 no-majority rows.
- `COORDINATOR_ONLY/returned_label_validation_report.json` - validation of the
  returned label files before adjudication.
- `COORDINATOR_ONLY/three_annotator_agreement_summary.json` - pre-adjudication
  three-annotator agreement summary.
- `COORDINATOR_ONLY/three_annotator_scoring_output/` - final consensus labels,
  per-query metrics, system summaries, pairwise tests, agreement metrics, and
  validation report.

## Headline final values

- Panel: 80 controlled queries.
- Labeled pool: 2,348 unique pooled records.
- Annotators: 3 independent annotators.
- No-majority rows adjudicated: 148.
- BioScouter hybrid final P@10: 75.8%.
- BioScouter hybrid strict P@10: 23.6%.
- BioScouter hybrid nDCG@10: 62.2%.
- Fleiss' kappa: 0.244.

Run the scorer from the repository root with:

```powershell
python reproducibility\paper\score_three_annotator_labels.py --packet-dir reproducibility\paper\independent_relevance_evaluation_80q_20260703
```
