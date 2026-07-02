# Return Checklist

Use this checklist before using the labels in the manuscript.

## Before Sending

- [ ] Send only `01_Evaluator_A/` to Evaluator A.
- [ ] Send only `02_Evaluator_B/` to Evaluator B.
- [ ] Ask them to work independently.
- [ ] Ask them not to edit any metadata columns.
- [ ] Ask them to return the completed CSV file.

## When Files Return

- [ ] Save Evaluator A's completed CSV as:
  `00_For_Ibrahim/returned_labels/BioScouter_Evaluator_A_Labeling_Sheet.csv`
- [ ] Save Evaluator B's completed CSV as:
  `00_For_Ibrahim/returned_labels/BioScouter_Evaluator_B_Labeling_Sheet.csv`
- [ ] Confirm both files have 300 labels.
- [ ] Confirm every label is one of `0`, `1`, or `2`.
- [ ] Run the scoring helper.
- [ ] Review disagreements.
- [ ] Create an adjudicated label file if you decide to report final adjudicated metrics.

## Manuscript Details To Report

- Number of queries: 30
- Results labeled per query: top 10
- Total rows per evaluator: 300
- Label scale: 0/1/2
- Evaluator background: bioinformatics/omics data-search experience
- Whether evaluators worked independently: yes
- Agreement statistics: report after scoring
- P@10/nDCG@10: report after scoring

