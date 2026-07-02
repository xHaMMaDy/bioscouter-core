# Scoring After The Evaluators Return Labels

## 1. Put Returned Files Here

Create this folder:

`00_For_Ibrahim/returned_labels/`

Place the completed files inside it with these exact names:

- `BioScouter_Evaluator_A_Labeling_Sheet.csv`
- `BioScouter_Evaluator_B_Labeling_Sheet.csv`

## 2. Run The Scoring Helper

From the packet root folder, run:

```powershell
python .\00_For_Ibrahim\score_returned_top10_labels.py `
  --annotator-a .\00_For_Ibrahim\returned_labels\BioScouter_Evaluator_A_Labeling_Sheet.csv `
  --annotator-b .\00_For_Ibrahim\returned_labels\BioScouter_Evaluator_B_Labeling_Sheet.csv `
  --out-dir .\00_For_Ibrahim\scoring_output
```

If Python is not on PATH, run it with the project Python:

```powershell
D:\AI Agents\VS Code\BioScouter\backend\venv\Scripts\python.exe .\00_For_Ibrahim\score_returned_top10_labels.py `
  --annotator-a .\00_For_Ibrahim\returned_labels\BioScouter_Evaluator_A_Labeling_Sheet.csv `
  --annotator-b .\00_For_Ibrahim\returned_labels\BioScouter_Evaluator_B_Labeling_Sheet.csv `
  --out-dir .\00_For_Ibrahim\scoring_output
```

## 3. Outputs

The helper writes:

- `label_validation.json` - missing/invalid label checks.
- `agreement_summary.json` - exact agreement, within-one agreement, Cohen's kappa.
- `per_query_metrics.csv` - P@10, strict P@10, and nDCG@10 per query for each evaluator.
- `system_summary.csv` - overall means across all 30 queries.
- `disagreements_for_adjudication.csv` - rows where Evaluator A and Evaluator B disagree.

## 4. How To Use In The Manuscript

Use the independent results only if both returned files are complete.

Report:

- two independent evaluators
- 30 benchmark queries
- 300 top-10 rows labeled per evaluator
- three-level relevance rubric
- exact agreement and kappa
- P@10 and nDCG@10

If disagreements are adjudicated, clearly say how adjudication was performed.

