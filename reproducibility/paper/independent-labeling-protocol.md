# Independent Relevance Labeling Protocol

This protocol is intended for the expanded BioScouter benchmark before final journal submission.

## Files

- `expanded-benchmark-queries.csv` - 30-query expanded benchmark panel.
- `expanded-benchmark-*/annotator_A_top10.csv` - top-10 results for annotator A.
- `expanded-benchmark-*/annotator_B_top10.csv` - top-10 results for annotator B.
- `score_relevance_labels.py` - scoring script for P@10, nDCG@10, agreement, and adjudication.

## Annotators

Use two independent annotators with bioinformatics or omics data-search experience. Annotators should label results independently without seeing each other's labels.

Suggested minimum:

- Annotator A: corresponding author or project author.
- Annotator B: supervisor, colleague, or domain expert not involved in creating the initial labels.

## Label Definitions

- `2` - Directly relevant to the requested omics type and biological condition. The record clearly matches the query's core intent.
- `1` - Partially relevant. The record is related but misses one requested constraint, such as organism, tissue, disease, assay, or sample threshold.
- `0` - Irrelevant, wrong omics type, wrong biological context, or insufficient metadata to judge relevance.

## Instructions

1. Open the assigned annotator CSV.
2. For every row, inspect the query, title, description, source, omics type, organism, tissue, disease, and source URL when needed.
3. Enter only `0`, `1`, or `2` in `annotator_label`.
4. Use `annotator_notes` for short explanations when a label is uncertain.
5. Do not edit query IDs, ranks, accessions, titles, or source metadata.
6. After both annotators finish, run:

```powershell
python score_relevance_labels.py --annotator-a expanded-benchmark-YYYYMMDD-HHMMSS\annotator_A_top10.csv --annotator-b expanded-benchmark-YYYYMMDD-HHMMSS\annotator_B_top10.csv
```

7. Resolve disagreements in `adjudication_sheet.csv`; report both raw agreement and adjudicated P@10/nDCG@10 in the manuscript if completed before submission.

## Reporting

Report:

- Number of queries.
- Number of top-10 results labeled.
- Annotator expertise.
- Label rubric.
- Mean P@10 and nDCG@10.
- Percent agreement and Cohen's kappa.
- How disagreements were adjudicated.

Do not claim independent validation unless both annotator files are completed and scored.
