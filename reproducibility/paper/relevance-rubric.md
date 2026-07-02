# Relevance Labeling Rubric

The v3 manuscript uses a three-level manual relevance rubric for Table S3.

## Labels

- `2` - directly relevant to the requested omics type and biological condition.
- `1` - partially relevant but missing one requested constraint, such as organism, tissue, assay, or sample threshold.
- `0` - irrelevant, wrong biological context, wrong omics type, or insufficiently described.

## Metrics

- `P@10` counts labels `1` or `2` among the top ten ranked results.
- `nDCG@10` uses graded labels `0`, `1`, and `2`.
- The current labels were assigned during manuscript preparation by the corresponding author and should be treated as an initial single-reviewer evaluation.

## Independent Reviewer-Ready Upgrade

For the expanded benchmark, use `run_expanded_benchmark.py` to generate
`annotator_A_top10.csv` and `annotator_B_top10.csv`. Two annotators should label
the top ten results for each query independently, then use
`score_relevance_labels.py` to calculate P@10, nDCG@10, percent agreement, and
Cohen's kappa. Disagreements should be resolved in the generated adjudication
sheet before reporting final independent-validation metrics.
