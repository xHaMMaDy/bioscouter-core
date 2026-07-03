# BioScouter Benchmark Upgrade Analysis

Controlled run: `controlled-evaluation-80-20260703-030544`

## Decision Gate

- PASS: at least 60 hybrid queries with results
- PASS: required internal modes present
- PASS: required baselines present
- PASS: traceable outputs present

Overall gate: **PASS**

## System Summary

| System | Requests | With results | Returned top-10 | Unique | Sources |
|---|---:|---:|---:|---:|---:|
| bioscouter_keyword | 80/80 | 70/80 | 667 | 583 | 9 |
| bioscouter_embedding | 80/80 | 70/80 | 667 | 616 | 10 |
| bioscouter_hybrid | 80/80 | 80/80 | 800 | 704 | 11 |
| native_source_api | 80/80 | 62/80 | 455 | 391 | 10 |
| omicsdi | 80/80 | 77/80 | 699 | 695 | 16 |

## Concept-Normalization Ablation

- Paired queries: 80
- More/fewer/same returned counts: 0/1/79
- Median candidate-set Jaccard overlap: 1.0

## Interpretation Boundary

The 80-query controlled outputs are not independently labeled. P@10, strict P@10, nDCG@10, MRR, and relevant-dataset superiority must not be inferred from returned counts or latency.

Use the separately completed and adjudicated 30-query, 300-record evaluation.
