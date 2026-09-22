# BioScouter Core 0.3.1

Release date: 22 September 2026

This release records the public scientific core and reproducibility materials used for the BioScouter manuscript revision. It updates the public author metadata to Ibrahim Abdelkarim Hammad, Amr M. Alhfnawy, and Sameh E. Hassanein, and preserves the controlled evaluation artifacts and scripts already included in the repository. Version 0.3.0 already exists on GitHub; this patch release contains the corrected author metadata and revision evidence.

## Reproducibility scope

- The controlled T01-T80 evaluation uses keyword candidates, local embedding reranking, and a hybrid union with a frozen local semantic corpus.
- The controlled run disabled production pgvector retrieval, automatic indexing, and concept expansion.
- The frozen local corpus contains 1,464 deduplicated records from 16 construction queries, 10 source labels, and eight omics categories.
- The L10 source-count reconciliation is included in `reproducibility/paper/source_count_reconciliation_20260922.csv`.
- Historical provider-level LLM cost, token, and latency telemetry is not part of the frozen benchmark artifacts.

## Production boundary

This repository does not contain the hosted frontend, authentication, accounts, credits, payment controls, deployment configuration, private observability, hosted Supabase RPC definitions, or production secrets. Those components are not required to inspect the public adapters, schemas, deterministic query handling, scoring scripts, or frozen evaluation outputs.

## Verification

The repository test suite passes before tagging this release. The exact release commit and the reproducibility archive hash are recorded in the accompanying release manifest used for the immutable archive.
