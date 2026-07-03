# BioScouter Core

BioScouter Core is the public scientific/reproducibility package for BioScouter. It contains the reusable parts needed to inspect, run, and reproduce the manuscript methods without exposing the production web platform.

## Included

- Public-source adapters for GEO, PRIDE, ENCODE, MetaboLights, Metabolomics Workbench, GDC, CPTAC/PDC, CellxGene, HCA, MGnify, MG-RAST, ENA, MassIVE, and GTEx.
- Unified multi-omics dataset schemas.
- Deterministic omics detection and curated concept expansion.
- Federated adapter orchestration, result merging, ranking, faceting, and metadata-readiness scoring.
- Benchmark queries, frozen benchmark outputs, independent relevance-evaluation materials, scoring scripts, and protocol notes.
- CLI and tests for validating the package.

## Not Included

This repository does not include the production BioScouter frontend, authentication, user accounts, credits, payments, admin panel, deployment configuration, private observability, LLM-provider orchestration, or production secrets.

That split is intentional: reviewers and readers can inspect the scientific adapters, normalization, ranking helpers, and evaluation scripts, while the deployed service at `bioscouter.com` remains a managed production platform.

## Install

```bash
python -m pip install -e ".[dev]"
```

Optional environment variables:

```bash
set BIOSCOUTER_NCBI_EMAIL=your.email@example.com
set BIOSCOUTER_NCBI_API_KEY=
set BIOSCOUTER_MGRAST_API_KEY=
```

NCBI requests work without an API key, but NCBI recommends setting a real email.

## Quick Search

```bash
bioscouter-core-search "TMT proteomics breast cancer" --source pride --max-results 10
```

JSON output:

```bash
bioscouter-core-search "single-cell lung cancer atlas" --max-results 5 --json
```

## Python Usage

```python
import asyncio
from bioscouter_core import BioScouterCoreSearch

async def main():
    search = BioScouterCoreSearch()
    try:
        results = await search.search("TMT proteomics breast cancer", sources=["pride"], max_results=10)
        print(results.total_results)
    finally:
        await search.close()

asyncio.run(main())
```

## Reproducibility

See:

- `reproducibility/paper/README.md`
- `reproducibility/paper/expanded-benchmark-queries.csv`
- `reproducibility/paper/expanded-benchmark-20260702-153134/expanded-benchmark-freeze.json`
- `reproducibility/paper/independent_relevance_evaluation_80q_20260703/README.md`
- `reproducibility/paper/score_three_annotator_labels.py`
- `docs/REPRODUCIBILITY.md`

## Code Availability Statement

Suggested manuscript wording is in `docs/CODE_AVAILABILITY.md`.

## Tests

```bash
python -m pytest
```
