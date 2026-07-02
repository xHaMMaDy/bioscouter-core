# Architecture

BioScouter Core exposes the scientific search path as a small Python package.

```mermaid
flowchart LR
    Q["User query"] --> O["Omics detection and optional concept expansion"]
    O --> S["Source selection"]
    S --> A["Public repository adapters"]
    A --> N["UnifiedDataset normalization"]
    N --> M["Cross-source merge and deduplication"]
    M --> R["Sort, facet, metadata-readiness scoring"]
    R --> E["UnifiedSearchResponse"]
```

## Package Layout

- `src/bioscouter_core/adapters`: public repository clients and normalizers.
- `src/bioscouter_core/models`: Pydantic schemas and source metadata.
- `src/bioscouter_core/services`: deterministic helper services.
- `src/bioscouter_core/orchestrator.py`: standalone federated search runner.
- `src/bioscouter_core/cli.py`: command-line interface.
- `reproducibility`: manuscript benchmark and relevance-labeling materials.

## Boundary

The core package does not own product behavior. Auth, credits, billing, admin settings, private LLM orchestration, frontend routes, and deployment configuration remain outside this repository.

The public package is sufficient to inspect how source APIs are queried, how records are normalized, how duplicated records are merged, and how frozen benchmark claims were generated and scored.

