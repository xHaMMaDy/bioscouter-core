"""Public scientific core for BioScouter.

This package contains the reusable dataset-discovery adapters, normalized
schemas, metadata-readiness scoring, deterministic query expansion, faceting,
and a lightweight federated-search orchestrator.
"""

from .models.unified import DataSource, OmicsType, UnifiedDataset, UnifiedSearchQuery, UnifiedSearchResponse
from .orchestrator import BioScouterCoreSearch, SearchConfig

__all__ = [
    "BioScouterCoreSearch",
    "DataSource",
    "OmicsType",
    "SearchConfig",
    "UnifiedDataset",
    "UnifiedSearchQuery",
    "UnifiedSearchResponse",
]

