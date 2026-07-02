# Models module
from .schemas import (
    Dataset,
    DownloadFileType,
    DownloadLink,
)
from .unified import (
    OmicsType,
    DataSource,
    AssayType,
    UnifiedDataset,
    UnifiedSearchQuery,
    UnifiedSearchStatus,
    UnifiedSearchResponse,
    detect_omics_types,
    get_sources_for_omics,
)

__all__ = [
    # Legacy GEO-specific models
    "Dataset",
    "DownloadFileType",
    "DownloadLink",
    # Unified multi-omics models
    "OmicsType",
    "DataSource",
    "AssayType",
    "UnifiedDataset",
    "UnifiedSearchQuery",
    "UnifiedSearchStatus",
    "UnifiedSearchResponse",
    "detect_omics_types",
    "get_sources_for_omics",
]
