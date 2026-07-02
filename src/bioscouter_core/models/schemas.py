"""Small legacy GEO schemas used by the public NCBI/GEO adapter.

The main public search API uses ``models.unified``. These classes are retained
only for compatibility with the NCBI helper that normalizes GEO records before
they are converted into ``UnifiedDataset``.
"""

from __future__ import annotations

from enum import Enum
from typing import Optional

from pydantic import BaseModel, Field


class DownloadFileType(str, Enum):
    """Types of downloadable files."""

    SOFT = "soft"
    MATRIX = "matrix"
    SUPPLEMENTARY = "supplementary"
    RAW = "raw"
    PROCESSED = "processed"
    METADATA = "metadata"
    OTHER = "other"


class DownloadLink(BaseModel):
    """A downloadable file link for a GEO dataset."""

    url: str = Field(..., description="Direct download URL")
    file_type: DownloadFileType = Field(default=DownloadFileType.OTHER, description="Type of file")
    file_name: Optional[str] = Field(default=None, description="Filename if known")
    file_size_bytes: Optional[int] = Field(default=None, description="File size in bytes if known")
    protocol: str = Field(default="https", description="Download protocol")
    description: Optional[str] = Field(default=None, description="Description of the file")
    needs_refresh: bool = Field(default=False, description="Whether URL needs refresh before download")


class Dataset(BaseModel):
    """A normalized GEO dataset result returned by the NCBI helper."""

    gse_id: str = Field(..., description="GEO Series accession, e.g. GSE12345")
    title: str = Field(..., description="Dataset title")
    summary: str = Field(default="", description="Dataset summary or abstract")
    organism: str = Field(default="Unknown", description="Organism studied")
    sample_count: Optional[int] = Field(default=None, description="Number of samples")
    sample_count_display: str = Field(default="Unknown", description="Display string for sample count")
    platform: str = Field(default="", description="Platform used")
    dataset_type: str = Field(default="", description="Dataset type")
    submission_date: Optional[str] = Field(default=None, description="Date submitted")
    geo_url: str = Field(default="", description="Direct link to GEO page")
    relevance_score: float = Field(default=0.0, ge=0.0, le=1.0, description="Computed relevance")
    pubmed_ids: list[str] = Field(default_factory=list, description="Associated PubMed IDs")
    sra_ids: list[str] = Field(default_factory=list, description="Associated SRA run IDs")
    quality_score: Optional[float] = Field(default=None, ge=0.0, le=1.0)
    quality_breakdown: Optional[dict] = Field(default=None)
    download_links: list[DownloadLink] = Field(default_factory=list, description="Direct download links")

    def model_post_init(self, __context) -> None:
        if not self.geo_url:
            self.geo_url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={self.gse_id}"
        self.sample_count_display = str(self.sample_count) if self.sample_count is not None else "Unknown"

