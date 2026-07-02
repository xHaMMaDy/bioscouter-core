"""
GEO (Gene Expression Omnibus) Adapter
Wraps the existing NCBI service to provide UnifiedDataset output.
"""

import asyncio
import re
from typing import List, Optional, Set

import structlog

from .base import BaseSourceAdapter, COMMON_STOPWORDS
from ..models.unified import (
    DataSource,
    OmicsType,
    AssayType,
    UnifiedDataset,
    GEOExtension,
    DownloadLink,
    DownloadFileType,
)
from ..services.ncbi import get_ncbi_service, NCBIService, NCBIError

logger = structlog.get_logger(__name__)

# GEO-specific stopwords (supplement to COMMON_STOPWORDS)
GEO_STOPWORDS: Set[str] = frozenset({
    'geo', 'ncbi', 'gse', 'gsm', 'gpl', 'expression', 'microarray',
    'chip', 'seq', 'sequencing', 'profiling', 'platform'
})


class GEOAdapter(BaseSourceAdapter):
    """
    Adapter for NCBI Gene Expression Omnibus (GEO).
    Normalizes GEO datasets to UnifiedDataset format.
    """
    
    def __init__(self, ncbi_service: Optional[NCBIService] = None):
        super().__init__()
        self._ncbi = ncbi_service or get_ncbi_service()
    
    @property
    def source(self) -> DataSource:
        return DataSource.GEO
    
    @property
    def supported_omics(self) -> List[OmicsType]:
        return [
            OmicsType.TRANSCRIPTOMICS,
            OmicsType.EPIGENOMICS,
            OmicsType.SINGLE_CELL,
        ]
    
    async def search(
        self,
        query: str,
        max_results: int = 100,
        organism: Optional[str] = None,
        soft_organism_filter: bool = True,
        **kwargs
    ) -> List[UnifiedDataset]:
        """
        Search GEO and return unified datasets.
        
        Args:
            query: Search query
            max_results: Maximum results to return (default increased to 100 for better recall)
            organism: Optional organism filter
            soft_organism_filter: If True, also search without organism to improve recall
            
        Returns:
            List of UnifiedDataset objects
        """
        try:
            # Use existing NCBI service to search with improved recall settings
            geo_datasets = await self._ncbi.search_and_fetch(
                query=query,
                organism=organism,
                max_results=max_results,
                soft_organism_filter=soft_organism_filter,
            )
            
            # Convert to unified format with relevance scoring
            unified = []
            for dataset in geo_datasets:
                try:
                    ud = self._to_unified(dataset, query)
                    if ud:
                        unified.append(ud)
                except Exception as e:
                    logger.warning("Failed to convert GEO dataset", error=str(e))
                    continue
            
            # Sort by relevance score
            unified.sort(key=lambda x: x.relevance_score, reverse=True)
            
            logger.info("GEO search complete", results=len(unified))
            return unified
            
        except NCBIError as e:
            logger.error("GEO search failed", error=str(e))
            raise
    
    async def get_dataset(self, accession: str) -> Optional[UnifiedDataset]:
        """
        Fetch a single GEO dataset by accession.
        
        Args:
            accession: GSE accession (e.g., "GSE123456")
            
        Returns:
            UnifiedDataset if found
        """
        try:
            dataset = await self._ncbi.fetch_by_gse_id(accession)
            if dataset:
                return self._to_unified(dataset, query=None)
            return None
        except NCBIError as e:
            logger.error("Failed to fetch GEO dataset", accession=accession, error=str(e))
            return None
    
    def _to_unified(self, dataset, query: Optional[str] = None) -> Optional[UnifiedDataset]:
        """
        Convert a GEO Dataset (old schema) to UnifiedDataset.
        
        Args:
            dataset: GEO Dataset object
            query: Original search query for relevance scoring
        """
        if not dataset or not dataset.gse_id:
            return None
        
        # Detect omics type and assay from dataset type field
        omics_type, assay_types = self._detect_omics_and_assay(dataset)
        
        # Build extension with GEO-specific data
        geo_ext = GEOExtension(
            gse_id=dataset.gse_id,
            platform=dataset.platform or None,
            platform_id=dataset.platform or None,
            series_type=dataset.dataset_type or None,
            sra_ids=dataset.sra_ids or [],
        )
        
        # Compute relevance score if query provided, otherwise use existing or compute basic score
        if query:
            relevance = self._compute_keyword_relevance(
                query, dataset.title, dataset.summary, GEO_STOPWORDS
            )
        elif dataset.relevance_score and dataset.relevance_score > 0:
            relevance = dataset.relevance_score
        else:
            relevance = 0.5  # Default for direct accession lookups
        
        # Generate download links for GEO dataset
        download_links = self._generate_download_links(dataset.gse_id, dataset.sra_ids)
        
        return UnifiedDataset(
            # Identifiers
            id=self.build_unified_id(dataset.gse_id),
            accession=dataset.gse_id,
            source=self.source,
            source_url=f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={dataset.gse_id}",
            
            # Classification
            omics_type=omics_type,
            assay_types=assay_types,
            
            # Basic metadata
            title=dataset.title or "",
            description=dataset.summary,
            organism=self._normalize_organism(dataset.organism),
            
            # Quantitative
            sample_count=dataset.sample_count or 0,
            sample_count_display=str(dataset.sample_count) if dataset.sample_count else "N/A",
            
            # Dates
            submission_date=dataset.submission_date,
            
            # Publications
            pubmed_ids=dataset.pubmed_ids or [],
            
            # AI/Search scores
            relevance_score=relevance,
            quality_score=dataset.quality_score,
            
            # Download links
            download_links=download_links,
            
            # Extensions
            extensions={"geo": geo_ext.model_dump()},
        )
    
    def _generate_download_links(self, gse_id: str, sra_ids: Optional[List[str]] = None) -> List[DownloadLink]:
        """
        Generate download links for a GEO dataset.
        
        GEO FTP structure:
        - ftp.ncbi.nlm.nih.gov/geo/series/{GSEnnn}nnn/{GSE}/
          - soft/ - SOFT formatted family file
          - matrix/ - Series matrix file (expression matrix)
          - suppl/ - Supplementary files (raw data, processed)
        
        Args:
            gse_id: GEO Series accession (e.g., "GSE123456")
            sra_ids: Optional list of SRA accessions
            
        Returns:
            List of DownloadLink objects
        """
        links = []
        
        # Extract GSE prefix for FTP path (e.g., "GSE123" from "GSE123456")
        # Pattern: GSEnnn becomes the folder name
        gse_prefix = gse_id[:len(gse_id) - 3] if len(gse_id) > 6 else gse_id[:3]
        ftp_base = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{gse_prefix}nnn/{gse_id}"
        
        # SOFT file - Contains all metadata and data in SOFT format
        links.append(DownloadLink(
            url=f"{ftp_base}/soft/{gse_id}_family.soft.gz",
            file_type=DownloadFileType.SOFT,
            file_name=f"{gse_id}_family.soft.gz",
            protocol="https",
            description="SOFT formatted family file (metadata + data)",
        ))
        
        # Series Matrix file - Tab-delimited expression matrix
        # Note: Some datasets have multiple matrix files, but most have just one
        links.append(DownloadLink(
            url=f"{ftp_base}/matrix/{gse_id}_series_matrix.txt.gz",
            file_type=DownloadFileType.MATRIX,
            file_name=f"{gse_id}_series_matrix.txt.gz",
            protocol="https",
            description="Series matrix file (expression values)",
        ))
        
        # Supplementary files directory
        links.append(DownloadLink(
            url=f"{ftp_base}/suppl/",
            file_type=DownloadFileType.SUPPLEMENTARY,
            protocol="https",
            description="Supplementary files (raw/processed data)",
        ))
        
        # SRA raw data link if available
        if sra_ids and len(sra_ids) > 0:
            # Link to SRA Run Selector for bulk download
            links.append(DownloadLink(
                url=f"https://www.ncbi.nlm.nih.gov/Traces/study/?acc={gse_id}",
                file_type=DownloadFileType.RAW,
                protocol="https",
                description=f"SRA Run Selector ({len(sra_ids)} runs) - Raw sequencing data",
            ))
        
        return links
    
    def _detect_omics_and_assay(self, dataset) -> tuple[OmicsType, List[AssayType]]:
        """
        Detect omics type and assay type from GEO dataset metadata.
        Uses dataset_type field and title/summary keywords.
        """
        text = f"{dataset.title or ''} {dataset.summary or ''} {dataset.dataset_type or ''}".lower()
        
        assay_types = []
        omics_type = OmicsType.TRANSCRIPTOMICS  # Default
        
        # Single-cell detection
        if any(kw in text for kw in ["single-cell", "single cell", "scrna", "10x genomics", "dropseq", "smart-seq"]):
            omics_type = OmicsType.SINGLE_CELL
            if "scrna" in text or "scrnaseq" in text:
                assay_types.append(AssayType.SCRNA_SEQ)
            elif "snrna" in text:
                assay_types.append(AssayType.SNRNA_SEQ)
            elif "spatial" in text:
                assay_types.append(AssayType.SPATIAL_TRANSCRIPTOMICS)
            else:
                assay_types.append(AssayType.SCRNA_SEQ)
            return omics_type, assay_types
        
        # Epigenomics detection
        if any(kw in text for kw in ["chip-seq", "chipseq", "atac-seq", "atacseq", "methylation", "bisulfite", "histone", "dnase"]):
            omics_type = OmicsType.EPIGENOMICS
            if "chip-seq" in text or "chipseq" in text:
                assay_types.append(AssayType.CHIP_SEQ)
            if "atac-seq" in text or "atacseq" in text:
                assay_types.append(AssayType.ATAC_SEQ)
            if "methylation" in text or "bisulfite" in text:
                assay_types.append(AssayType.BISULFITE_SEQ)
            if "dnase" in text:
                assay_types.append(AssayType.DNASE_SEQ)
            if "cut&run" in text or "cut-and-run" in text:
                assay_types.append(AssayType.CUT_AND_RUN)
            if "cut&tag" in text or "cut-and-tag" in text:
                assay_types.append(AssayType.CUT_AND_TAG)
            if not assay_types:
                assay_types.append(AssayType.CHIP_SEQ)
            return omics_type, assay_types
        
        # Transcriptomics (RNA-seq vs microarray)
        omics_type = OmicsType.TRANSCRIPTOMICS
        
        if "rna-seq" in text or "rnaseq" in text or "rna seq" in text:
            assay_types.append(AssayType.RNA_SEQ)
        elif "mrna-seq" in text or "mrna seq" in text:
            assay_types.append(AssayType.MRNA_SEQ)
        elif "microarray" in text or "array" in text or "affymetrix" in text or "agilent" in text:
            assay_types.append(AssayType.MICROARRAY)
        else:
            # Check dataset type field
            ds_type = (dataset.dataset_type or "").lower()
            if "expression profiling by array" in ds_type:
                assay_types.append(AssayType.MICROARRAY)
            elif "expression profiling by high throughput sequencing" in ds_type:
                assay_types.append(AssayType.RNA_SEQ)
            else:
                assay_types.append(AssayType.RNA_SEQ)  # Default assumption
        
        return omics_type, assay_types


# Factory function
_geo_adapter: Optional[GEOAdapter] = None


def get_geo_adapter() -> GEOAdapter:
    """Get or create GEO adapter instance."""
    global _geo_adapter
    if _geo_adapter is None:
        _geo_adapter = GEOAdapter()
    return _geo_adapter
