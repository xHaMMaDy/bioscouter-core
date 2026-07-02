"""
ENCODE Adapter
Searches the ENCODE (Encyclopedia of DNA Elements) project for epigenomics datasets.
API Docs: https://www.encodeproject.org/help/rest-api/
"""

import asyncio
import re
from typing import List, Optional, Any, Set
from datetime import datetime
from urllib.parse import urlencode

import httpx
import structlog

from .base import (
    BaseSourceAdapter,
    COMMON_STOPWORDS,
    DEFAULT_HTTP_TIMEOUT,
    build_default_limits,
)
from ..models.unified import (
    DataSource,
    OmicsType,
    AssayType,
    UnifiedDataset,
    ENCODEExtension,
    DownloadLink,
    DownloadFileType,
)

logger = structlog.get_logger(__name__)

# ENCODE API base URL
ENCODE_API_BASE = "https://www.encodeproject.org"

# ENCODE-specific stopwords (supplement to COMMON_STOPWORDS)
ENCODE_STOPWORDS: Set[str] = frozenset({
    'encode', 'epigenomics', 'epigenomic', 'chromatin', 'histone',
    'chip', 'seq', 'sequencing', 'experiment', 'biosample'
})


class ENCODEError(Exception):
    """ENCODE API error."""
    pass


class ENCODEAdapter(BaseSourceAdapter):
    """
    Adapter for ENCODE Project (epigenomics database).
    Uses ENCODE REST API.
    """
    
    def __init__(self, timeout: float = 30.0):
        super().__init__()
        self.timeout = timeout
        self._client: Optional[httpx.AsyncClient] = None
    
    @property
    def source(self) -> DataSource:
        return DataSource.ENCODE
    
    @property
    def supported_omics(self) -> List[OmicsType]:
        return [OmicsType.EPIGENOMICS, OmicsType.TRANSCRIPTOMICS]
    
    async def _get_client(self) -> httpx.AsyncClient:
        """Get or create HTTP client."""
        if self._client is None or self._client.is_closed:
            self._client = httpx.AsyncClient(
                timeout=self.timeout,
                limits=build_default_limits(),
                headers={
                    "Accept": "application/json",
                    "User-Agent": "BioScouter/1.0",
                }
            )
        return self._client
    
    async def close(self):
        """Close HTTP client."""
        if self._client and not self._client.is_closed:
            await self._client.aclose()
            self._client = None
    
    def _build_search_query(self, query: str, organism: Optional[str] = None) -> dict:
        """Build ENCODE search parameters."""
        params = {
            "type": "Experiment",
            "status": "released",
            "format": "json",
            "limit": "100",
        }
        
        # Parse query for common terms
        query_lower = query.lower()
        
        # Detect assay type from query
        assay_mapping = {
            "chip-seq": "ChIP-seq",
            "chipseq": "ChIP-seq",
            "chip seq": "ChIP-seq",
            "atac-seq": "ATAC-seq",
            "atacseq": "ATAC-seq",
            "atac seq": "ATAC-seq",
            "dnase-seq": "DNase-seq",
            "dnaseseq": "DNase-seq",
            "rna-seq": "RNA-seq",
            "rnaseq": "RNA-seq",
            "rna seq": "RNA-seq",
            "hi-c": "Hi-C",
            "hic": "Hi-C",
            "methylation": "WGBS",
            "bisulfite": "WGBS",
            "cage": "CAGE",
            "mint-chip": "Mint-ChIP-seq",
            "cut&run": "CUT&RUN",
            "cut&tag": "CUT&Tag",
        }
        
        for term, assay in assay_mapping.items():
            if term in query_lower:
                params["assay_title"] = assay
                break
        
        # Detect histone marks
        histone_marks = [
            "H3K4me1", "H3K4me3", "H3K27ac", "H3K27me3", "H3K36me3",
            "H3K9me3", "H3K9ac", "H3K79me2", "H4K20me1"
        ]
        for mark in histone_marks:
            if mark.lower() in query_lower:
                params["target.label"] = mark
                break
        
        # Detect transcription factors
        tf_pattern = r'\b([A-Z][A-Z0-9]+)\b'
        tf_matches = re.findall(tf_pattern, query)
        for tf in tf_matches:
            if len(tf) >= 2 and tf not in ["RNA", "DNA", "SEQ", "CHIP", "ATAC"]:
                params["target.label"] = tf
                break
        
        # Add organism filter
        if organism:
            organism_mapping = {
                "homo sapiens": "Homo sapiens",
                "human": "Homo sapiens",
                "mus musculus": "Mus musculus",
                "mouse": "Mus musculus",
            }
            params["replicates.library.biosample.donor.organism.scientific_name"] = (
                organism_mapping.get(organism.lower(), organism)
            )
        
        # Add general search term
        params["searchTerm"] = query
        
        return params
    
    async def _fetch_with_retry(self, url: str) -> httpx.Response:
        """Fetch URL with retry logic. Delegates to the shared base helper."""
        client = await self._get_client()
        return await super()._fetch_with_retry(client, url)
    
    def _detect_omics_type(self, experiment: dict) -> OmicsType:
        """Detect omics type from experiment metadata."""
        assay = (experiment.get("assay_title") or "").lower()
        
        if any(term in assay for term in ["rna-seq", "rnaseq", "rna seq", "cage", "rampage"]):
            return OmicsType.TRANSCRIPTOMICS
        
        return OmicsType.EPIGENOMICS
    
    def _detect_assay_type(self, experiment: dict) -> AssayType:
        """Detect assay type from experiment metadata."""
        assay = (experiment.get("assay_title") or "").lower()
        
        if "chip-seq" in assay or "chip seq" in assay:
            return AssayType.CHIP_SEQ
        if "atac-seq" in assay or "atac seq" in assay:
            return AssayType.ATAC_SEQ
        if "dnase-seq" in assay or "dnase seq" in assay:
            return AssayType.DNASE_SEQ
        if "rna-seq" in assay or "rna seq" in assay:
            return AssayType.RNA_SEQ
        if "hi-c" in assay or "hic" in assay:
            return AssayType.HI_C
        if "wgbs" in assay or "bisulfite" in assay or "methylation" in assay:
            return AssayType.WGBS
        if "mint-chip" in assay:
            return AssayType.MINT_CHIP
        if "cut&run" in assay or "cutnrun" in assay or "cut-and-run" in assay:
            return AssayType.CUT_AND_RUN
        if "cut&tag" in assay or "cutntag" in assay or "cut-and-tag" in assay:
            return AssayType.CUT_AND_TAG
        
        return AssayType.OTHER
    
    async def search(
        self,
        query: str,
        max_results: int = 50,
        organism: Optional[str] = None,
        **kwargs
    ) -> List[UnifiedDataset]:
        """
        Search ENCODE for epigenomics datasets.
        
        Args:
            query: Search query
            max_results: Maximum results to return
            organism: Optional organism filter
            
        Returns:
            List of UnifiedDataset objects
            
        Raises:
            ENCODEError: On API or network failures
        """
        try:
            params = self._build_search_query(query, organism)
            
            url = f"{ENCODE_API_BASE}/search/?{urlencode(params, doseq=True)}"
            
            logger.info("Searching ENCODE", query=query, params=params)
            
            # Use retry-enabled fetch
            response = await self._fetch_with_retry(url)
            response.raise_for_status()
            
            data = response.json()
            
            # Parse experiments from response
            experiments = data.get("@graph", [])
            
            logger.info("ENCODE response", experiment_count=len(experiments))
            
            unified = []
            for experiment in experiments[:max_results]:
                try:
                    ud = self._to_unified(experiment, query=query)
                    if ud:
                        unified.append(ud)
                except Exception as e:
                    logger.warning("Failed to parse ENCODE experiment", error=str(e))
                    continue
            
            # Sort by relevance score
            unified.sort(key=lambda x: x.relevance_score, reverse=True)
            
            logger.info("ENCODE search complete", results=len(unified))
            return unified
            
        except httpx.HTTPStatusError as e:
            logger.error("ENCODE API error", status=e.response.status_code, error=str(e))
            raise ENCODEError(f"ENCODE API error: {e.response.status_code}") from e
        except httpx.RequestError as e:
            logger.error("ENCODE request failed", error=str(e))
            raise ENCODEError(f"ENCODE request failed: {str(e)}") from e
    
    async def get_dataset(self, accession: str) -> Optional[UnifiedDataset]:
        """Fetch a single ENCODE experiment by accession."""
        try:
            # Normalize accession
            if not accession.upper().startswith("ENCSR"):
                accession = f"ENCSR{accession}"
            
            url = f"{ENCODE_API_BASE}/experiments/{accession}/?format=json"
            
            client = await self._get_client()
            response = await client.get(url)
            
            if response.status_code == 404:
                return None
            
            response.raise_for_status()
            experiment = response.json()
            
            return self._to_unified(experiment)
            
        except httpx.HTTPStatusError as e:
            if e.response.status_code == 404:
                return None
            logger.error("ENCODE API error", accession=accession, error=str(e))
            return None
        except Exception as e:
            logger.error("Failed to fetch ENCODE experiment", accession=accession, error=str(e))
            return None
    
    def _to_unified(self, experiment: dict, query: Optional[str] = None) -> Optional[UnifiedDataset]:
        """Convert an ENCODE experiment to UnifiedDataset."""
        try:
            accession = experiment.get("accession") or ""
            if not accession:
                return None
            
            # Build title from experiment details
            assay = experiment.get("assay_title") or ""
            target = experiment.get("target", {})
            target_label = target.get("label", "") if isinstance(target, dict) else ""
            biosample = experiment.get("biosample_ontology", {})
            biosample_term = biosample.get("term_name", "") if isinstance(biosample, dict) else ""
            
            title = f"{assay}"
            if target_label:
                title += f" of {target_label}"
            if biosample_term:
                title += f" in {biosample_term}"
            
            if not title.strip():
                title = f"ENCODE Experiment {accession}"
            
            # Description
            description = experiment.get("description") or ""
            
            # Organism
            organism = []
            replicates = experiment.get("replicates", [])
            for rep in replicates:
                library = rep.get("library", {})
                biosample_data = library.get("biosample", {}) if isinstance(library, dict) else {}
                donor = biosample_data.get("donor", {}) if isinstance(biosample_data, dict) else {}
                org = donor.get("organism", {}) if isinstance(donor, dict) else {}
                org_name = org.get("scientific_name", "") if isinstance(org, dict) else ""
                if org_name and org_name not in organism:
                    organism.append(org_name)
            
            if not organism:
                # Try alternate path
                org_simple = experiment.get("organism", [])
                if isinstance(org_simple, list):
                    organism = [o.get("scientific_name", str(o)) if isinstance(o, dict) else str(o) for o in org_simple]
            
            # Sample count (number of replicates)
            sample_count = len(replicates) if replicates else 0
            
            # Parse dates (keep as strings)
            submission_date = experiment.get("date_released") or experiment.get("date_created")
            release_date = experiment.get("date_released")
            
            # Detect types
            omics_type = self._detect_omics_type(experiment)
            assay_type = self._detect_assay_type(experiment)
            
            # Compute relevance (uses base class method with ENCODE stopwords)
            relevance = self._compute_keyword_relevance(query, title, description, ENCODE_STOPWORDS) if query else 0.5
            
            # Lab info
            lab = experiment.get("lab", {})
            lab_name = lab.get("title", "") if isinstance(lab, dict) else ""
            
            award = experiment.get("award", {})
            award_name = award.get("project", "") if isinstance(award, dict) else ""
            
            # Create extension
            extension = ENCODEExtension(
                experiment_id=accession,
                target=target_label,
                biosample_term=biosample_term,
                biosample_type=biosample.get("classification", "") if isinstance(biosample, dict) else "",
                lab=lab_name,
                award=award_name,
                status=experiment.get("status", "released"),
                replicates=sample_count,
            )
            
            # Build download links - start with individual files from the experiment
            download_links = []
            
            # Extract actual data files from the experiment
            files = experiment.get("files", [])
            for file_info in files:
                if isinstance(file_info, dict):
                    file_accession = file_info.get("accession", "")
                    file_status = file_info.get("status", "")
                    file_href = file_info.get("href", "")
                    file_size = file_info.get("file_size")
                    file_type_str = file_info.get("file_type", "")
                    output_type = file_info.get("output_type", "")
                    assembly = file_info.get("assembly", "")
                    
                    # Only include released files
                    if file_status != "released" or not file_href:
                        continue
                    
                    # Build download URL
                    download_url = f"https://www.encodeproject.org{file_href}"
                    
                    # Build description
                    desc_parts = [file_accession]
                    if file_type_str:
                        desc_parts.append(f"({file_type_str})")
                    if output_type:
                        desc_parts.append(f"- {output_type}")
                    if assembly:
                        desc_parts.append(f"[{assembly}]")
                    description = " ".join(desc_parts)
                    
                    # Determine file type category
                    file_type_lower = file_type_str.lower()
                    if "fastq" in file_type_lower:
                        dl_file_type = DownloadFileType.RAW
                    elif "bam" in file_type_lower or "sam" in file_type_lower:
                        dl_file_type = DownloadFileType.RAW
                    elif "bed" in file_type_lower or "peak" in file_type_lower:
                        dl_file_type = DownloadFileType.OTHER
                    elif "bigwig" in file_type_lower or "bigbed" in file_type_lower:
                        dl_file_type = DownloadFileType.OTHER
                    else:
                        dl_file_type = DownloadFileType.OTHER
                    
                    download_links.append(DownloadLink(
                        url=download_url,
                        file_type=dl_file_type,
                        description=description,
                        file_size_bytes=file_size,
                        protocol="https",
                    ))
            
            # Add utility links (metadata and batch download)
            download_links.append(DownloadLink(
                url=f"https://www.encodeproject.org/metadata/?type=Experiment&accession={accession}",
                file_type=DownloadFileType.METADATA,
                description="Download metadata TSV",
                protocol="https",
            ))
            download_links.append(DownloadLink(
                url=f"https://www.encodeproject.org/batch_download/?type=Experiment&accession={accession}",
                file_type=DownloadFileType.RAW,
                description="Batch download (all files)",
                protocol="https",
            ))
            
            return UnifiedDataset(
                id=self.build_unified_id(accession),
                accession=accession,
                source=self.source,
                source_url=f"https://www.encodeproject.org/experiments/{accession}/",
                title=title,
                description=description,
                organism=organism,
                sample_count=sample_count,
                sample_count_display=str(sample_count) if sample_count else "N/A",
                submission_date=submission_date,
                release_date=release_date,
                omics_type=omics_type,
                assay_types=[assay_type] if assay_type else [],
                relevance_score=relevance,
                download_links=download_links,
                extensions={"encode": extension.model_dump()},
            )
            
        except Exception as e:
            logger.warning("Failed to parse ENCODE experiment", error=str(e))
            return None
    
    def _create_source_url(self, accession: str) -> str:
        """Create URL to ENCODE experiment page."""
        return f"https://www.encodeproject.org/experiments/{accession}/"


# Singleton instance
_encode_adapter: Optional[ENCODEAdapter] = None


def get_encode_adapter() -> ENCODEAdapter:
    """Get or create ENCODE adapter instance."""
    global _encode_adapter
    if _encode_adapter is None:
        _encode_adapter = ENCODEAdapter()
    return _encode_adapter
