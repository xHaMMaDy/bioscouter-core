"""
CellxGene Adapter
Searches CellxGene Data Portal for single-cell RNA-seq datasets.
API Docs: https://api.cellxgene.cziscience.com/curation/ui/
"""

import asyncio
import re
from typing import List, Optional, Any, Set
from datetime import datetime, timezone
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
    CellxGeneExtension,
    DownloadLink,
    DownloadFileType,
)

logger = structlog.get_logger(__name__)

# CellxGene API base URL (using the Discover API)
CELLXGENE_API_BASE = "https://api.cellxgene.cziscience.com"

# CellxGene-specific stopwords (supplement to COMMON_STOPWORDS)
CELLXGENE_STOPWORDS: Set[str] = frozenset({
    'single', 'cell', 'cells', 'scrna', 'scrnaseq', 'rna', 'seq',
    'sequencing', 'atlas', 'cellxgene', 'czi'
})

# Cache duration for collections (seconds)
COLLECTIONS_CACHE_TTL = 300  # 5 minutes


class CellxGeneError(Exception):
    """CellxGene API error."""
    pass


class CellxGeneAdapter(BaseSourceAdapter):
    """
    Adapter for CellxGene Data Portal.
    Searches for single-cell RNA-seq datasets from the Chan Zuckerberg Initiative.
    """
    
    def __init__(self, timeout: float = 30.0):
        super().__init__()
        self.timeout = timeout
        self._client: Optional[httpx.AsyncClient] = None
    
    @property
    def source(self) -> DataSource:
        return DataSource.CELLXGENE
    
    @property
    def supported_omics(self) -> List[OmicsType]:
        return [OmicsType.SINGLE_CELL, OmicsType.TRANSCRIPTOMICS]
    
    async def _get_client(self) -> httpx.AsyncClient:
        """Get or create HTTP client."""
        if self._client is None or self._client.is_closed:
            self._client = httpx.AsyncClient(
                timeout=self.timeout,
                limits=build_default_limits(),
                headers={
                    "Accept": "application/json",
                    "Content-Type": "application/json",
                    "User-Agent": "BioScouter/1.0",
                }
            )
        return self._client
    
    async def close(self):
        """Close HTTP client."""
        if self._client and not self._client.is_closed:
            await self._client.aclose()
            self._client = None
    
    def _parse_organism(self, organism: Optional[str]) -> Optional[str]:
        """Convert organism name to CellxGene organism ID format."""
        if not organism:
            return None
        
        organism_lower = organism.lower()
        organism_mapping = {
            "human": "NCBITaxon:9606",
            "homo sapiens": "NCBITaxon:9606",
            "mouse": "NCBITaxon:10090",
            "mus musculus": "NCBITaxon:10090",
        }
        
        return organism_mapping.get(organism_lower)
    
    def _parse_tissue(self, query: str) -> Optional[str]:
        """Detect tissue from query for filtering."""
        query_lower = query.lower()
        
        tissue_terms = {
            "brain": "brain",
            "heart": "heart",
            "liver": "liver",
            "lung": "lung",
            "kidney": "kidney",
            "blood": "blood",
            "pbmc": "blood",
            "bone marrow": "bone marrow",
            "skin": "skin",
            "intestine": "intestine",
            "colon": "large intestine",
            "gut": "intestine",
            "pancreas": "pancreas",
            "spleen": "spleen",
            "thymus": "thymus",
            "lymph node": "lymph node",
            "retina": "retina",
            "eye": "eye",
            "muscle": "skeletal muscle tissue",
            "adipose": "adipose tissue",
            "placenta": "placenta",
            "embryo": "embryo",
            "organoid": "organoid",
            "tumor": "tumor",
            "cancer": "tumor",
        }
        
        for term, tissue in tissue_terms.items():
            if term in query_lower:
                return tissue
        
        return None
    
    def _parse_disease(self, query: str) -> Optional[str]:
        """Detect disease from query for filtering."""
        query_lower = query.lower()
        
        disease_terms = {
            "covid": "COVID-19",
            "covid-19": "COVID-19",
            "sars-cov-2": "COVID-19",
            "alzheimer": "Alzheimer disease",
            "parkinson": "Parkinson disease",
            "cancer": "cancer",
            "tumor": "cancer",
            "diabetes": "diabetes mellitus",
            "lupus": "systemic lupus erythematosus",
            "arthritis": "rheumatoid arthritis",
            "asthma": "asthma",
            "fibrosis": "fibrosis",
            "cirrhosis": "liver cirrhosis",
            "autism": "autism spectrum disorder",
            "multiple sclerosis": "multiple sclerosis",
            "ms": "multiple sclerosis",
            "healthy": "normal",
            "normal": "normal",
        }
        
        for term, disease in disease_terms.items():
            if term in query_lower:
                return disease
        
        return None
    
    async def _fetch_with_retry(self, url: str) -> httpx.Response:
        """Fetch URL with retry logic. Delegates to the shared base helper."""
        client = await self._get_client()
        return await super()._fetch_with_retry(client, url)
    
    async def _get_collections_cached(self) -> List[dict]:
        """
        Get collections with caching to avoid repeated large fetches.
        Addresses memory bloat issue by caching collections list.
        """
        import time
        
        # Check if we have a valid cached copy
        if (hasattr(self, '_collections_cache') and 
            hasattr(self, '_collections_cache_time') and
            time.time() - self._collections_cache_time < COLLECTIONS_CACHE_TTL):
            return self._collections_cache
        
        # Fetch fresh collections
        url = f"{CELLXGENE_API_BASE}/curation/v1/collections"
        response = await self._fetch_with_retry(url)
        response.raise_for_status()
        
        collections = response.json()
        
        # Cache the result
        self._collections_cache = collections
        self._collections_cache_time = time.time()
        
        return collections
    
    def _detect_assay_type(self, assay: Optional[str]) -> AssayType:
        """Detect assay type from CellxGene assay field."""
        if not assay:
            return AssayType.OTHER
        
        assay_lower = assay.lower()
        
        if "10x" in assay_lower or "chromium" in assay_lower:
            return AssayType.RNA_SEQ
        if "smart-seq" in assay_lower:
            return AssayType.SMART_SEQ
        if "cite-seq" in assay_lower:
            return AssayType.CITE_SEQ
        if "drop-seq" in assay_lower:
            return AssayType.RNA_SEQ
        if "multiome" in assay_lower:
            return AssayType.MULTIOME
        if "atac" in assay_lower:
            return AssayType.ATAC_SEQ
        if "rna" in assay_lower or "scrna" in assay_lower:
            return AssayType.RNA_SEQ
        
        return AssayType.OTHER
    
    async def search(
        self,
        query: str,
        max_results: int = 50,
        organism: Optional[str] = None,
        **kwargs
    ) -> List[UnifiedDataset]:
        """
        Search CellxGene for single-cell datasets.
        Uses the collections API with caching.
        
        Args:
            query: Search query
            max_results: Maximum results to return
            organism: Optional organism filter
            
        Returns:
            List of UnifiedDataset objects
            
        Raises:
            CellxGeneError: On API or network failures
        """
        try:
            logger.info("Searching CellxGene", query=query, max_results=max_results)
            
            # Get collections with caching (addresses memory bloat issue)
            collections = await self._get_collections_cached()
            
            logger.info("CellxGene response", collection_count=len(collections))
            
            # Parse tissue and disease from query for filtering
            target_tissue = self._parse_tissue(query)
            target_disease = self._parse_disease(query)
            target_organism = self._parse_organism(organism)
            
            # Filter and convert collections
            unified = []
            for collection in collections:
                try:
                    # Filter by visibility
                    if collection.get("visibility") != "PUBLIC":
                        continue
                    
                    ud = self._to_unified(
                        collection, 
                        query=query,
                        target_tissue=target_tissue,
                        target_disease=target_disease,
                        target_organism=target_organism
                    )
                    if ud:
                        unified.append(ud)
                except Exception as e:
                    logger.warning("Failed to parse CellxGene collection", error=str(e))
                    continue
            
            # Sort by relevance score
            unified.sort(key=lambda x: x.relevance_score, reverse=True)
            
            # Return top results
            unified = unified[:max_results]
            
            logger.info("CellxGene search complete", results=len(unified))
            return unified
            
        except httpx.HTTPStatusError as e:
            logger.error("CellxGene API error", status=e.response.status_code, error=str(e))
            raise CellxGeneError(f"CellxGene API error: {e.response.status_code}") from e
        except httpx.RequestError as e:
            logger.error("CellxGene request failed", error=str(e))
            raise CellxGeneError(f"CellxGene request failed: {str(e)}") from e
    
    async def get_dataset(self, accession: str) -> Optional[UnifiedDataset]:
        """Fetch a single CellxGene collection by ID."""
        try:
            url = f"{CELLXGENE_API_BASE}/curation/v1/collections/{accession}"
            
            client = await self._get_client()
            response = await client.get(url)
            
            if response.status_code == 404:
                return None
            
            response.raise_for_status()
            collection = response.json()
            
            return self._to_unified(collection)
            
        except httpx.HTTPStatusError as e:
            if e.response.status_code == 404:
                return None
            logger.error("CellxGene API error", accession=accession, error=str(e))
            return None
        except Exception as e:
            logger.error("Failed to fetch CellxGene collection", accession=accession, error=str(e))
            return None
    
    def _to_unified(
        self, 
        collection: dict, 
        query: Optional[str] = None,
        target_tissue: Optional[str] = None,
        target_disease: Optional[str] = None,
        target_organism: Optional[str] = None
    ) -> Optional[UnifiedDataset]:
        """Convert a CellxGene collection to UnifiedDataset."""
        try:
            collection_id = collection.get("collection_id") or collection.get("id") or ""
            if not collection_id:
                return None
            
            # Extract basic info
            name = collection.get("name") or ""
            description = collection.get("description") or ""
            
            # Get datasets in collection to extract metadata
            datasets = collection.get("datasets", [])
            
            # Aggregate metadata from datasets
            organisms = set()
            tissues = set()
            diseases = set()
            assays = set()
            total_cells = 0
            cell_types = set()
            download_links = []
            
            for dataset in datasets:
                # Organisms
                for org in dataset.get("organism", []):
                    if isinstance(org, dict):
                        organisms.add(org.get("label", ""))
                    elif isinstance(org, str):
                        organisms.add(org)
                
                # Tissues
                for tissue in dataset.get("tissue", []):
                    if isinstance(tissue, dict):
                        tissues.add(tissue.get("label", ""))
                    elif isinstance(tissue, str):
                        tissues.add(tissue)
                
                # Diseases
                for disease in dataset.get("disease", []):
                    if isinstance(disease, dict):
                        diseases.add(disease.get("label", ""))
                    elif isinstance(disease, str):
                        diseases.add(disease)
                
                # Assays
                for assay in dataset.get("assay", []):
                    if isinstance(assay, dict):
                        assays.add(assay.get("label", ""))
                    elif isinstance(assay, str):
                        assays.add(assay)
                
                # Cell count
                total_cells += dataset.get("cell_count", 0)
                
                # Cell types
                for ct in dataset.get("cell_type", []):
                    if isinstance(ct, dict):
                        cell_types.add(ct.get("label", ""))
                    elif isinstance(ct, str):
                        cell_types.add(ct)
                
                # Download links from assets
                dataset_title = dataset.get("title", "")
                for asset in dataset.get("assets", []):
                    url = asset.get("url")
                    if url:
                        file_type = asset.get("filetype", "").lower()
                        file_size = asset.get("filesize")
                        
                        # Determine file type enum
                        if file_type == "h5ad":
                            dl_type = DownloadFileType.MATRIX
                            description = "AnnData format (.h5ad) - gene expression matrix with metadata"
                        elif file_type == "rds":
                            dl_type = DownloadFileType.MATRIX
                            description = "Seurat/R format (.rds) - gene expression matrix"
                        else:
                            dl_type = DownloadFileType.PROCESSED
                            description = f"{file_type.upper()} file"
                        
                        # Extract filename from URL
                        file_name = url.split("/")[-1] if "/" in url else None
                        
                        download_links.append(DownloadLink(
                            url=url,
                            file_type=dl_type,
                            file_name=file_name,
                            file_size_bytes=file_size,
                            protocol="https",
                            description=f"{dataset_title} - {description}" if dataset_title else description,
                        ))

            # Fallback to collection-level count if datasets are missing counts
            if not total_cells:
                total_cells = collection.get("cell_count", 0)
            
            # Clean up sets
            organisms = {o for o in organisms if o}
            tissues = {t for t in tissues if t}
            diseases = {d for d in diseases if d}
            assays = {a for a in assays if a}
            cell_types = {c for c in cell_types if c}
            
            # Filter by target criteria
            match_score = 0.0
            
            if target_tissue:
                if any(target_tissue.lower() in t.lower() for t in tissues):
                    match_score += 0.3
                elif not tissues:
                    pass  # Don't penalize if no tissue info
                else:
                    match_score -= 0.1
            
            if target_disease:
                if any(target_disease.lower() in d.lower() for d in diseases):
                    match_score += 0.3
                elif "normal" in diseases or "healthy" in str(diseases).lower():
                    if target_disease.lower() in ["normal", "healthy"]:
                        match_score += 0.2
            
            if target_organism:
                org_found = False
                for org in organisms:
                    if "human" in org.lower() or "sapiens" in org.lower():
                        if "9606" in target_organism:
                            org_found = True
                            match_score += 0.2
                            break
                    elif "mouse" in org.lower() or "musculus" in org.lower():
                        if "10090" in target_organism:
                            org_found = True
                            match_score += 0.2
                            break
            
            # Compute base relevance (uses base class method with CellxGene stopwords)
            base_relevance = self._compute_keyword_relevance(query, name, description, CELLXGENE_STOPWORDS) if query else 0.5
            relevance = max(0.3, min(1.0, base_relevance + match_score))
            
            # Detect assay type
            primary_assay = list(assays)[0] if assays else None
            assay_type = self._detect_assay_type(primary_assay)
            
            # Build summary if description is short
            summary = description
            if len(summary) < 50 and tissues:
                summary = f"{name}. Tissues: {', '.join(list(tissues)[:3])}."
            
            # Publisher metadata
            publisher = collection.get("publisher_metadata") or {}
            pub_date_raw = publisher.get("published_at") or collection.get("created_at")
            
            # Convert release_date to string (CellxGene returns Unix timestamp as float)
            pub_date = None
            if pub_date_raw:
                if isinstance(pub_date_raw, (int, float)):
                    # Convert Unix timestamp to ISO format string
                    pub_date = datetime.fromtimestamp(pub_date_raw, tz=timezone.utc).strftime("%Y-%m-%d")
                elif isinstance(pub_date_raw, str):
                    pub_date = pub_date_raw
            
            # DOI links
            links = collection.get("links", [])
            doi = None
            for link in links:
                if link.get("link_type") == "DOI":
                    doi = link.get("link_url")
                    break
            
            # Create extension
            extension = CellxGeneExtension(
                collection_id=collection_id,
                collection_name=name,
                tissue=list(tissues)[0] if tissues else None,
                disease=list(diseases)[0] if diseases else None,
                assay=primary_assay,
                cell_type_count=len(cell_types),
                cell_count=total_cells,
                is_primary_data=collection.get("is_primary_data", True),
                publisher_metadata=publisher if publisher else None,
            )
            
            return UnifiedDataset(
                id=self.build_unified_id(collection_id),
                accession=collection_id,
                source=self.source,
                source_url=f"https://cellxgene.cziscience.com/collections/{collection_id}",
                title=name,
                description=summary,
                organism=list(organisms) if organisms else ["Homo sapiens"],
                sample_count=total_cells,
                sample_count_display=str(total_cells) if total_cells else "N/A",
                omics_type=OmicsType.SINGLE_CELL,
                assay_types=[assay_type] if assay_type else [],
                release_date=pub_date,
                relevance_score=relevance,
                download_links=download_links,
                extensions={"cellxgene": extension.model_dump()},
            )
            
        except Exception as e:
            logger.warning("Failed to parse CellxGene collection", error=str(e))
            return None
    
    def _create_source_url(self, accession: str) -> str:
        """Create URL to CellxGene collection page."""
        return f"https://cellxgene.cziscience.com/collections/{accession}"


# Singleton instance
_cellxgene_adapter: Optional[CellxGeneAdapter] = None


def get_cellxgene_adapter() -> CellxGeneAdapter:
    """Get or create CellxGene adapter instance."""
    global _cellxgene_adapter
    if _cellxgene_adapter is None:
        _cellxgene_adapter = CellxGeneAdapter()
    return _cellxgene_adapter
