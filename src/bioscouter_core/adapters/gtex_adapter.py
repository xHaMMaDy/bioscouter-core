"""
GTEx (Genotype-Tissue Expression) Adapter

Provides access to GTEx tissue-level expression data through the GTEx Portal API v2.
Uses tissue-as-dataset approach where each of 54 tissues is returned as a searchable dataset.

Features:
- Hierarchical tissue organization (e.g., Brain → Brain - Cortex)
- Gene search with TPM sorting (find tissues where gene is highly expressed)
- Multi-gene search with expression heatmap
- eQTL data integration
- Support for GTEx v8 and v10 datasets
- Tissue metadata caching (forever) and expression data caching (24h)

API Documentation: https://gtexportal.org/api/v2/
"""

import asyncio
import hashlib
from dataclasses import dataclass, field
from datetime import datetime, date, timedelta
from typing import Any, Dict, List, Optional, Tuple
from enum import Enum

import httpx
import structlog

from .base import build_default_limits
from bioscouter_core.models.unified import (
    AssayType,
    CurationLevel,
    DataSource,
    OmicsType,
    UnifiedDataset,
    normalize_organism,
)

logger = structlog.get_logger(__name__)

# API Configuration
GTEX_API_BASE_URL = "https://gtexportal.org/api/v2"
GTEX_PORTAL_URL = "https://gtexportal.org/home"

# Cache settings
TISSUE_CACHE_TTL = float('inf')  # Never expires (static data)
EXPRESSION_CACHE_TTL = 86400  # 24 hours

# Rate limiting
RATE_LIMIT_DELAY = 0.2  # 200ms between requests


class GTExVersion(str, Enum):
    """GTEx dataset versions."""
    V8 = "gtex_v8"
    V10 = "gtex_v10"


class GTExSearchError(Exception):
    """Custom exception for GTEx search errors."""
    pass


@dataclass
class GTExExtension:
    """GTEx-specific metadata extension."""
    tissue_site: str  # Parent tissue group (e.g., "Brain")
    tissue_site_detail: str  # Specific tissue (e.g., "Brain - Cortex")
    tissue_site_detail_id: str  # API ID
    sample_count: int = 0
    subject_count: int = 0
    gtex_version: str = "v10"
    
    # Expression data (optional, for gene searches)
    top_genes: List[Dict[str, Any]] = field(default_factory=list)  # [{gene, tpm, rank}]
    query_gene_tpm: Optional[float] = None  # TPM for searched gene
    query_gene_rank: Optional[int] = None  # Rank of searched gene in tissue
    
    # eQTL data
    eqtl_count: int = 0
    eqtl_genes: List[str] = field(default_factory=list)  # Genes with eQTLs in this tissue
    
    # Tissue hierarchy info
    tissue_group: Optional[str] = None
    tissue_color: Optional[str] = None  # Hex color for visualization
    
    def model_dump(self) -> dict:
        return {
            "tissue_site": self.tissue_site,
            "tissue_site_detail": self.tissue_site_detail,
            "tissue_site_detail_id": self.tissue_site_detail_id,
            "sample_count": self.sample_count,
            "subject_count": self.subject_count,
            "gtex_version": self.gtex_version,
            "top_genes": self.top_genes,
            "query_gene_tpm": self.query_gene_tpm,
            "query_gene_rank": self.query_gene_rank,
            "eqtl_count": self.eqtl_count,
            "eqtl_genes": self.eqtl_genes,
            "tissue_group": self.tissue_group,
            "tissue_color": self.tissue_color,
        }


# Tissue hierarchy for grouping
TISSUE_HIERARCHY = {
    "Adipose Tissue": ["Adipose - Subcutaneous", "Adipose - Visceral (Omentum)"],
    "Adrenal Gland": ["Adrenal Gland"],
    "Artery": ["Artery - Aorta", "Artery - Coronary", "Artery - Tibial"],
    "Blood": ["Cells - EBV-transformed lymphocytes", "Cells - Cultured fibroblasts", "Whole Blood"],
    "Brain": [
        "Brain - Amygdala", "Brain - Anterior cingulate cortex (BA24)",
        "Brain - Caudate (basal ganglia)", "Brain - Cerebellar Hemisphere",
        "Brain - Cerebellum", "Brain - Cortex", "Brain - Frontal Cortex (BA9)",
        "Brain - Hippocampus", "Brain - Hypothalamus", "Brain - Nucleus accumbens (basal ganglia)",
        "Brain - Putamen (basal ganglia)", "Brain - Spinal cord (cervical c-1)",
        "Brain - Substantia nigra"
    ],
    "Breast": ["Breast - Mammary Tissue"],
    "Colon": ["Colon - Sigmoid", "Colon - Transverse"],
    "Esophagus": ["Esophagus - Gastroesophageal Junction", "Esophagus - Mucosa", "Esophagus - Muscularis"],
    "Heart": ["Heart - Atrial Appendage", "Heart - Left Ventricle"],
    "Kidney": ["Kidney - Cortex", "Kidney - Medulla"],
    "Liver": ["Liver"],
    "Lung": ["Lung"],
    "Muscle": ["Muscle - Skeletal"],
    "Nerve": ["Nerve - Tibial"],
    "Ovary": ["Ovary"],
    "Pancreas": ["Pancreas"],
    "Pituitary": ["Pituitary"],
    "Prostate": ["Prostate"],
    "Skin": ["Skin - Not Sun Exposed (Suprapubic)", "Skin - Sun Exposed (Lower leg)"],
    "Small Intestine": ["Small Intestine - Terminal Ileum"],
    "Spleen": ["Spleen"],
    "Stomach": ["Stomach"],
    "Testis": ["Testis"],
    "Thyroid": ["Thyroid"],
    "Uterus": ["Uterus"],
    "Vagina": ["Vagina"],
}

# Tissue colors for visualization (matching GTEx portal)
TISSUE_COLORS = {
    "Adipose Tissue": "#FF6600",
    "Adrenal Gland": "#33DD33",
    "Artery": "#FF5555",
    "Blood": "#FF00BB",
    "Brain": "#EEEE00",
    "Breast": "#33FFFF",
    "Colon": "#DDAA77",
    "Esophagus": "#8B7355",
    "Heart": "#9900FF",
    "Kidney": "#227777",
    "Liver": "#AAFF00",
    "Lung": "#99FF99",
    "Muscle": "#AAAAFF",
    "Nerve": "#FFD700",
    "Ovary": "#FFAAFF",
    "Pancreas": "#995522",
    "Pituitary": "#AAFF66",
    "Prostate": "#DDDDDD",
    "Skin": "#0000FF",
    "Small Intestine": "#DD7799",
    "Spleen": "#778855",
    "Stomach": "#FFDD99",
    "Testis": "#AAAAAA",
    "Thyroid": "#006600",
    "Uterus": "#FF66AA",
    "Vagina": "#FF9999",
}


class GTExAdapter:
    """
    Adapter for GTEx Portal API.
    
    Implements tissue-as-dataset search where each tissue is returned as a UnifiedDataset.
    Supports text search on tissue names and gene-specific expression queries.
    """
    
    # Source identifier
    source = DataSource.GTEX
    source_name = "GTEx"
    
    def __init__(self):
        self._client: Optional[httpx.AsyncClient] = None
        self._last_request_time: float = 0
        
        # Caches
        self._tissues_cache: Optional[List[dict]] = None
        self._tissues_cache_time: float = 0
        self._expression_cache: Dict[str, Tuple[dict, float]] = {}  # key -> (data, timestamp)
        
    async def _get_client(self) -> httpx.AsyncClient:
        """Get or create HTTP client."""
        if self._client is None or self._client.is_closed:
            self._client = httpx.AsyncClient(
                timeout=httpx.Timeout(30.0),
                limits=build_default_limits(),
                follow_redirects=True,
                headers={
                    "Accept": "application/json",
                    "User-Agent": "BioScouter/1.0 (GTEx Adapter)",
                }
            )
        return self._client
    
    async def _rate_limit(self):
        """Enforce rate limiting."""
        import time
        now = time.time()
        elapsed = now - self._last_request_time
        if elapsed < RATE_LIMIT_DELAY:
            await asyncio.sleep(RATE_LIMIT_DELAY - elapsed)
        self._last_request_time = time.time()
    
    async def _fetch_tissues(self, version: str = "v10") -> List[dict]:
        """
        Fetch all tissue metadata from GTEx API.
        Caches forever since tissue list is static.
        """
        import time
        now = time.time()
        
        # Check cache
        if self._tissues_cache is not None:
            return self._tissues_cache
        
        await self._rate_limit()
        client = await self._get_client()
        
        dataset_id = GTExVersion.V10.value if version == "v10" else GTExVersion.V8.value
        
        try:
            logger.info("Fetching GTEx tissue metadata", version=version)
            
            response = await client.get(
                f"{GTEX_API_BASE_URL}/dataset/tissueSiteDetail",
                params={"datasetId": dataset_id}
            )
            response.raise_for_status()
            
            data = response.json()
            tissues = data.get("data", []) if isinstance(data, dict) else data
            
            # Enrich with hierarchy info
            for tissue in tissues:
                tissue_name = tissue.get("tissueSiteDetail", "")
                for group, members in TISSUE_HIERARCHY.items():
                    if tissue_name in members:
                        tissue["tissueGroup"] = group
                        tissue["tissueColor"] = TISSUE_COLORS.get(group, "#888888")
                        break
            
            self._tissues_cache = tissues
            self._tissues_cache_time = now
            
            logger.info("GTEx tissues cached", count=len(tissues))
            
            return tissues
            
        except Exception as e:
            logger.error("Failed to fetch GTEx tissues", error=str(e))
            # Return cached data if available, even if stale
            if self._tissues_cache:
                logger.warning("Using stale GTEx tissue cache")
                return self._tissues_cache
            raise GTExSearchError(f"Failed to fetch GTEx tissues: {str(e)}")
    
    async def _fetch_gene_expression(
        self,
        gene_symbol: str,
        version: str = "v10"
    ) -> Dict[str, dict]:
        """
        Fetch expression data for a gene across all tissues.
        Returns dict mapping tissue_id -> {median_tpm, ...}
        """
        import time
        cache_key = f"{gene_symbol}_{version}"
        
        # Check cache
        if cache_key in self._expression_cache:
            data, cached_time = self._expression_cache[cache_key]
            if time.time() - cached_time < EXPRESSION_CACHE_TTL:
                return data
        
        await self._rate_limit()
        client = await self._get_client()
        
        dataset_id = GTExVersion.V10.value if version == "v10" else GTExVersion.V8.value
        
        try:
            logger.info("Fetching GTEx expression", gene=gene_symbol, version=version)
            
            # First, search for gene to get gencodeId
            gene_response = await client.get(
                f"{GTEX_API_BASE_URL}/reference/gene",
                params={"geneId": gene_symbol, "datasetId": dataset_id}
            )
            
            if gene_response.status_code != 200:
                return {}
            
            gene_data = gene_response.json()
            genes = gene_data.get("data", [])
            
            if not genes:
                return {}
            
            gencode_id = genes[0].get("gencodeId")
            if not gencode_id:
                return {}
            
            # Fetch expression data
            await self._rate_limit()
            expr_response = await client.get(
                f"{GTEX_API_BASE_URL}/expression/geneExpression",
                params={
                    "datasetId": dataset_id,
                    "gencodeId": gencode_id,
                }
            )
            
            if expr_response.status_code != 200:
                return {}
            
            expr_data = expr_response.json()
            expressions = expr_data.get("data", [])
            
            # Map by tissue
            result = {}
            for expr in expressions:
                tissue_id = expr.get("tissueSiteDetailId")
                if tissue_id:
                    result[tissue_id] = {
                        "median_tpm": expr.get("median", 0),
                        "gene_symbol": gene_symbol,
                        "gencode_id": gencode_id,
                    }
            
            # Cache
            self._expression_cache[cache_key] = (result, time.time())
            
            return result
            
        except Exception as e:
            logger.error("Failed to fetch GTEx expression", gene=gene_symbol, error=str(e))
            return {}
    
    async def _fetch_tissue_eqtls(
        self,
        tissue_id: str,
        version: str = "v10"
    ) -> dict:
        """Fetch eQTL summary for a tissue."""
        # Note: eQTL queries can be slow; consider caching
        await self._rate_limit()
        client = await self._get_client()
        
        dataset_id = GTExVersion.V10.value if version == "v10" else GTExVersion.V8.value
        
        try:
            response = await client.get(
                f"{GTEX_API_BASE_URL}/association/singleTissueEqtl",
                params={
                    "datasetId": dataset_id,
                    "tissueSiteDetailId": tissue_id,
                    "page": 0,
                    "itemsPerPage": 1,  # Just to get count
                }
            )
            
            if response.status_code == 200:
                data = response.json()
                return {
                    "count": data.get("numResults", 0),
                }
            return {"count": 0}
            
        except Exception:
            return {"count": 0}
    
    def _tissue_matches_query(self, tissue: dict, query_terms: List[str]) -> bool:
        """Check if tissue matches search query."""
        searchable = " ".join([
            tissue.get("tissueSite", "") or "",
            tissue.get("tissueSiteDetail", "") or "",
            tissue.get("tissueGroup", "") or "",
        ]).lower()
        
        # Any term matches
        return any(term in searchable for term in query_terms)
    
    def _tissue_to_dataset(
        self,
        tissue: dict,
        version: str = "v10",
        gene_expression: Optional[dict] = None,
        include_eqtl: bool = True,
        eqtl_data: Optional[dict] = None,
    ) -> UnifiedDataset:
        """Convert GTEx tissue metadata to UnifiedDataset."""
        tissue_id = tissue.get("tissueSiteDetailId", "")
        tissue_name = tissue.get("tissueSiteDetail", "Unknown Tissue")
        tissue_site = tissue.get("tissueSite", "")
        
        # Build description
        sample_count = tissue.get("rnaSeqAndGenotypeSampleCount", 0) or tissue.get("rnaSeqSampleCount", 0) or 0
        subject_count = tissue.get("rnaSeqAndGenotypeSubjectCount", 0) or tissue.get("rnaSeqSubjectCount", 0) or 0
        
        description = f"GTEx {version.upper()} expression data for {tissue_name}. "
        description += f"Contains RNA-seq data from {sample_count} samples across {subject_count} subjects."
        
        if gene_expression and gene_expression.get("median_tpm"):
            gene = gene_expression.get("gene_symbol", "")
            tpm = gene_expression.get("median_tpm", 0)
            description += f" {gene} median TPM: {tpm:.2f}."
        
        # Source URL
        source_url = f"{GTEX_PORTAL_URL}/tissue/{tissue_id}"
        
        # Build extension
        extension = GTExExtension(
            tissue_site=tissue_site,
            tissue_site_detail=tissue_name,
            tissue_site_detail_id=tissue_id,
            sample_count=sample_count,
            subject_count=subject_count,
            gtex_version=version,
            tissue_group=tissue.get("tissueGroup"),
            tissue_color=tissue.get("tissueColor"),
            eqtl_count=eqtl_data.get("count", 0) if eqtl_data else 0,
        )
        
        # Add gene expression data if available
        if gene_expression:
            extension.query_gene_tpm = gene_expression.get("median_tpm")
            extension.top_genes = [{
                "gene": gene_expression.get("gene_symbol"),
                "tpm": gene_expression.get("median_tpm"),
            }]
        
        return UnifiedDataset(
            id=f"gtex:{tissue_id}:{version}",
            accession=tissue_id,
            source=DataSource.GTEX,
            source_url=source_url,
            omics_type=OmicsType.TRANSCRIPTOMICS,
            assay_types=[AssayType.RNA_SEQ],
            title=f"GTEx {tissue_name}",
            description=description,
            organism=["Homo sapiens"],
            sample_count=sample_count,
            sample_count_display=f"{sample_count} samples, {subject_count} subjects",
            tissue=[tissue_site, tissue_name] if tissue_site != tissue_name else [tissue_name],
            curation_level=CurationLevel.CURATED,
            extensions={"gtex": extension.model_dump()},
        )
    
    async def search(
        self,
        query: str,
        max_results: int = 54,
        organism: Optional[str] = None,  # Ignored - GTEx is human only
        version: str = "v10",
        include_eqtl: bool = True,
    ) -> List[UnifiedDataset]:
        """
        Search GTEx for tissues matching query.
        
        Supports:
        - Tissue name search: "brain", "heart", "liver"
        - Gene expression search: "BRCA1" (returns tissues sorted by TPM)
        - Combined: "BRCA1 brain" (gene in brain tissues)
        
        Args:
            query: Search query (tissue names, gene symbols, or both)
            max_results: Maximum tissues to return (max 54)
            organism: Ignored (GTEx is human-only)
            version: GTEx version (v8, v10)
            include_eqtl: Include eQTL counts
            
        Returns:
            List of UnifiedDataset, one per matching tissue
        """
        logger.info("Searching GTEx", query=query, max_results=max_results, version=version)
        
        # Parse query for potential gene symbols (uppercase, 2-10 chars)
        query_terms = [t.lower() for t in query.split() if len(t) > 1]
        potential_genes = [t.upper() for t in query.split() if t.isupper() or (len(t) >= 2 and t[0].isupper())]
        
        try:
            # Fetch all tissues
            tissues = await self._fetch_tissues(version)
            
            if not tissues:
                logger.warning("No GTEx tissues found")
                return []
            
            # Check if any query term is a gene symbol
            gene_expression_map: Dict[str, dict] = {}
            gene_found = None
            
            for gene_candidate in potential_genes[:3]:  # Check up to 3 potential genes
                expr_data = await self._fetch_gene_expression(gene_candidate, version)
                if expr_data:
                    gene_expression_map = expr_data
                    gene_found = gene_candidate
                    logger.info("GTEx gene found", gene=gene_found, tissues_with_expression=len(expr_data))
                    break
            
            # Filter tissues
            matching_tissues = []
            
            for tissue in tissues:
                tissue_id = tissue.get("tissueSiteDetailId", "")
                
                # If gene search, only include tissues with expression
                if gene_found:
                    if tissue_id not in gene_expression_map:
                        continue
                
                # Text filter on tissue names (if not pure gene search)
                tissue_terms = [t for t in query_terms if t != gene_found.lower()] if gene_found else query_terms
                
                if tissue_terms and not self._tissue_matches_query(tissue, tissue_terms):
                    continue
                
                # Include tissue
                gene_expr = gene_expression_map.get(tissue_id) if gene_found else None
                matching_tissues.append((tissue, gene_expr))
            
            # Sort by expression if gene search, otherwise by sample count
            if gene_found:
                matching_tissues.sort(
                    key=lambda x: x[1].get("median_tpm", 0) if x[1] else 0,
                    reverse=True
                )
            else:
                matching_tissues.sort(
                    key=lambda x: x[0].get("rnaSeqAndGenotypeSampleCount", 0) or 0,
                    reverse=True
                )
            
            # Limit results
            matching_tissues = matching_tissues[:max_results]
            
            # Convert to datasets
            datasets = []
            for tissue, gene_expr in matching_tissues:
                # Optionally fetch eQTL data (can be slow)
                eqtl_data = None
                if include_eqtl:
                    # Skip eQTL fetch for speed; use cached count if available
                    eqtl_data = {"count": tissue.get("eqtlCount", 0)}
                
                dataset = self._tissue_to_dataset(
                    tissue,
                    version=version,
                    gene_expression=gene_expr,
                    include_eqtl=include_eqtl,
                    eqtl_data=eqtl_data,
                )
                
                # Calculate relevance score
                if gene_found and gene_expr:
                    # Gene search: higher TPM = higher relevance
                    tpm = gene_expr.get("median_tpm", 0)
                    max_tpm = max(
                        (e.get("median_tpm", 0) for e in gene_expression_map.values()),
                        default=1
                    )
                    dataset.relevance_score = min(tpm / max_tpm, 1.0) if max_tpm > 0 else 0.5
                else:
                    # Text search: base relevance
                    dataset.relevance_score = 0.5
                
                datasets.append(dataset)
            
            logger.info("GTEx search complete", results=len(datasets))
            
            return datasets
            
        except Exception as e:
            logger.error("GTEx search failed", error=str(e))
            raise GTExSearchError(f"GTEx search failed: {str(e)}")
    
    async def get_tissue_hierarchy(self, version: str = "v10") -> Dict[str, List[dict]]:
        """
        Get tissues organized by hierarchy for UI display.
        
        Returns:
            Dict mapping tissue group -> list of tissue details
        """
        tissues = await self._fetch_tissues(version)
        
        hierarchy = {}
        for tissue in tissues:
            group = tissue.get("tissueGroup") or tissue.get("tissueSite", "Other")
            if group not in hierarchy:
                hierarchy[group] = []
            hierarchy[group].append({
                "id": tissue.get("tissueSiteDetailId"),
                "name": tissue.get("tissueSiteDetail"),
                "sample_count": tissue.get("rnaSeqAndGenotypeSampleCount", 0),
                "color": tissue.get("tissueColor", "#888888"),
            })
        
        # Sort tissues within each group by sample count
        for group in hierarchy:
            hierarchy[group].sort(key=lambda x: x["sample_count"], reverse=True)
        
        return hierarchy
    
    async def compare_tissues(
        self,
        tissue_ids: List[str],
        genes: List[str],
        version: str = "v10"
    ) -> Dict[str, Dict[str, float]]:
        """
        Get expression matrix for comparing tissues and genes.
        
        Args:
            tissue_ids: List of tissue IDs to compare
            genes: List of gene symbols
            version: GTEx version
            
        Returns:
            Dict mapping gene -> {tissue_id: tpm}
        """
        result = {}
        
        for gene in genes[:10]:  # Limit to 10 genes
            expr_data = await self._fetch_gene_expression(gene, version)
            result[gene] = {
                tid: expr_data.get(tid, {}).get("median_tpm", 0)
                for tid in tissue_ids
                if tid in expr_data
            }
        
        return result
    
    async def get_dataset(self, accession: str) -> Optional[UnifiedDataset]:
        """
        Fetch a single GTEx tissue dataset by accession.
        
        Args:
            accession: GTEx tissue ID (e.g., "Brain_Cortex", "Muscle_Skeletal")
            
        Returns:
            UnifiedDataset if found, None otherwise
        """
        try:
            logger.info("Fetching GTEx dataset", accession=accession)
            
            # Fetch tissue data
            tissues = await self._fetch_tissues()
            
            if not tissues:
                logger.warning("Could not fetch GTEx tissue data")
                return None
            
            # Find matching tissue (case-insensitive)
            accession_lower = accession.lower().replace("_", " ").replace("-", " ")
            
            matching_tissue = None
            for tissue in tissues:
                tissue_id = tissue.get("tissueSiteDetailId", "")
                tissue_name = tissue.get("tissueSiteDetail", "")
                
                if (tissue_id.lower() == accession.lower() or
                    tissue_id.lower().replace("_", " ") == accession_lower or
                    tissue_name.lower() == accession_lower):
                    matching_tissue = tissue
                    break
            
            if not matching_tissue:
                logger.warning("GTEx tissue not found", accession=accession)
                return None
            
            # Convert to unified dataset
            dataset = self._tissue_to_dataset(matching_tissue)
            
            logger.info("GTEx dataset fetched", accession=accession)
            return dataset
            
        except Exception as e:
            logger.error("Failed to fetch GTEx dataset", 
                        accession=accession, error=str(e))
            return None
    
    async def close(self):
        """Close HTTP client."""
        if self._client and not self._client.is_closed:
            await self._client.aclose()


# Singleton instance
_adapter_instance: Optional[GTExAdapter] = None


def get_gtex_adapter() -> GTExAdapter:
    """Get singleton GTEx adapter instance."""
    global _adapter_instance
    if _adapter_instance is None:
        _adapter_instance = GTExAdapter()
    return _adapter_instance
