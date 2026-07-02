"""
GDC (Genomic Data Commons) Adapter
Searches the GDC Data Portal for cancer genomics datasets (TCGA, TARGET, etc.).
API Docs: https://docs.gdc.cancer.gov/API/Users_Guide/
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
    GDCExtension,
    DownloadLink,
    DownloadFileType,
)

logger = structlog.get_logger(__name__)

# GDC API base URL
GDC_API_BASE = "https://api.gdc.cancer.gov"

# GDC-specific stopwords (supplement to COMMON_STOPWORDS)
GDC_STOPWORDS: Set[str] = frozenset({
    'cancer', 'tcga', 'gdc', 'genomic', 'genomics', 'tumor', 'tumour',
    'target', 'cptac', 'project', 'program'
})


class GDCError(Exception):
    """GDC API error."""
    pass


class GDCAdapter(BaseSourceAdapter):
    """
    Adapter for GDC Data Portal (TCGA, TARGET, CPTAC).
    Uses GDC REST API.
    """
    
    def __init__(self, timeout: float = 30.0):
        super().__init__()
        self.timeout = timeout
        self._client: Optional[httpx.AsyncClient] = None
    
    @property
    def source(self) -> DataSource:
        return DataSource.GDC
    
    @property
    def supported_omics(self) -> List[OmicsType]:
        return [OmicsType.MULTI_OMICS, OmicsType.GENOMICS, OmicsType.TRANSCRIPTOMICS]
    
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
    
    def _build_filters(self, query: str, organism: Optional[str] = None) -> Optional[dict]:
        """Build GDC API filters from query."""
        filters = {
            "op": "and",
            "content": []
        }
        
        # Note: GDC API doesn't support 'state' filter - projects are returned based on availability
        
        query_lower = query.lower()
        
        # Detect cancer types
        cancer_mapping = {
            "breast": "BRCA",
            "lung": "LUAD",
            "lung adenocarcinoma": "LUAD",
            "lung squamous": "LUSC",
            "colon": "COAD",
            "colorectal": "COAD",
            "prostate": "PRAD",
            "ovarian": "OV",
            "ovary": "OV",
            "kidney": "KIRC",
            "renal": "KIRC",
            "liver": "LIHC",
            "hepatocellular": "LIHC",
            "stomach": "STAD",
            "gastric": "STAD",
            "bladder": "BLCA",
            "melanoma": "SKCM",
            "skin": "SKCM",
            "thyroid": "THCA",
            "pancreatic": "PAAD",
            "pancreas": "PAAD",
            "glioblastoma": "GBM",
            "glioma": "LGG",
            "brain": "GBM",
            "leukemia": "LAML",
            "aml": "LAML",
            "cervical": "CESC",
            "head and neck": "HNSC",
            "esophageal": "ESCA",
            "sarcoma": "SARC",
            "testicular": "TGCT",
            "uterine": "UCEC",
            "endometrial": "UCEC",
        }
        
        detected_project = None
        for term, project in cancer_mapping.items():
            if term in query_lower:
                detected_project = f"TCGA-{project}"
                break
        
        if detected_project:
            filters["content"].append({
                "op": "=",
                "content": {
                    "field": "project_id",
                    "value": detected_project
                }
            })
        
        # Detect data types
        data_type_mapping = {
            "rna-seq": "Gene Expression Quantification",
            "rnaseq": "Gene Expression Quantification",
            "transcriptome": "Gene Expression Quantification",
            "gene expression": "Gene Expression Quantification",
            "mutation": "Masked Somatic Mutation",
            "somatic": "Masked Somatic Mutation",
            "variant": "Masked Somatic Mutation",
            "snv": "Masked Somatic Mutation",
            "copy number": "Copy Number Segment",
            "cnv": "Copy Number Segment",
            "methylation": "Methylation Beta Value",
            "mirna": "miRNA Expression Quantification",
            "microrna": "miRNA Expression Quantification",
            "proteomics": "Protein Expression Quantification",
            "protein": "Protein Expression Quantification",
        }
        
        for term, data_type in data_type_mapping.items():
            if term in query_lower:
                filters["content"].append({
                    "op": "=",
                    "content": {
                        "field": "data_type",
                        "value": data_type
                    }
                })
                break
        
        # Detect programs
        if "tcga" in query_lower:
            filters["content"].append({
                "op": "=",
                "content": {
                    "field": "program.name",
                    "value": "TCGA"
                }
            })
        elif "target" in query_lower:
            filters["content"].append({
                "op": "=",
                "content": {
                    "field": "program.name",
                    "value": "TARGET"
                }
            })
        elif "cptac" in query_lower:
            filters["content"].append({
                "op": "=",
                "content": {
                    "field": "program.name",
                    "value": "CPTAC"
                }
            })
        
        # If no specific filters detected, return None to fetch all and filter client-side
        return filters if filters["content"] else None
    
    async def _fetch_with_retry(self, url: str, params: Optional[dict] = None, json_data: Optional[dict] = None) -> httpx.Response:
        """Fetch URL with retry logic. Delegates to the shared base helper.

        Picks ``POST`` automatically when a JSON body is provided, otherwise
        ``GET``. The base helper handles timeout / network / 5xx retries
        with exponential backoff.
        """
        client = await self._get_client()
        method = "POST" if json_data else "GET"
        return await super()._fetch_with_retry(
            client, url, method=method, params=params, json=json_data
        )
    
    def _detect_omics_type(self, project: dict) -> OmicsType:
        """Detect omics type from project data."""
        data_categories = project.get("summary", {}).get("data_categories", [])
        categories = [dc.get("data_category", "") for dc in data_categories if dc]
        
        if len(categories) > 2:
            return OmicsType.MULTI_OMICS
        
        category_str = " ".join(categories).lower()
        
        if "transcriptome" in category_str or "gene expression" in category_str:
            return OmicsType.TRANSCRIPTOMICS
        if "sequencing" in category_str or "genomic" in category_str:
            return OmicsType.GENOMICS
        if "proteomic" in category_str:
            return OmicsType.PROTEOMICS
        
        return OmicsType.MULTI_OMICS
    
    def _detect_assay_type(self, project: dict) -> AssayType:
        """Detect assay type from project data."""
        strategies = project.get("summary", {}).get("experimental_strategies", [])
        strategy_str = " ".join([s.get("experimental_strategy", "") for s in strategies]).lower()
        
        if "rna-seq" in strategy_str or "rna seq" in strategy_str:
            return AssayType.RNA_SEQ
        if "wgs" in strategy_str or "whole genome" in strategy_str:
            return AssayType.WGS
        if "wes" in strategy_str or "whole exome" in strategy_str:
            return AssayType.WES
        if "methylation" in strategy_str:
            return AssayType.WGBS
        if "mirna" in strategy_str:
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
        Search GDC for cancer genomics projects.
        
        Args:
            query: Search query
            max_results: Maximum results to return
            organism: Optional organism filter
            
        Returns:
            List of UnifiedDataset objects
            
        Raises:
            GDCError: On API or network failures
        """
        try:
            # Build request - fetch more if we need to filter client-side
            filters = self._build_filters(query, organism)
            fetch_size = min(max_results, 100) if filters else min(max_results * 2, 200)
            
            params = {
                "from": 0,
                "size": fetch_size,
                "pretty": "false",
            }
            
            url = f"{GDC_API_BASE}/projects"
            
            logger.info("Searching GDC", query=query, max_results=max_results, has_filters=bool(filters))
            
            # Use GET with filters as JSON parameter
            if filters:
                import json
                params["filters"] = json.dumps(filters)
            
            # Use retry-enabled fetch
            response = await self._fetch_with_retry(url, params=params)
            response.raise_for_status()
            
            data = response.json()
            
            # Parse projects from response
            hits = data.get("data", {}).get("hits", [])
            
            logger.info("GDC response", project_count=len(hits))
            
            # Prepare query terms for client-side filtering
            query_lower = query.lower()
            query_terms = [w for w in query_lower.split() if len(w) > 2]
            
            unified = []
            for project in hits:
                if len(unified) >= max_results:
                    break
                    
                try:
                    ud = self._to_unified(project, query=query)
                    if not ud:
                        continue
                    
                    # If no server-side filters, do client-side keyword matching
                    if not filters:
                        title_lower = ud.title.lower()
                        desc_lower = (ud.description or "").lower()
                        combined = f"{title_lower} {desc_lower}"
                        
                        # Check if any query term matches
                        if not any(term in combined for term in query_terms):
                            # Check for cancer-related terms that might be in project_id
                            project_id = project.get("project_id", "").lower()
                            if not any(term in project_id for term in query_terms):
                                continue
                    
                    unified.append(ud)
                except Exception as e:
                    logger.warning("Failed to parse GDC project", error=str(e))
                    continue
            
            # Sort by relevance score
            unified.sort(key=lambda x: x.relevance_score, reverse=True)
            
            logger.info("GDC search complete", results=len(unified))
            return unified
            
        except httpx.HTTPStatusError as e:
            logger.error("GDC API error", status=e.response.status_code, error=str(e))
            raise GDCError(f"GDC API error: {e.response.status_code}") from e
        except httpx.RequestError as e:
            logger.error("GDC request failed", error=str(e))
            raise GDCError(f"GDC request failed: {str(e)}") from e
    
    async def get_dataset(self, accession: str) -> Optional[UnifiedDataset]:
        """Fetch a single GDC project by ID."""
        try:
            url = f"{GDC_API_BASE}/projects/{accession}"
            
            client = await self._get_client()
            response = await client.get(url)
            
            if response.status_code == 404:
                return None
            
            response.raise_for_status()
            data = response.json()
            project = data.get("data", {})
            
            return self._to_unified(project)
            
        except httpx.HTTPStatusError as e:
            if e.response.status_code == 404:
                return None
            logger.error("GDC API error", accession=accession, error=str(e))
            return None
        except Exception as e:
            logger.error("Failed to fetch GDC project", accession=accession, error=str(e))
            return None
    
    def _to_unified(self, project: dict, query: Optional[str] = None) -> Optional[UnifiedDataset]:
        """Convert a GDC project to UnifiedDataset."""
        try:
            project_id = project.get("project_id") or ""
            if not project_id:
                return None
            
            # Extract project info
            name = project.get("name") or project_id
            program = project.get("program", {})
            program_name = program.get("name", "") if isinstance(program, dict) else ""
            
            title = f"{program_name}: {name}" if program_name else name
            
            # Primary site
            primary_site = project.get("primary_site") or []
            if isinstance(primary_site, str):
                primary_site = [primary_site]
            
            # Disease type
            disease_type = project.get("disease_type") or []
            if isinstance(disease_type, str):
                disease_type = [disease_type]
            
            # Build description
            description_parts = []
            if disease_type:
                description_parts.append(f"Disease: {', '.join(disease_type)}")
            if primary_site:
                description_parts.append(f"Primary site: {', '.join(primary_site)}")
            
            description = ". ".join(description_parts) if description_parts else f"{program_name} project"
            
            # Summary stats
            summary = project.get("summary", {})
            case_count = summary.get("case_count", 0)
            file_count = summary.get("file_count", 0)
            
            # Data categories
            data_categories = summary.get("data_categories", [])
            categories = [dc.get("data_category", "") for dc in data_categories if dc]
            
            # Experimental strategies
            exp_strategies = summary.get("experimental_strategies", [])
            strategies = [es.get("experimental_strategy", "") for es in exp_strategies if es]
            
            # Detect types
            omics_type = self._detect_omics_type(project)
            assay_type = self._detect_assay_type(project)
            
            # Compute relevance (uses base class method with GDC stopwords)
            relevance = self._compute_keyword_relevance(query, title, description, GDC_STOPWORDS) if query else 0.5
            
            # Create extension
            extension = GDCExtension(
                project_id=project_id,
                program=program_name,
                primary_site=primary_site[0] if primary_site else None,
                data_categories=categories,
                experimental_strategies=strategies,
                case_count=case_count,
                file_count=file_count,
            )
            
            # Date fallbacks (projects API exposes updated_datetime but not always release date)
            # Note: "released" field is a boolean, not a date
            submission_date = project.get("updated_datetime")
            release_date = submission_date  # Use updated_datetime as release_date
            
            # Ensure release_date is a string, not bool or other type
            if release_date is not None and not isinstance(release_date, str):
                release_date = None
            
            # Build download links for GDC
            download_links = [
                DownloadLink(
                    url=f"https://portal.gdc.cancer.gov/projects/{project_id}",
                    file_type=DownloadFileType.OTHER,
                    description="GDC project page (file repository)",
                    protocol="https",
                ),
                DownloadLink(
                    url=f"https://portal.gdc.cancer.gov/repository?facetTab=files&filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.project_id%22%2C%22value%22%3A%5B%22{project_id}%22%5D%7D%7D%5D%7D",
                    file_type=DownloadFileType.RAW,
                    description="Browse files in GDC repository",
                    protocol="https",
                ),
                DownloadLink(
                    url=f"https://api.gdc.cancer.gov/files?filters=%7B%22op%22%3A%22%3D%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.project_id%22%2C%22value%22%3A%22{project_id}%22%7D%7D&format=json",
                    file_type=DownloadFileType.METADATA,
                    description="GDC API - file list (JSON)",
                    protocol="https",
                ),
            ]

            return UnifiedDataset(
                id=self.build_unified_id(project_id),
                accession=project_id,
                source=self.source,
                source_url=f"https://portal.gdc.cancer.gov/projects/{project_id}",
                title=title,
                description=description,
                organism=["Homo sapiens"],  # GDC is human-only
                sample_count=case_count,
                sample_count_display=str(case_count) if case_count else "N/A",
                submission_date=submission_date,
                release_date=release_date,
                omics_type=omics_type,
                assay_types=[assay_type] if assay_type else [],
                relevance_score=relevance,
                download_links=download_links,
                extensions={"gdc": extension.model_dump()},
            )
            
        except Exception as e:
            logger.warning("Failed to parse GDC project", error=str(e))
            return None
    
    def _create_source_url(self, accession: str) -> str:
        """Create URL to GDC project page."""
        return f"https://portal.gdc.cancer.gov/projects/{accession}"


# Singleton instance
_gdc_adapter: Optional[GDCAdapter] = None


def get_gdc_adapter() -> GDCAdapter:
    """Get or create GDC adapter instance."""
    global _gdc_adapter
    if _gdc_adapter is None:
        _gdc_adapter = GDCAdapter()
    return _gdc_adapter
