"""
MG-RAST Adapter
Searches the MG-RAST database for metagenomics and metatranscriptomics datasets.
API Docs: https://api.mg-rast.org/api.html
"""

import asyncio
import time
from typing import List, Optional, Dict, Set, Any

import httpx
import structlog

from .base import BaseSourceAdapter, COMMON_STOPWORDS, MAX_RETRY_ATTEMPTS, RETRY_BASE_DELAY, build_default_limits
from ..models.unified import (
    DataSource,
    OmicsType,
    AssayType,
    UnifiedDataset,
    MGRASTExtension,
    DownloadLink,
    DownloadFileType,
)
from ..core.config import get_settings

logger = structlog.get_logger(__name__)

# MG-RAST API base URL
MGRAST_API_BASE = "https://api.mg-rast.org"

# Cache settings
CACHE_TTL = 3600  # 1 hour cache for project listings

# Rate limiting - MG-RAST recommends conservative rate limits
RATE_LIMIT_REQUESTS_PER_SECOND = 2
RATE_LIMIT_SEMAPHORE_SIZE = 2

# Result limits
MAX_METAGENOME_RESULTS = 500  # Cap for metagenome-level searches
DEFAULT_PAGE_SIZE = 50

# MG-RAST-specific stopwords
MGRAST_STOPWORDS: Set[str] = frozenset({
    'metagenomics', 'metagenomic', 'metagenome', 'microbiome', 'microbiota',
    '16s', 'amplicon', 'shotgun', 'environmental', 'sample', 'samples',
    'mgrast', 'mg-rast', 'project', 'analysis'
})


class MGRASTError(Exception):
    """MG-RAST API error."""
    pass


class MGRASTAdapter(BaseSourceAdapter):
    """
    Adapter for MG-RAST (Metagenomics RAST).
    Uses MG-RAST REST API.
    """
    
    def __init__(self, timeout: float = 30.0):
        super().__init__()
        self.timeout = timeout
        self._client: Optional[httpx.AsyncClient] = None
        # Cache with TTL
        self._projects_cache: Dict[str, Any] = {}
        self._projects_cache_time: float = 0
        # Rate limiting semaphore
        self._rate_semaphore = asyncio.Semaphore(RATE_LIMIT_SEMAPHORE_SIZE)
        self._last_request_time: float = 0
    
    @property
    def source(self) -> DataSource:
        return DataSource.MGRAST
    
    @property
    def supported_omics(self) -> List[OmicsType]:
        # Supports both metagenomics and transcriptomics (for metatranscriptomics)
        return [OmicsType.METAGENOMICS, OmicsType.TRANSCRIPTOMICS]
    
    def _get_headers(self) -> dict:
        """Get HTTP headers including optional API key."""
        headers = {
            "Accept": "application/json",
            "User-Agent": "BioScouter/1.0",
        }
        
        # Add API key if configured
        settings = get_settings()
        if settings.mgrast_api_key:
            headers["Authorization"] = f"mgrast {settings.mgrast_api_key}"
        
        return headers
    
    async def _get_client(self) -> httpx.AsyncClient:
        """Get or create HTTP client."""
        if self._client is None or self._client.is_closed:
            self._client = httpx.AsyncClient(
                timeout=self.timeout,
                limits=build_default_limits(),
                headers=self._get_headers()
            )
        return self._client
    
    async def _rate_limited_request(self, method: str, url: str, **kwargs) -> httpx.Response:
        """Make a rate-limited HTTP request."""
        async with self._rate_semaphore:
            # Ensure minimum delay between requests (more conservative for MG-RAST)
            now = time.time()
            min_delay = 1.0 / RATE_LIMIT_REQUESTS_PER_SECOND
            elapsed = now - self._last_request_time
            if elapsed < min_delay:
                await asyncio.sleep(min_delay - elapsed)
            
            client = await self._get_client()
            self._last_request_time = time.time()
            
            if method.upper() == "GET":
                return await client.get(url, **kwargs)
            elif method.upper() == "POST":
                return await client.post(url, **kwargs)
            else:
                raise ValueError(f"Unsupported HTTP method: {method}")
    
    async def _fetch_with_retry(self, url: str) -> httpx.Response:
        """Fetch URL with retry logic.

        Like MGnify, MG-RAST funnels requests through its own rate limiter,
        so we keep a local retry loop that mirrors the shared base helper
        rather than delegating directly.
        """
        last_exception = None
        for attempt in range(MAX_RETRY_ATTEMPTS):
            try:
                response = await self._rate_limited_request("GET", url)
                if response.status_code >= 500 and attempt < MAX_RETRY_ATTEMPTS - 1:
                    delay = RETRY_BASE_DELAY * (2 ** attempt)
                    logger.warning(
                        "MG-RAST server error, retrying",
                        status=response.status_code, attempt=attempt + 1, delay=delay,
                    )
                    await asyncio.sleep(delay)
                    continue
                return response
            except (httpx.TimeoutException, httpx.NetworkError) as e:
                last_exception = e
                if attempt < MAX_RETRY_ATTEMPTS - 1:
                    delay = RETRY_BASE_DELAY * (2 ** attempt)
                    logger.warning(
                        "MG-RAST request failed, retrying",
                        error=str(e), attempt=attempt + 1, delay=delay,
                    )
                    await asyncio.sleep(delay)
                    continue
                break
        if last_exception is not None:
            raise last_exception
        raise MGRASTError("Request failed after retries")
    
    async def close(self):
        """Close HTTP client."""
        if self._client and not self._client.is_closed:
            await self._client.aclose()
            self._client = None
    
    def _simplify_query(self, query: str) -> str:
        """Simplify a natural language query for MG-RAST search."""
        all_stopwords = COMMON_STOPWORDS | MGRAST_STOPWORDS
        
        words = query.lower().split()
        filtered = [w for w in words if w not in all_stopwords and len(w) > 2]
        
        result = ' '.join(filtered[:5]) if filtered else query.split()[0]
        logger.debug("Simplified MG-RAST query", original=query, simplified=result)
        return result
    
    def _compute_keyword_relevance(self, query: str, title: str, description: str) -> float:
        """Compute relevance score based on keyword matching."""
        if not query:
            return 0.5
        
        query_lower = query.lower()
        all_stopwords = COMMON_STOPWORDS | MGRAST_STOPWORDS
        
        query_words = [w.strip() for w in query_lower.split() if len(w.strip()) > 2]
        keywords = [w for w in query_words if w not in all_stopwords]
        
        if not keywords:
            keywords = query_words[:3]
        
        if not keywords:
            return 0.3
        
        title_lower = (title or "").lower()
        desc_lower = (description or "").lower()
        
        title_matches = sum(1 for kw in keywords if kw in title_lower)
        desc_matches = sum(1 for kw in keywords if kw in desc_lower)
        
        total_keywords = len(keywords)
        title_score = title_matches / total_keywords if total_keywords > 0 else 0
        desc_score = desc_matches / total_keywords if total_keywords > 0 else 0
        
        base_score = (title_score * 0.6) + (desc_score * 0.4)
        
        if query_lower in title_lower:
            base_score = min(1.0, base_score + 0.3)
        elif query_lower in desc_lower:
            base_score = min(1.0, base_score + 0.15)
        
        return max(0.3, min(1.0, base_score + 0.2))
    
    def _detect_assay_types(self, metagenome: dict) -> List[AssayType]:
        """Detect assay types from metagenome metadata."""
        assay_types = []
        
        seq_type = (metagenome.get("sequence_type") or "").lower()
        project_name = (metagenome.get("project_name") or "").lower()
        project_desc = (metagenome.get("project_description") or "").lower()
        combined = f"{seq_type} {project_name} {project_desc}"
        
        # Detect metatranscriptomics
        if "metatranscriptom" in combined or seq_type == "mt":
            assay_types.append(AssayType.METATRANSCRIPTOMICS)
        
        # Detect amplicon sequencing
        if "amplicon" in seq_type or "16s" in combined:
            assay_types.append(AssayType.AMPLICON_16S)
        if "18s" in combined:
            assay_types.append(AssayType.AMPLICON_18S)
        if "its" in combined:
            assay_types.append(AssayType.AMPLICON_ITS)
        
        # Detect shotgun metagenomics
        if "wgs" in seq_type or "shotgun" in combined or seq_type == "wgs":
            assay_types.append(AssayType.SHOTGUN_METAGENOMICS)
        
        # Default to shotgun if nothing detected
        if not assay_types:
            assay_types.append(AssayType.SHOTGUN_METAGENOMICS)
        
        return assay_types
    
    def _is_metatranscriptomics(self, metagenome: dict) -> bool:
        """Check if metagenome is metatranscriptomics."""
        seq_type = (metagenome.get("sequence_type") or "").lower()
        project_name = (metagenome.get("project_name") or "").lower()
        combined = f"{seq_type} {project_name}"
        
        return "metatranscriptom" in combined or seq_type == "mt"
    
    async def search(
        self,
        query: str,
        max_results: int = 50,
        organism: Optional[str] = None,
        biome: Optional[str] = None,
        env_package: Optional[str] = None,
        search_level: str = "project",
        **kwargs
    ) -> List[UnifiedDataset]:
        """
        Search MG-RAST for metagenomics/metatranscriptomics datasets.
        
        Args:
            query: Search query
            max_results: Maximum results to return
            organism: Optional organism filter
            biome: Optional biome filter
            env_package: Optional environment package filter
            search_level: "project" (default) or "metagenome" for granular search
            
        Returns:
            List of UnifiedDataset objects
        """
        try:
            logger.info("Searching MG-RAST", query=query, max_results=max_results,
                       biome=biome, search_level=search_level)
            
            simplified_query = self._simplify_query(query)
            
            if search_level == "metagenome":
                return await self._search_metagenomes(simplified_query, max_results, organism, biome, env_package, query)
            else:
                return await self._search_projects(simplified_query, max_results, organism, biome, env_package, query)
                
        except httpx.HTTPStatusError as e:
            logger.error("MG-RAST API error", status=e.response.status_code, error=str(e))
            # Return empty list instead of raising to allow partial results from MGnify
            return []
        except httpx.RequestError as e:
            logger.error("MG-RAST request failed", error=str(e))
            # Return empty list instead of raising to allow partial results from MGnify
            return []
        except Exception as e:
            logger.error("MG-RAST unexpected error", error=str(e))
            return []
    
    async def _search_projects(
        self,
        simplified_query: str,
        max_results: int,
        organism: Optional[str],
        biome: Optional[str],
        env_package: Optional[str],
        original_query: str
    ) -> List[UnifiedDataset]:
        """Search at project level."""
        unified = []
        offset = 0
        limit = min(max_results, DEFAULT_PAGE_SIZE)
        
        while len(unified) < max_results:
            # Build search URL using MG-RAST project search
            # Note: verbosity=full is not allowed for list endpoints
            url = f"{MGRAST_API_BASE}/project?limit={limit}&offset={offset}"
            
            # Add search query
            if simplified_query:
                url += f"&match=all&name={simplified_query}"
            
            # Add biome/env filters if specified
            if biome:
                url += f"&biome={biome}"
            if env_package:
                url += f"&env_package_type={env_package}"
            
            response = await self._fetch_with_retry(url)
            
            if response.status_code == 404:
                break
            
            if response.status_code != 200:
                logger.warning("MG-RAST project search returned non-200", status=response.status_code)
                break
            
            data = response.json()
            projects = data.get("data", [])
            
            if not projects:
                break
            
            for project in projects:
                if len(unified) >= max_results:
                    break
                
                # Filter by organism if specified
                if organism:
                    # Check in project metadata
                    project_name = (project.get("name") or "").lower()
                    project_desc = (project.get("description") or "").lower()
                    if organism.lower() not in f"{project_name} {project_desc}":
                        continue
                
                ud = self._project_to_unified(project, original_query)
                if ud:
                    unified.append(ud)
            
            # Check if there are more results
            total = data.get("total_count", 0)
            offset += limit
            if offset >= total:
                break
        
        # Sort by relevance
        unified.sort(key=lambda x: x.relevance_score, reverse=True)
        
        logger.info("MG-RAST project search complete", results=len(unified))
        return unified[:max_results]
    
    async def _search_metagenomes(
        self,
        simplified_query: str,
        max_results: int,
        organism: Optional[str],
        biome: Optional[str],
        env_package: Optional[str],
        original_query: str
    ) -> List[UnifiedDataset]:
        """Search at metagenome level (more granular, capped at MAX_METAGENOME_RESULTS)."""
        unified = []
        offset = 0
        limit = min(max_results, DEFAULT_PAGE_SIZE)
        effective_max = min(max_results, MAX_METAGENOME_RESULTS)
        
        while len(unified) < effective_max:
            # Build search URL using MG-RAST metagenome search
            # Note: verbosity=full is not allowed for list endpoints
            url = f"{MGRAST_API_BASE}/metagenome?limit={limit}&offset={offset}"
            
            # Add search query
            if simplified_query:
                url += f"&match=all&name={simplified_query}"
            
            # Add biome/env filters if specified
            if biome:
                url += f"&biome={biome}"
            if env_package:
                url += f"&env_package_type={env_package}"
            
            response = await self._fetch_with_retry(url)
            
            if response.status_code == 404:
                break
            
            if response.status_code != 200:
                logger.warning("MG-RAST metagenome search returned non-200", status=response.status_code)
                break
            
            data = response.json()
            metagenomes = data.get("data", [])
            
            if not metagenomes:
                break
            
            for mg in metagenomes:
                if len(unified) >= effective_max:
                    break
                
                # Filter by organism if specified
                if organism:
                    mg_name = (mg.get("name") or "").lower()
                    if organism.lower() not in mg_name:
                        continue
                
                ud = self._metagenome_to_unified(mg, original_query)
                if ud:
                    unified.append(ud)
            
            total = data.get("total_count", 0)
            offset += limit
            if offset >= total:
                break
        
        unified.sort(key=lambda x: x.relevance_score, reverse=True)
        
        logger.info("MG-RAST metagenome search complete", results=len(unified))
        return unified[:effective_max]
    
    def _project_to_unified(self, project: dict, query: str) -> Optional[UnifiedDataset]:
        """Convert MG-RAST project to UnifiedDataset."""
        try:
            project_id = project.get("id", "")
            
            project_name = project.get("name") or project_id
            project_desc = project.get("description") or ""
            
            # Get metagenome count
            metagenome_count = len(project.get("metagenomes", [])) or project.get("metagenome_count", 0)
            
            # Determine if metatranscriptomics
            is_metatranscriptomic = self._is_metatranscriptomics({
                "project_name": project_name,
                "project_description": project_desc
            })
            assay_types = self._detect_assay_types({
                "project_name": project_name,
                "project_description": project_desc
            })
            
            omics_type = OmicsType.TRANSCRIPTOMICS if is_metatranscriptomic else OmicsType.METAGENOMICS
            
            # Get metadata
            metadata = project.get("metadata", {})
            env_biome = metadata.get("env_biome", {}).get("value") if isinstance(metadata.get("env_biome"), dict) else metadata.get("env_biome")
            env_feature = metadata.get("env_feature", {}).get("value") if isinstance(metadata.get("env_feature"), dict) else metadata.get("env_feature")
            env_material = metadata.get("env_material", {}).get("value") if isinstance(metadata.get("env_material"), dict) else metadata.get("env_material")
            env_package = metadata.get("env_package", {}).get("value") if isinstance(metadata.get("env_package"), dict) else metadata.get("env_package")
            
            # Try to find secondary accession (SRA/ENA)
            secondary_accession = metadata.get("external_id") or metadata.get("ncbi_id")
            
            # Compute relevance
            relevance = self._compute_keyword_relevance(query, project_name, project_desc)
            
            # Generate download links
            download_links = self._generate_download_links(project_id, project)
            
            # Build extension
            extension = MGRASTExtension(
                project_id=project_id,
                secondary_accession=secondary_accession,
                biome=env_biome,
                feature=env_feature,
                material=env_material,
                env_package=env_package,
                sequence_type=project.get("sequence_type"),
                bp_count=project.get("bp_count"),
                sequence_count=project.get("sequence_count"),
                top_phyla=[],
                functional_categories=[],
                total_size_bytes=None,
                is_metatranscriptomics=is_metatranscriptomic,
            )
            
            return UnifiedDataset(
                id=self.build_unified_id(project_id),
                accession=project_id,
                source=self.source,
                source_url=f"https://www.mg-rast.org/mgmain.html?mgpage=project&project={project_id}",
                omics_type=omics_type,
                assay_types=assay_types,
                title=project_name,
                description=project_desc,
                organism=[],  # MG-RAST projects cover many organisms
                sample_count=metagenome_count,
                sample_count_display=f"{metagenome_count:,}" if metagenome_count else "0",
                submission_date=project.get("created_on"),
                release_date=project.get("created_on"),
                relevance_score=relevance,
                match_reasons=[f"MG-RAST project matching '{query}'"],
                download_links=download_links,
                extensions={"mgrast": extension.model_dump()},
            )
            
        except Exception as e:
            logger.error("Failed to convert MG-RAST project", project_id=project.get("id"), error=str(e))
            return None
    
    def _metagenome_to_unified(self, mg: dict, query: str) -> Optional[UnifiedDataset]:
        """Convert MG-RAST metagenome to UnifiedDataset."""
        try:
            mg_id = mg.get("id", "")
            
            mg_name = mg.get("name") or mg_id
            project_name = mg.get("project_name") or ""
            
            # Determine if metatranscriptomics
            is_metatranscriptomic = self._is_metatranscriptomics(mg)
            assay_types = self._detect_assay_types(mg)
            
            omics_type = OmicsType.TRANSCRIPTOMICS if is_metatranscriptomic else OmicsType.METAGENOMICS
            
            # Get metadata
            metadata = mg.get("metadata", {})
            env_biome = metadata.get("env_biome")
            env_feature = metadata.get("env_feature")
            env_material = metadata.get("env_material")
            env_package = mg.get("env_package_type")
            
            # Try to find secondary accession
            secondary_accession = metadata.get("external_id") or metadata.get("ncbi_id")
            
            # Compute relevance
            relevance = self._compute_keyword_relevance(query, mg_name, project_name)
            
            # Build extension
            extension = MGRASTExtension(
                mgm_id=mg_id,
                project_id=mg.get("project_id") or "",
                secondary_accession=secondary_accession,
                biome=env_biome,
                feature=env_feature,
                material=env_material,
                env_package=env_package,
                sequence_type=mg.get("sequence_type"),
                bp_count=mg.get("bp_count"),
                sequence_count=mg.get("sequence_count"),
                is_metatranscriptomics=is_metatranscriptomic,
            )
            
            # Generate download links
            download_links = [
                DownloadLink(
                    url=f"https://www.mg-rast.org/mgmain.html?mgpage=overview&metagenome={mg_id}",
                    file_type=DownloadFileType.OTHER,
                    protocol="https",
                    description="MG-RAST metagenome page",
                ),
                DownloadLink(
                    url=f"{MGRAST_API_BASE}/download/{mg_id}?file=050.1",
                    file_type=DownloadFileType.RAW,
                    protocol="https",
                    description="Uploaded sequences (FASTA)",
                ),
            ]
            
            return UnifiedDataset(
                id=self.build_unified_id(mg_id),
                accession=mg_id,
                source=self.source,
                source_url=f"https://www.mg-rast.org/mgmain.html?mgpage=overview&metagenome={mg_id}",
                omics_type=omics_type,
                assay_types=assay_types,
                title=mg_name,
                description=f"Part of project: {project_name}" if project_name else "",
                organism=[],
                sample_count=1,
                sample_count_display="1",
                relevance_score=relevance,
                match_reasons=[f"MG-RAST metagenome matching '{query}'"],
                download_links=download_links,
                extensions={"mgrast": extension.model_dump()},
            )
            
        except Exception as e:
            logger.error("Failed to convert MG-RAST metagenome", mg_id=mg.get("id"), error=str(e))
            return None
    
    def _generate_download_links(self, project_id: str, project: dict) -> List[DownloadLink]:
        """Generate download links for an MG-RAST project."""
        links = []
        
        # Project page
        links.append(DownloadLink(
            url=f"https://www.mg-rast.org/mgmain.html?mgpage=project&project={project_id}",
            file_type=DownloadFileType.OTHER,
            protocol="https",
            description="MG-RAST project page",
        ))
        
        # API download endpoint for project data
        links.append(DownloadLink(
            url=f"{MGRAST_API_BASE}/project/{project_id}?verbosity=full",
            file_type=DownloadFileType.METADATA,
            protocol="https",
            description="Project metadata (JSON)",
        ))
        
        # Add links for individual metagenomes if available
        metagenomes = project.get("metagenomes", [])[:5]  # Limit to first 5
        for mg in metagenomes:
            mg_id = mg[0] if isinstance(mg, list) else mg.get("metagenome_id", mg)
            if mg_id:
                links.append(DownloadLink(
                    url=f"{MGRAST_API_BASE}/download/{mg_id}?file=050.1",
                    file_type=DownloadFileType.RAW,
                    protocol="https",
                    description=f"Sequences for {mg_id} (FASTA)",
                    needs_refresh=True,
                ))
        
        return links
    
    async def get_dataset(self, accession: str) -> Optional[UnifiedDataset]:
        """Fetch a single project or metagenome by accession."""
        try:
            # Determine if it's a project (mgp) or metagenome (mgm)
            if accession.startswith("mgp") or accession.startswith("mgp"):
                url = f"{MGRAST_API_BASE}/project/{accession}?verbosity=full"
                response = await self._fetch_with_retry(url)
                
                if response.status_code == 404:
                    return None
                
                if response.status_code != 200:
                    return None
                
                data = response.json()
                return self._project_to_unified(data, accession)
            else:
                # Assume metagenome
                url = f"{MGRAST_API_BASE}/metagenome/{accession}?verbosity=full"
                response = await self._fetch_with_retry(url)
                
                if response.status_code == 404:
                    return None
                
                if response.status_code != 200:
                    return None
                
                data = response.json()
                return self._metagenome_to_unified(data, accession)
            
        except Exception as e:
            logger.error("Failed to fetch MG-RAST dataset", accession=accession, error=str(e))
            return None


# Singleton instance
_mgrast_adapter: Optional[MGRASTAdapter] = None


def get_mgrast_adapter() -> MGRASTAdapter:
    """Get or create MG-RAST adapter instance."""
    global _mgrast_adapter
    if _mgrast_adapter is None:
        _mgrast_adapter = MGRASTAdapter()
    return _mgrast_adapter
