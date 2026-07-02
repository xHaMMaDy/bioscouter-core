"""
MGnify (EBI Metagenomics) Adapter
Searches the MGnify database for metagenomics and metatranscriptomics datasets.
API Docs: https://www.ebi.ac.uk/metagenomics/api/v1/
"""

import asyncio
import time
from typing import List, Optional, Dict, Set, Any
from datetime import datetime

import httpx
import structlog

from .base import BaseSourceAdapter, COMMON_STOPWORDS, MAX_RETRY_ATTEMPTS, RETRY_BASE_DELAY, build_default_limits
from ..models.unified import (
    DataSource,
    OmicsType,
    AssayType,
    UnifiedDataset,
    MGnifyExtension,
    DownloadLink,
    DownloadFileType,
)

logger = structlog.get_logger(__name__)

# MGnify API base URL
MGNIFY_API_BASE = "https://www.ebi.ac.uk/metagenomics/api/v1"

# Cache settings
CACHE_TTL = 3600  # 1 hour cache for study lists and biome tree
BIOME_CACHE_TTL = 3600  # 1 hour for biome tree (rarely changes)

# Rate limiting
RATE_LIMIT_REQUESTS_PER_SECOND = 5
RATE_LIMIT_SEMAPHORE_SIZE = 5

# Result limits
MAX_SAMPLE_RESULTS = 500  # Cap for sample-level searches
DEFAULT_PAGE_SIZE = 50

# MGnify-specific stopwords
MGNIFY_STOPWORDS: Set[str] = frozenset({
    'metagenomics', 'metagenomic', 'metagenome', 'microbiome', 'microbiota',
    '16s', 'amplicon', 'shotgun', 'environmental', 'sample', 'samples'
})


class MGnifyError(Exception):
    """MGnify API error."""
    pass


class BiomeNode:
    """Represents a node in the biome hierarchy tree."""
    def __init__(self, id: str, name: str, lineage: str, count: int = 0):
        self.id = id
        self.name = name
        self.lineage = lineage
        self.count = count
        self.children: List['BiomeNode'] = []
    
    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return {
            "id": self.id,
            "name": self.name,
            "lineage": self.lineage,
            "count": self.count,
            "children": [child.to_dict() for child in self.children]
        }


class MGnifyAdapter(BaseSourceAdapter):
    """
    Adapter for MGnify (EBI Metagenomics).
    Uses MGnify REST API v1.
    """
    
    def __init__(self, timeout: float = 30.0):
        super().__init__()
        self.timeout = timeout
        self._client: Optional[httpx.AsyncClient] = None
        # Caches with TTL
        self._studies_cache: Dict[str, Any] = {}
        self._studies_cache_time: float = 0
        self._biome_tree_cache: Optional[List[BiomeNode]] = None
        self._biome_tree_cache_time: float = 0
        # Rate limiting semaphore
        self._rate_semaphore = asyncio.Semaphore(RATE_LIMIT_SEMAPHORE_SIZE)
        self._last_request_time: float = 0
    
    @property
    def source(self) -> DataSource:
        return DataSource.MGNIFY
    
    @property
    def supported_omics(self) -> List[OmicsType]:
        # Supports both metagenomics and transcriptomics (for metatranscriptomics)
        return [OmicsType.METAGENOMICS, OmicsType.TRANSCRIPTOMICS]
    
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
    
    async def _rate_limited_request(self, method: str, url: str, **kwargs) -> httpx.Response:
        """Make a rate-limited HTTP request."""
        async with self._rate_semaphore:
            # Ensure minimum delay between requests
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

        MGnify wraps every request in a rate limiter, so we can't fully
        delegate to :meth:`BaseSourceAdapter._fetch_with_retry` (which
        expects a plain client). The loop below mirrors the shared
        helper's behavior — exponential backoff on timeouts, network
        errors, and 5xx — but goes through ``_rate_limited_request`` for
        each attempt.
        """
        last_exception = None
        for attempt in range(MAX_RETRY_ATTEMPTS):
            try:
                response = await self._rate_limited_request("GET", url)
                if response.status_code >= 500 and attempt < MAX_RETRY_ATTEMPTS - 1:
                    delay = RETRY_BASE_DELAY * (2 ** attempt)
                    logger.warning(
                        "MGnify server error, retrying",
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
                        "MGnify request failed, retrying",
                        error=str(e), attempt=attempt + 1, delay=delay,
                    )
                    await asyncio.sleep(delay)
                    continue
                break
        if last_exception is not None:
            raise last_exception
        raise MGnifyError("Request failed after retries")
    
    async def close(self):
        """Close HTTP client."""
        if self._client and not self._client.is_closed:
            await self._client.aclose()
            self._client = None
    
    async def get_biome_tree(self) -> List[dict]:
        """
        Fetch the full biome hierarchy tree from MGnify.
        Returns a nested tree structure for the biome filter UI.
        Cached for 1 hour.
        """
        now = time.time()
        
        # Return cached tree if still valid
        if self._biome_tree_cache and (now - self._biome_tree_cache_time) < BIOME_CACHE_TTL:
            logger.debug("Using cached biome tree")
            return [node.to_dict() for node in self._biome_tree_cache]
        
        try:
            # Fetch all biomes from MGnify
            url = f"{MGNIFY_API_BASE}/biomes?page_size=1000"
            response = await self._fetch_with_retry(url)
            response.raise_for_status()
            data = response.json()
            
            biomes_data = data.get("data", [])
            
            # Build tree structure from flat list
            root_nodes: List[BiomeNode] = []
            nodes_by_lineage: Dict[str, BiomeNode] = {}
            
            for biome in biomes_data:
                attrs = biome.get("attributes", {})
                biome_id = biome.get("id", "")
                lineage = attrs.get("lineage", "")
                
                # Parse lineage to get name and parent
                parts = lineage.split(":")
                name = parts[-1] if parts else biome_id
                
                # Get study count
                samples_count = attrs.get("samples-count", 0)
                
                node = BiomeNode(
                    id=biome_id,
                    name=name,
                    lineage=lineage,
                    count=samples_count
                )
                nodes_by_lineage[lineage] = node
                
                # Find parent
                if len(parts) > 1:
                    parent_lineage = ":".join(parts[:-1])
                    parent_node = nodes_by_lineage.get(parent_lineage)
                    if parent_node:
                        parent_node.children.append(node)
                    else:
                        # Parent not found yet, add to roots for now
                        root_nodes.append(node)
                else:
                    root_nodes.append(node)
            
            # Re-organize orphans after all nodes are created
            final_roots = []
            for node in root_nodes:
                parts = node.lineage.split(":")
                if len(parts) > 1:
                    parent_lineage = ":".join(parts[:-1])
                    parent_node = nodes_by_lineage.get(parent_lineage)
                    if parent_node and node not in parent_node.children:
                        parent_node.children.append(node)
                        continue
                final_roots.append(node)
            
            self._biome_tree_cache = final_roots
            self._biome_tree_cache_time = now
            
            logger.info("Refreshed MGnify biome tree cache", root_count=len(final_roots))
            return [node.to_dict() for node in final_roots]
            
        except Exception as e:
            logger.error("Failed to fetch biome tree", error=str(e))
            # Return empty list on error, don't crash
            return []
    
    def _simplify_query(self, query: str) -> str:
        """Simplify a natural language query for MGnify search."""
        all_stopwords = COMMON_STOPWORDS | MGNIFY_STOPWORDS
        
        words = query.lower().split()
        filtered = [w for w in words if w not in all_stopwords and len(w) > 2]
        
        result = ' '.join(filtered[:5]) if filtered else query.split()[0]
        logger.debug("Simplified MGnify query", original=query, simplified=result)
        return result
    
    def _compute_keyword_relevance(self, query: str, title: str, description: str) -> float:
        """Compute relevance score based on keyword matching."""
        if not query:
            return 0.5
        
        query_lower = query.lower()
        all_stopwords = COMMON_STOPWORDS | MGNIFY_STOPWORDS
        
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
    
    def _detect_assay_types(self, study: dict) -> List[AssayType]:
        """Detect assay types from study metadata."""
        assay_types = []
        
        experiment_type = (study.get("experiment_type") or "").lower()
        biomes = study.get("biomes", [])
        study_name = (study.get("study_name") or "").lower()
        study_abstract = (study.get("study_abstract") or "").lower()
        combined = f"{experiment_type} {study_name} {study_abstract}"
        
        # Detect metatranscriptomics
        if "metatranscriptom" in combined or "rna" in experiment_type:
            assay_types.append(AssayType.METATRANSCRIPTOMICS)
        
        # Detect amplicon sequencing types
        if "amplicon" in experiment_type or "16s" in combined:
            assay_types.append(AssayType.AMPLICON_16S)
        if "18s" in combined:
            assay_types.append(AssayType.AMPLICON_18S)
        if "its" in combined and ("fung" in combined or "amplicon" in combined):
            assay_types.append(AssayType.AMPLICON_ITS)
        
        # Detect shotgun metagenomics
        if "metagenomic" in experiment_type or "shotgun" in combined or "wgs" in combined:
            assay_types.append(AssayType.SHOTGUN_METAGENOMICS)
        
        # Default to shotgun if nothing detected and not amplicon
        if not assay_types:
            if "amplicon" not in combined:
                assay_types.append(AssayType.SHOTGUN_METAGENOMICS)
            else:
                assay_types.append(AssayType.AMPLICON_16S)
        
        return assay_types
    
    def _is_metatranscriptomics(self, study: dict) -> bool:
        """Check if study is metatranscriptomics."""
        experiment_type = (study.get("experiment_type") or "").lower()
        study_name = (study.get("study_name") or "").lower()
        study_abstract = (study.get("study_abstract") or "").lower()
        combined = f"{experiment_type} {study_name} {study_abstract}"
        
        return "metatranscriptom" in combined or "rna" in experiment_type
    
    async def search(
        self,
        query: str,
        max_results: int = 50,
        organism: Optional[str] = None,
        biome: Optional[str] = None,
        search_level: str = "study",
        **kwargs
    ) -> List[UnifiedDataset]:
        """
        Search MGnify for metagenomics/metatranscriptomics datasets.
        
        Args:
            query: Search query
            max_results: Maximum results to return
            organism: Optional organism filter
            biome: Optional biome lineage filter (e.g., "root:Environmental:Aquatic:Marine")
            search_level: "study" (default) or "sample" for granular search
            
        Returns:
            List of UnifiedDataset objects
        """
        try:
            logger.info("Searching MGnify", query=query, max_results=max_results, 
                       biome=biome, search_level=search_level)
            
            simplified_query = self._simplify_query(query)
            
            if search_level == "sample":
                return await self._search_samples(simplified_query, max_results, organism, biome, query)
            else:
                return await self._search_studies(simplified_query, max_results, organism, biome, query)
                
        except httpx.HTTPStatusError as e:
            logger.error("MGnify API error", status=e.response.status_code, error=str(e))
            raise MGnifyError(f"MGnify API error: {e.response.status_code}") from e
        except httpx.RequestError as e:
            logger.error("MGnify request failed", error=str(e))
            raise MGnifyError(f"MGnify request failed: {str(e)}") from e
    
    async def _search_studies(
        self,
        simplified_query: str,
        max_results: int,
        organism: Optional[str],
        biome: Optional[str],
        original_query: str
    ) -> List[UnifiedDataset]:
        """Search at study level."""
        unified = []
        page = 1
        page_size = min(max_results, DEFAULT_PAGE_SIZE)
        
        while len(unified) < max_results:
            # Build search URL
            params = {
                "search": simplified_query,
                "page": page,
                "page_size": page_size,
                "ordering": "-samples_count"  # Most samples first
            }
            
            if biome:
                params["lineage"] = biome
            
            url = f"{MGNIFY_API_BASE}/studies?{'&'.join(f'{k}={v}' for k, v in params.items())}"
            
            response = await self._fetch_with_retry(url)
            
            if response.status_code == 404:
                break
            
            response.raise_for_status()
            data = response.json()
            
            studies = data.get("data", [])
            
            if not studies:
                break
            
            for study in studies:
                if len(unified) >= max_results:
                    break
                
                attrs = study.get("attributes", {})
                
                # Filter by organism if specified
                if organism:
                    # MGnify doesn't have direct organism filter, check in abstract
                    abstract = (attrs.get("study-abstract") or "").lower()
                    if organism.lower() not in abstract:
                        continue
                
                ud = self._study_to_unified(study, original_query)
                if ud:
                    unified.append(ud)
            
            # Check if there are more pages
            links = data.get("links", {})
            if not links.get("next"):
                break
            
            page += 1
        
        # Sort by relevance
        unified.sort(key=lambda x: x.relevance_score, reverse=True)
        
        logger.info("MGnify study search complete", results=len(unified))
        return unified[:max_results]
    
    async def _search_samples(
        self,
        simplified_query: str,
        max_results: int,
        organism: Optional[str],
        biome: Optional[str],
        original_query: str
    ) -> List[UnifiedDataset]:
        """Search at sample level (more granular, capped at MAX_SAMPLE_RESULTS)."""
        unified = []
        page = 1
        page_size = min(max_results, DEFAULT_PAGE_SIZE)
        effective_max = min(max_results, MAX_SAMPLE_RESULTS)
        
        while len(unified) < effective_max:
            params = {
                "search": simplified_query,
                "page": page,
                "page_size": page_size,
            }
            
            if biome:
                params["lineage"] = biome
            
            url = f"{MGNIFY_API_BASE}/samples?{'&'.join(f'{k}={v}' for k, v in params.items())}"
            
            response = await self._fetch_with_retry(url)
            
            if response.status_code == 404:
                break
            
            response.raise_for_status()
            data = response.json()
            
            samples = data.get("data", [])
            
            if not samples:
                break
            
            for sample in samples:
                if len(unified) >= effective_max:
                    break
                
                attrs = sample.get("attributes", {})
                
                # Filter by organism if specified
                if organism:
                    sample_desc = (attrs.get("sample-desc") or "").lower()
                    if organism.lower() not in sample_desc:
                        continue
                
                ud = self._sample_to_unified(sample, original_query)
                if ud:
                    unified.append(ud)
            
            links = data.get("links", {})
            if not links.get("next"):
                break
            
            page += 1
        
        unified.sort(key=lambda x: x.relevance_score, reverse=True)
        
        logger.info("MGnify sample search complete", results=len(unified))
        return unified[:effective_max]
    
    def _study_to_unified(self, study: dict, query: str) -> Optional[UnifiedDataset]:
        """Convert MGnify study to UnifiedDataset."""
        try:
            study_id = study.get("id", "")
            attrs = study.get("attributes", {})
            
            study_name = attrs.get("study-name") or attrs.get("bioproject") or study_id
            study_abstract = attrs.get("study-abstract") or ""
            
            # Get sample count
            sample_count = attrs.get("samples-count", 0)
            
            # Get biomes
            biomes = []
            relationships = study.get("relationships", {})
            biomes_rel = relationships.get("biomes", {})
            biomes_data = biomes_rel.get("data", [])
            for b in biomes_data:
                if isinstance(b, dict):
                    biomes.append(b.get("id", ""))
            
            # Detect experiment type and assay types
            is_metatranscriptomic = self._is_metatranscriptomics({"study_name": study_name, "study_abstract": study_abstract})
            assay_types = self._detect_assay_types({"study_name": study_name, "study_abstract": study_abstract})
            
            # Determine omics type
            omics_type = OmicsType.TRANSCRIPTOMICS if is_metatranscriptomic else OmicsType.METAGENOMICS
            
            # Get secondary accession (ENA/SRA) for deduplication
            secondary_accession = attrs.get("bioproject") or attrs.get("secondary-accession")
            
            # Compute relevance
            relevance = self._compute_keyword_relevance(query, study_name, study_abstract)
            
            # Generate download links
            download_links = self._generate_download_links(study_id, attrs)
            
            # Build extension
            extension = MGnifyExtension(
                study_id=study_id,
                secondary_accession=secondary_accession,
                biomes=biomes,
                pipeline_version=attrs.get("pipeline-version"),
                sample_count=sample_count,
                analysis_count=attrs.get("analyses-count", 0),
                experiment_type=attrs.get("experiment-type"),
                top_phyla=[],  # Would need additional API call to fetch
                functional_categories=[],
                total_size_bytes=None,  # Would need additional API call
                is_metatranscriptomics=is_metatranscriptomic,
            )
            
            return UnifiedDataset(
                id=self.build_unified_id(study_id),
                accession=study_id,
                source=self.source,
                source_url=f"https://www.ebi.ac.uk/metagenomics/studies/{study_id}",
                omics_type=omics_type,
                assay_types=assay_types,
                title=study_name,
                description=study_abstract,
                organism=[],  # MGnify studies cover many organisms
                sample_count=sample_count,
                sample_count_display=f"{sample_count:,}" if sample_count else "0",
                submission_date=attrs.get("last-update"),
                release_date=attrs.get("last-update"),
                relevance_score=relevance,
                match_reasons=[f"MGnify study matching '{query}'"],
                download_links=download_links,
                extensions={"mgnify": extension.model_dump()},
            )
            
        except Exception as e:
            logger.error("Failed to convert MGnify study", study_id=study.get("id"), error=str(e))
            return None
    
    def _sample_to_unified(self, sample: dict, query: str) -> Optional[UnifiedDataset]:
        """Convert MGnify sample to UnifiedDataset."""
        try:
            sample_id = sample.get("id", "")
            attrs = sample.get("attributes", {})
            
            sample_name = attrs.get("sample-name") or attrs.get("sample-alias") or sample_id
            sample_desc = attrs.get("sample-desc") or ""
            
            # Get biome
            biomes = []
            biome_lineage = attrs.get("environment-biome") or ""
            if biome_lineage:
                biomes.append(biome_lineage)
            
            # Compute relevance
            relevance = self._compute_keyword_relevance(query, sample_name, sample_desc)
            
            # Build extension
            extension = MGnifyExtension(
                study_id=sample_id,
                secondary_accession=attrs.get("accession"),
                biomes=biomes,
                sample_count=1,
                analysis_count=attrs.get("analyses-count", 0),
                is_metatranscriptomics=False,
            )
            
            return UnifiedDataset(
                id=self.build_unified_id(sample_id),
                accession=sample_id,
                source=self.source,
                source_url=f"https://www.ebi.ac.uk/metagenomics/samples/{sample_id}",
                omics_type=OmicsType.METAGENOMICS,
                assay_types=[AssayType.SHOTGUN_METAGENOMICS],
                title=sample_name,
                description=sample_desc,
                organism=[],
                sample_count=1,
                sample_count_display="1",
                relevance_score=relevance,
                match_reasons=[f"MGnify sample matching '{query}'"],
                download_links=[],
                extensions={"mgnify": extension.model_dump()},
            )
            
        except Exception as e:
            logger.error("Failed to convert MGnify sample", sample_id=sample.get("id"), error=str(e))
            return None
    
    def _generate_download_links(self, study_id: str, attrs: dict) -> List[DownloadLink]:
        """Generate download links for an MGnify study."""
        links = []
        
        # Study page
        links.append(DownloadLink(
            url=f"https://www.ebi.ac.uk/metagenomics/studies/{study_id}",
            file_type=DownloadFileType.OTHER,
            protocol="https",
            description="MGnify study page",
        ))
        
        # Analysis results download
        links.append(DownloadLink(
            url=f"https://www.ebi.ac.uk/metagenomics/api/v1/studies/{study_id}/downloads",
            file_type=DownloadFileType.PROCESSED,
            protocol="https",
            description="Analysis results (TSV, JSON)",
            needs_refresh=True,
        ))
        
        # ENA link if available
        bioproject = attrs.get("bioproject")
        if bioproject:
            links.append(DownloadLink(
                url=f"https://www.ebi.ac.uk/ena/browser/view/{bioproject}",
                file_type=DownloadFileType.RAW,
                protocol="https",
                description=f"Raw reads on ENA ({bioproject})",
            ))
        
        return links
    
    async def get_dataset(self, accession: str) -> Optional[UnifiedDataset]:
        """Fetch a single study by accession."""
        try:
            url = f"{MGNIFY_API_BASE}/studies/{accession}"
            response = await self._fetch_with_retry(url)
            
            if response.status_code == 404:
                return None
            
            response.raise_for_status()
            data = response.json()
            study = data.get("data", {})
            
            return self._study_to_unified(study, accession)
            
        except Exception as e:
            logger.error("Failed to fetch MGnify study", accession=accession, error=str(e))
            return None


# Singleton instance
_mgnify_adapter: Optional[MGnifyAdapter] = None


def get_mgnify_adapter() -> MGnifyAdapter:
    """Get or create MGnify adapter instance."""
    global _mgnify_adapter
    if _mgnify_adapter is None:
        _mgnify_adapter = MGnifyAdapter()
    return _mgnify_adapter
