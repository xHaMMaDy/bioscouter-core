"""
MetaboLights Adapter
Searches the MetaboLights database for metabolomics datasets.
API Docs: https://www.ebi.ac.uk/metabolights/ws/
"""

import asyncio
import re
import time
from typing import List, Optional, Any, Dict, Set
from datetime import datetime
from urllib.parse import urlencode

import httpx
import structlog

from .base import BaseSourceAdapter, COMMON_STOPWORDS, build_default_limits
from ..models.unified import (
    DataSource,
    OmicsType,
    AssayType,
    UnifiedDataset,
    MetaboLightsExtension,
    DownloadLink,
    DownloadFileType,
)

logger = structlog.get_logger(__name__)

# MetaboLights API base URL
METABOLIGHTS_API_BASE = "https://www.ebi.ac.uk/metabolights/ws"

# Cache settings
STUDY_LIST_CACHE_TTL = 3600  # 1 hour - refresh study list every hour to catch new studies
MAX_STUDY_LIST_SIZE = 5000  # Limit cached study IDs to prevent memory bloat

# Batch processing settings
BATCH_SIZE = 15  # Number of concurrent study fetches
SEARCH_MULTIPLIER = 4  # Check 4x max_results to account for filtering

# MetaboLights-specific stopwords
METABOLIGHTS_STOPWORDS: Set[str] = frozenset({
    'metabolomics', 'metabolomic', 'metabolome', 'metabolite', 'metabolites',
    'mass', 'spectrometry', 'ms', 'nmr', 'profiling'
})


class MetaboLightsError(Exception):
    """MetaboLights API error."""
    pass


class MetaboLightsAdapter(BaseSourceAdapter):
    """
    Adapter for MetaboLights (metabolomics database).
    Uses MetaboLights REST API.
    """
    
    def __init__(self, timeout: float = 30.0):
        super().__init__()
        self.timeout = timeout
        self._client: Optional[httpx.AsyncClient] = None
        # Cache for study list with TTL
        self._study_list_cache: List[str] = []
        self._study_list_cache_time: float = 0
    
    @property
    def source(self) -> DataSource:
        return DataSource.METABOLIGHTS
    
    @property
    def supported_omics(self) -> List[OmicsType]:
        return [OmicsType.METABOLOMICS]
    
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
    
    async def _get_study_list(self) -> List[str]:
        """Get list of all study IDs with caching."""
        now = time.time()
        
        # Return cached list if still valid
        if self._study_list_cache and (now - self._study_list_cache_time) < STUDY_LIST_CACHE_TTL:
            logger.debug("Using cached study list", count=len(self._study_list_cache))
            return self._study_list_cache
        
        # Fetch fresh list
        url = f"{METABOLIGHTS_API_BASE}/studies"
        client = await self._get_client()
        response = await client.get(url)
        response.raise_for_status()
        
        data = response.json()
        
        # MetaboLights returns a list of study IDs
        study_ids = data.get("content", []) if isinstance(data, dict) else data
        if isinstance(study_ids, dict):
            study_ids = list(study_ids.keys())
        
        # Update cache (limit size to prevent memory bloat)
        self._study_list_cache = study_ids[:MAX_STUDY_LIST_SIZE]
        self._study_list_cache_time = now
        
        logger.info("Refreshed MetaboLights study list cache", count=len(self._study_list_cache))
        return study_ids
    
    async def close(self):
        """Close HTTP client."""
        if self._client and not self._client.is_closed:
            await self._client.aclose()
            self._client = None
    
    def _simplify_query(self, query: str) -> str:
        """
        Simplify a natural language query for MetaboLights search.
        """
        # Remove common metabolomics terms that don't help search
        stopwords = {
            'metabolomics', 'metabolomic', 'metabolome', 'metabolite', 'metabolites',
            'mass', 'spectrometry', 'ms', 'nmr', 'data', 'dataset', 'datasets',
            'analysis', 'study', 'experiment', 'experiments', 'research',
            'large', 'high', 'quality', 'based', 'using', 'profiling'
        }
        
        words = query.lower().split()
        filtered = [w for w in words if w not in stopwords and len(w) > 2]
        
        result = ' '.join(filtered[:4]) if filtered else query.split()[0]
        logger.info("Simplified MetaboLights query", original=query, simplified=result)
        return result
    
    def _compute_keyword_relevance(self, query: str, title: str, description: str) -> float:
        """Compute relevance score based on keyword matching."""
        if not query:
            return 0.5
        
        query_lower = query.lower()
        stopwords = {'the', 'a', 'an', 'and', 'or', 'in', 'on', 'at', 'to', 'for', 
                    'of', 'with', 'by', 'from', 'data', 'dataset', 'study', 'analysis',
                    'metabolomics', 'metabolomic', 'using', 'based'}
        
        query_words = [w.strip() for w in query_lower.split() if len(w.strip()) > 2]
        keywords = [w for w in query_words if w not in stopwords]
        
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
    
    def _detect_assay_type(self, study: dict) -> AssayType:
        """Detect assay type from study metadata."""
        title = (study.get("title") or "").lower()
        description = (study.get("description") or "").lower()
        combined = f"{title} {description}"
        
        # MS-based techniques
        if "lc-ms" in combined or "lcms" in combined or "liquid chromatography" in combined:
            return AssayType.LC_MS
        if "gc-ms" in combined or "gcms" in combined or "gas chromatography" in combined:
            return AssayType.GC_MS
        if "ce-ms" in combined or "capillary electrophoresis" in combined:
            return AssayType.CE_MS
        if any(term in combined for term in ["mass spec", "ms", "mass-spec", "orbitrap", "qtof"]):
            return AssayType.LC_MS  # Default MS type
        
        # NMR-based techniques
        if "nmr" in combined or "nuclear magnetic resonance" in combined:
            return AssayType.NMR
        
        return AssayType.OTHER
    
    async def _fetch_with_retry(self, url: str) -> httpx.Response:
        """Fetch URL with retry logic. Delegates to the shared base helper."""
        client = await self._get_client()
        return await super()._fetch_with_retry(client, url)
    
    async def search(
        self,
        query: str,
        max_results: int = 50,
        organism: Optional[str] = None,
        **kwargs
    ) -> List[UnifiedDataset]:
        """
        Search MetaboLights for metabolomics datasets.
        Uses parallel fetching for better performance.
        
        Args:
            query: Search query
            max_results: Maximum results to return
            organism: Optional organism filter
            
        Returns:
            List of UnifiedDataset objects
            
        Raises:
            MetaboLightsError: On API or network failures
        """
        try:
            logger.info("Searching MetaboLights", query=query, max_results=max_results)
            
            # Get cached study list
            study_ids = await self._get_study_list()
            
            logger.info("MetaboLights total studies", count=len(study_ids))
            
            # Fetch details for studies and filter by query
            simplified_query = self._simplify_query(query)
            query_terms = simplified_query.lower().split()
            
            unified = []
            
            # Fetch studies in parallel batches for better performance
            max_to_check = min(len(study_ids), max_results * SEARCH_MULTIPLIER)
            
            for i in range(0, max_to_check, BATCH_SIZE):
                if len(unified) >= max_results:
                    break
                    
                batch = study_ids[i:i + BATCH_SIZE]
                
                # Fetch batch in parallel with error isolation
                tasks = [self._fetch_study_summary(sid) for sid in batch]
                results = await asyncio.gather(*tasks, return_exceptions=True)
                
                for study in results:
                    if len(unified) >= max_results:
                        break
                        
                    # Skip errors and None results
                    if study is None or isinstance(study, Exception):
                        if isinstance(study, Exception):
                            logger.debug("Study fetch failed in batch", error=str(study))
                        continue
                    
                    # Filter by organism if specified
                    if organism:
                        study_organism = study.get("organism", [])
                        if isinstance(study_organism, str):
                            study_organism = [study_organism]
                        if not any(organism.lower() in org.lower() for org in study_organism):
                            continue
                    
                    # Check if query matches
                    title = (study.get("title") or "").lower()
                    description = (study.get("description") or "").lower()
                    combined = f"{title} {description}"
                    
                    if any(term in combined for term in query_terms):
                        ud = self._to_unified(study, query=query)
                        if ud:
                            unified.append(ud)
            
            # Sort by relevance score
            unified.sort(key=lambda x: x.relevance_score, reverse=True)
            
            logger.info("MetaboLights search complete", results=len(unified))
            return unified[:max_results]
            
        except httpx.HTTPStatusError as e:
            logger.error("MetaboLights API error", status=e.response.status_code, error=str(e))
            raise MetaboLightsError(f"MetaboLights API error: {e.response.status_code}") from e
        except httpx.RequestError as e:
            logger.error("MetaboLights request failed", error=str(e))
            raise MetaboLightsError(f"MetaboLights request failed: {str(e)}") from e
    
    async def _fetch_study_summary(self, study_id: str) -> Optional[dict]:
        """
        Fetch summary information for a study.
        Returns None on error instead of raising to allow batch processing to continue.
        """
        try:
            url = f"{METABOLIGHTS_API_BASE}/studies/{study_id}"
            client = await self._get_client()
            response = await client.get(url)
            
            if response.status_code == 404:
                return None
            
            response.raise_for_status()
            data = response.json()
            
            # Handle different API response structures
            # New structure: data is nested in isaInvestigation.studies[0]
            if "isaInvestigation" in data:
                isa = data.get("isaInvestigation", {})
                studies = isa.get("studies", [])
                study_info = studies[0] if studies else {}
                
                # Extract people/contacts for organism info
                people = study_info.get("people", [])
                
                return {
                    "accession": isa.get("identifier") or study_id,
                    "title": study_info.get("title", ""),
                    "description": study_info.get("description", ""),
                    "submissionDate": study_info.get("submissionDate"),
                    "releaseDate": study_info.get("publicReleaseDate"),
                    "organism": self._extract_organisms(study_info),
                    "sampleCount": self._count_samples(study_info),
                    "studyFactors": self._extract_factors(study_info),
                    "platform": self._extract_platform(study_info),
                }
            
            # Legacy/alternative structure - fields at top level
            return {
                "accession": data.get("accession") or data.get("studyIdentifier") or study_id,
                "title": data.get("title") or data.get("studyTitle") or "",
                "description": data.get("description") or data.get("studyDescription") or "",
                "submissionDate": data.get("submissionDate"),
                "releaseDate": data.get("releaseDate"),
                "organism": data.get("organism", []),
                "sampleCount": data.get("sampleCount") or data.get("numberOfSamples") or 0,
                "studyFactors": data.get("studyFactors") or data.get("factors") or [],
                "platform": data.get("platform") or data.get("instrument"),
            }
            
        except Exception as e:
            # Return None instead of raising to allow batch processing to continue
            logger.debug("Failed to fetch study", study_id=study_id, error=str(e))
            return None
    
    def _extract_organisms(self, study_info: dict) -> List[str]:
        """Extract organism names from study info."""
        organisms = []
        
        # Try to get from characteristics in samples
        materials = study_info.get("materials", {})
        if not isinstance(materials, dict):
            materials = {}
        samples = materials.get("samples", [])
        if not isinstance(samples, list):
            return organisms

        for sample in samples[:5]:  # Check first 5 samples
            if not isinstance(sample, dict):
                continue
            characteristics = sample.get("characteristics", [])
            if not isinstance(characteristics, list):
                continue
            for char in characteristics:
                if not isinstance(char, dict):
                    continue
                category = char.get("category", {})
                if isinstance(category, dict) and category.get("annotationValue", "").lower() == "organism":
                    value_obj = char.get("value", {})
                    if isinstance(value_obj, dict):
                        org_value = value_obj.get("annotationValue", "")
                        if org_value and org_value not in organisms:
                            organisms.append(org_value)
        
        # Fallback: direct organism field if present
        if not organisms:
            org_field = study_info.get("organism") or study_info.get("organisms")
            if isinstance(org_field, list):
                for o in org_field:
                    if isinstance(o, dict):
                        val = o.get("name") or o.get("annotationValue") or ""
                        if val and val not in organisms:
                            organisms.append(val)
                    elif isinstance(o, str) and o not in organisms:
                        organisms.append(o)
            elif isinstance(org_field, str) and org_field not in organisms:
                organisms.append(org_field)

        return organisms
    
    def _extract_platform(self, study_info: dict) -> Optional[str]:
        """Extract analytical platform from study info."""
        assays = study_info.get("assays", [])
        if isinstance(assays, list) and assays:
            tech_platform = assays[0].get("technologyPlatform", "")
            if tech_platform:
                return tech_platform
            measurement_type = assays[0].get("measurementType", {}).get("annotationValue", "")
            if measurement_type:
                return measurement_type
        return None
    
    def _count_samples(self, study_info: dict) -> int:
        """Count samples in study, handling various response structures."""
        materials = study_info.get("materials", {})
        if not isinstance(materials, dict):
            return int(study_info.get("numberOfSamples") or 0)
        samples = materials.get("samples", [])
        if isinstance(samples, list):
            return len(samples)
        # Fallback to reported number if present
        count = study_info.get("numberOfSamples") or 0
        try:
            return int(count)
        except Exception:
            return 0
    
    def _extract_factors(self, study_info: dict) -> List[str]:
        """Extract study factors, handling various response structures."""
        factors = study_info.get("factors", [])
        if not isinstance(factors, list):
            return []
        result = []
        for f in factors:
            if isinstance(f, dict):
                name = f.get("factorName", "")
                if name and isinstance(name, str):
                    result.append(name)
        return result
    
    async def get_dataset(self, accession: str) -> Optional[UnifiedDataset]:
        """Fetch a single MetaboLights study by accession with file details."""
        try:
            # Normalize accession
            if not accession.upper().startswith("MTBLS"):
                accession = f"MTBLS{accession}"
            
            url = f"{METABOLIGHTS_API_BASE}/studies/{accession}"
            
            client = await self._get_client()
            response = await client.get(url)
            
            if response.status_code == 404:
                return None
            
            response.raise_for_status()
            data = response.json()
            
            # Parse the API response (same logic as _fetch_study_summary)
            if "isaInvestigation" in data:
                isa = data.get("isaInvestigation", {})
                studies = isa.get("studies", [])
                study_info = studies[0] if studies else {}
                
                study = {
                    "accession": isa.get("identifier") or accession,
                    "title": study_info.get("title", ""),
                    "description": study_info.get("description", ""),
                    "submissionDate": study_info.get("submissionDate"),
                    "releaseDate": study_info.get("publicReleaseDate"),
                    "organism": self._extract_organisms(study_info),
                    "sampleCount": self._count_samples(study_info),
                    "studyFactors": self._extract_factors(study_info),
                    "platform": self._extract_platform(study_info),
                }
            else:
                # Legacy/alternative structure
                study = {
                    "accession": data.get("accession") or data.get("studyIdentifier") or accession,
                    "title": data.get("title") or data.get("studyTitle") or "",
                    "description": data.get("description") or data.get("studyDescription") or "",
                    "submissionDate": data.get("submissionDate"),
                    "releaseDate": data.get("releaseDate"),
                    "organism": data.get("organism", []),
                    "sampleCount": data.get("sampleCount") or data.get("numberOfSamples") or 0,
                    "studyFactors": data.get("studyFactors") or data.get("factors") or [],
                    "platform": data.get("platform") or data.get("instrument"),
                }
            
            # Fetch file tree to get actual files
            file_list = await self._fetch_file_tree(accession)
            
            return self._to_unified(study, file_list=file_list)
            
        except httpx.HTTPStatusError as e:
            if e.response.status_code == 404:
                return None
            logger.error("MetaboLights API error", accession=accession, error=str(e))
            return None
        except Exception as e:
            logger.error("Failed to fetch MetaboLights study", accession=accession, error=str(e))
            return None
    
    async def _fetch_file_tree(self, accession: str) -> List[dict]:
        """Fetch file tree for a MetaboLights study."""
        try:
            url = f"{METABOLIGHTS_API_BASE}/studies/{accession}/files/tree"
            client = await self._get_client()
            response = await client.get(url)
            
            if response.status_code != 200:
                return []
            
            data = response.json()
            return data.get("study", [])
        except Exception as e:
            logger.debug("Failed to fetch file tree", accession=accession, error=str(e))
            return []
    
    def _to_unified(self, study: dict, query: Optional[str] = None, file_list: Optional[List[dict]] = None) -> Optional[UnifiedDataset]:
        """Convert a MetaboLights study to UnifiedDataset."""
        try:
            accession = study.get("accession") or study.get("studyIdentifier") or ""
            if not accession:
                return None
            
            title = study.get("title") or study.get("studyTitle") or f"MetaboLights Study {accession}"
            description = study.get("description") or study.get("studyDescription") or ""
            
            # Parse organism
            organism = study.get("organism", [])
            if isinstance(organism, str):
                organism = [organism]
            elif isinstance(organism, list):
                organism = [str(o.get("name", o) if isinstance(o, dict) else o) for o in organism]
            
            # Parse dates - keep as string for the model
            submission_date = study.get("submissionDate") or study.get("releaseDate")
            # submission_date should stay as string (model expects string, not datetime)
            
            # Sample count
            sample_count = study.get("sampleCount") or study.get("numberOfSamples") or 0
            if isinstance(sample_count, str):
                try:
                    sample_count = int(sample_count)
                except (ValueError, TypeError):
                    sample_count = 0
            
            # Detect assay type
            assay_type = self._detect_assay_type(study)
            
            # Compute relevance
            relevance = self._compute_keyword_relevance(query, title, description) if query else 0.5
            
            # Study factors (experimental conditions)
            factors = study.get("studyFactors") or study.get("factors") or []
            if isinstance(factors, list):
                factors = [str(f.get("name", f) if isinstance(f, dict) else f) for f in factors]
            
            # Platform/instrument
            platform = study.get("platform") or study.get("instrument")
            if isinstance(platform, list):
                platform = platform[0] if platform else None
            
            # Create extension
            extension = MetaboLightsExtension(
                mtbls_id=accession,
                metabolites_identified=study.get("metabolitesIdentified") or 0,
                study_design=study.get("studyDesign"),
                analytical_platform=platform,
                study_factors=factors,
            )
            
            # Generate download links (with file list if available)
            download_links = self._generate_download_links(accession, file_list=file_list)
            
            return UnifiedDataset(
                id=self.build_unified_id(accession),
                accession=accession,
                source=self.source,
                source_url=f"https://www.ebi.ac.uk/metabolights/{accession}",
                title=title,
                description=description,
                organism=organism,
                sample_count=sample_count,
                sample_count_display=str(sample_count) if sample_count else "N/A",
                submission_date=submission_date,
                omics_type=OmicsType.METABOLOMICS,
                assay_types=[assay_type] if assay_type else [],
                relevance_score=relevance,
                download_links=download_links,
                extensions={"metabolights": extension.model_dump()},
            )
            
        except Exception as e:
            logger.warning("Failed to parse MetaboLights study", error=str(e))
            return None
    
    def _create_source_url(self, accession: str) -> str:
        """Create URL to MetaboLights study page."""
        return f"https://www.ebi.ac.uk/metabolights/{accession}"
    
    def _generate_download_links(self, accession: str, file_list: Optional[List[dict]] = None) -> List[DownloadLink]:
        """
        Generate download links for a MetaboLights study.
        
        MetaboLights FTP structure:
        - ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/{MTBLS}/
        
        Args:
            accession: MTBLS accession (e.g., "MTBLS123")
            file_list: Optional list of files from the file tree API
            
        Returns:
            List of DownloadLink objects
        """
        links = []
        ftp_base = f"ftp://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/{accession}/"
        
        # If we have a file list, add individual ISA metadata files first
        if file_list:
            # Separate files and directories
            metadata_files = []
            raw_dir_exists = False
            
            for item in file_list:
                is_dir = item.get("directory", False)
                file_name = item.get("file", "")
                file_type = item.get("type", "")
                
                if is_dir:
                    if file_name.upper() == "FILES":
                        raw_dir_exists = True
                else:
                    # ISA metadata files: i_Investigation.txt, s_*.txt, a_*.txt, m_*.tsv
                    if file_type.startswith("metadata") or file_name.startswith(("i_", "s_", "a_", "m_")):
                        metadata_files.append({
                            "name": file_name,
                            "type": file_type,
                        })
            
            # Add ISA metadata files as individual download links
            for mf in metadata_files:
                fname = mf["name"]
                ftype = mf["type"]
                
                # Determine description based on file prefix
                if fname.startswith("i_"):
                    desc = "Investigation file (study metadata)"
                elif fname.startswith("s_"):
                    desc = "Sample sheet"
                elif fname.startswith("a_"):
                    desc = "Assay file"
                elif fname.startswith("m_"):
                    desc = "Metabolite assignment file (MAF)"
                else:
                    desc = "Metadata file"
                
                links.append(DownloadLink(
                    url=f"{ftp_base}{fname}",
                    file_type=DownloadFileType.METADATA,
                    file_name=fname,
                    protocol="ftp",
                    description=desc,
                ))
            
            # Add raw data directory link if it exists
            if raw_dir_exists:
                links.append(DownloadLink(
                    url=f"{ftp_base}FILES/",
                    file_type=DownloadFileType.RAW,
                    protocol="ftp",
                    description="Raw/derived data files (FTP directory)",
                ))
        
        # Always add utility links
        
        # ISA-Tab ZIP download
        links.append(DownloadLink(
            url=f"https://www.ebi.ac.uk/metabolights/ws/studies/{accession}/download",
            file_type=DownloadFileType.METADATA,
            protocol="https",
            description="ISA-Tab metadata archive (ZIP)",
        ))
        
        # HTTPS file browser (for viewing in browser)
        links.append(DownloadLink(
            url=f"https://www.ebi.ac.uk/metabolights/{accession}/files",
            file_type=DownloadFileType.OTHER,
            protocol="https",
            description="MetaboLights file browser",
        ))
        
        # FTP archive - full study
        links.append(DownloadLink(
            url=ftp_base,
            file_type=DownloadFileType.RAW,
            protocol="ftp",
            description="FTP archive (all study files)",
        ))
        
        return links


# Singleton instance
_metabolights_adapter: Optional[MetaboLightsAdapter] = None


def get_metabolights_adapter() -> MetaboLightsAdapter:
    """Get or create MetaboLights adapter instance."""
    global _metabolights_adapter
    if _metabolights_adapter is None:
        _metabolights_adapter = MetaboLightsAdapter()
    return _metabolights_adapter
