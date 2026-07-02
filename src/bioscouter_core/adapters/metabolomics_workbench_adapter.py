"""
Metabolomics Workbench Adapter
Searches NIH Metabolomics Workbench for metabolomics datasets.
API Docs: https://www.metabolomicsworkbench.org/tools/MWRestAPIv1.0.pdf
"""

import asyncio
import time
import re
from typing import List, Optional, Any, Dict, Set
from datetime import datetime

import httpx
import structlog

from .base import BaseSourceAdapter, COMMON_STOPWORDS, build_default_limits
from ..models.unified import (
    DataSource,
    OmicsType,
    AssayType,
    CurationLevel,
    UnifiedDataset,
    MetabolomicsWorkbenchExtension,
    DownloadLink,
    DownloadFileType,
)

logger = structlog.get_logger(__name__)

# Metabolomics Workbench REST API base URL
MW_API_BASE = "https://www.metabolomicsworkbench.org/rest"

# Cache settings
STUDY_LIST_CACHE_TTL = 3600  # 1 hour
MAX_STUDY_LIST_SIZE = 3000  # Limit cached study IDs

# Batch processing settings
BATCH_SIZE = 20  # Number of concurrent study fetches (configurable per plan)
SEARCH_MULTIPLIER = 4  # Check 4x max_results to account for filtering

# MW-specific stopwords
MW_STOPWORDS: Set[str] = frozenset({
    'metabolomics', 'metabolomic', 'metabolome', 'metabolite', 'metabolites',
    'mass', 'spectrometry', 'ms', 'nmr', 'profiling', 'workbench'
})


class MetabolomicsWorkbenchError(Exception):
    """Metabolomics Workbench API error."""
    pass


class MetabolomicsWorkbenchAdapter(BaseSourceAdapter):
    """
    Adapter for NIH Metabolomics Workbench.
    Uses MW REST API to search and retrieve metabolomics studies.
    Community-submitted datasets with validation.
    """
    
    def __init__(
        self,
        timeout: float = 6.0,
        batch_size: int = BATCH_SIZE,
        search_timeout: float = 10.0,
    ):
        super().__init__()
        self.timeout = timeout
        self.batch_size = batch_size
        self.search_timeout = search_timeout
        self._client: Optional[httpx.AsyncClient] = None
        
        # Cache for study list
        self._study_list_cache: List[Dict] = []
        self._study_list_cache_time: float = 0
    
    @property
    def source(self) -> DataSource:
        return DataSource.METABOLOMICS_WORKBENCH
    
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
    
    async def close(self):
        """Close HTTP client."""
        if self._client and not self._client.is_closed:
            await self._client.aclose()
            self._client = None
    
    async def _fetch_with_retry(self, url: str) -> httpx.Response:
        """Fetch URL with retry logic. Delegates to the shared base helper."""
        client = await self._get_client()
        return await super()._fetch_with_retry(client, url)
    
    async def _get_study_list_cached(self) -> List[Dict]:
        """
        Get list of all available studies with caching.
        Uses /study/study_id/ST/available endpoint.
        """
        current_time = time.time()
        
        # Return cached list if still valid
        if (self._study_list_cache and 
            current_time - self._study_list_cache_time < STUDY_LIST_CACHE_TTL):
            logger.debug("Using cached MW study list", count=len(self._study_list_cache))
            return self._study_list_cache
        
        logger.info("Fetching MW study list")
        
        # Fetch all studies using the /all endpoint
        # MW REST API format: context/input_specification/input_item/output_specification
        url = f"{MW_API_BASE}/study/study_id/ST/summary"
        response = await self._fetch_with_retry(url)
        
        if response.status_code != 200:
            logger.error("MW study list fetch failed", status=response.status_code)
            # Try alternative endpoint
            url = f"{MW_API_BASE}/study/study_id/ST"
            response = await self._fetch_with_retry(url)
            if response.status_code != 200:
                raise MetabolomicsWorkbenchError(f"Failed to fetch study list: {response.status_code}")
        
        try:
            text = response.text
            # Check if response is HTML (error page) or actual JSON
            if text.strip().startswith('<'):
                logger.warning("MW returned HTML instead of JSON, trying alternative endpoint")
                # Try different endpoint structure
                url = f"{MW_API_BASE}/study/study_title/ALL/summary"
                response = await self._fetch_with_retry(url)
                if response.status_code != 200 or response.text.strip().startswith('<'):
                    logger.error("MW API not returning JSON", response_preview=text[:200])
                    return []
                text = response.text
            
            data = response.json()
        except Exception as e:
            logger.error("MW response parse failed", error=str(e))
            raise MetabolomicsWorkbenchError(f"Failed to parse MW response: {e}")
        
        # MW returns dict with numbered keys like {"1": {...}, "2": {...}} or list
        studies = []
        if isinstance(data, list):
            studies = data
        elif isinstance(data, dict):
            # Could be {"1": study1, "2": study2} format (common for MW)
            if all(k.isdigit() for k in list(data.keys())[:5] if isinstance(k, str)):
                studies = list(data.values())
            elif "study_id" in data:
                studies = [data]
            else:
                # Try to find studies array in response
                for key in ["studies", "data", "results"]:
                    if key in data and isinstance(data[key], list):
                        studies = data[key]
                        break
                # If still empty, assume dict values are studies
                if not studies and data:
                    studies = list(data.values())
        
        # Update cache
        self._study_list_cache = studies[:MAX_STUDY_LIST_SIZE]
        self._study_list_cache_time = current_time
        
        logger.info("MW study list cached", count=len(self._study_list_cache))
        return self._study_list_cache
    
    async def _fetch_study_details(self, study_id: str) -> Optional[Dict]:
        """Fetch detailed information for a specific study."""
        try:
            url = f"{MW_API_BASE}/study/study_id/{study_id}/summary"
            response = await self._fetch_with_retry(url)
            
            if response.status_code != 200:
                logger.debug("MW study details fetch failed", study_id=study_id, status=response.status_code)
                return None
            
            data = response.json()
            
            # Handle various response formats
            if isinstance(data, list) and len(data) > 0:
                return data[0]
            elif isinstance(data, dict):
                return data
            
            return None
            
        except Exception as e:
            logger.warning("Failed to fetch MW study details", study_id=study_id, error=str(e))
            return None
    
    async def _fetch_batch_details(self, study_ids: List[str]) -> List[Dict]:
        """Fetch details for multiple studies in parallel batches."""
        all_details = []
        
        for i in range(0, len(study_ids), self.batch_size):
            batch = study_ids[i:i + self.batch_size]
            
            # Create tasks for parallel fetch
            tasks = [self._fetch_study_details(sid) for sid in batch]
            results = await asyncio.gather(*tasks, return_exceptions=True)
            
            # Collect successful results
            for result in results:
                if isinstance(result, dict) and result:
                    all_details.append(result)
                elif isinstance(result, Exception):
                    logger.debug("Batch fetch exception", error=str(result))
            
            # Brief pause between batches to avoid overwhelming the API
            if i + self.batch_size < len(study_ids):
                await asyncio.sleep(0.1)
        
        return all_details
    
    def _detect_assay_type(self, study: Dict) -> AssayType:
        """Detect assay type from study metadata."""
        analysis_type = (study.get("analysis_type") or "").lower()
        ms_type = (study.get("ms_type") or "").lower()
        chromatography = (study.get("chromatography") or "").lower()
        title = (study.get("title") or "").lower()
        description = (study.get("study_summary") or study.get("description") or "").lower()
        combined = f"{analysis_type} {ms_type} {chromatography} {title} {description}"
        
        # NMR-based techniques
        if "nmr" in combined or "nuclear magnetic resonance" in combined:
            return AssayType.NMR
        
        # MS-based techniques with chromatography type
        if "gc" in combined or "gas chromatography" in combined:
            return AssayType.GC_MS
        if "ce" in combined or "capillary electrophoresis" in combined:
            return AssayType.CE_MS
        if "lc" in combined or "liquid chromatography" in combined:
            return AssayType.LC_MS
        
        # Default to LC-MS for MS-based studies
        if "ms" in combined or "mass spec" in combined:
            return AssayType.LC_MS
        
        return AssayType.OTHER
    
    def _extract_secondary_accessions(self, study: Dict) -> List[str]:
        """Extract MetaboLights/other accessions for cross-source deduplication."""
        accessions = []
        
        # Check for MetaboLights accession
        mtbls_acc = study.get("metabo_lights") or study.get("metabolights_id")
        if mtbls_acc:
            accessions.append(mtbls_acc)
        
        # Check other fields that might contain accessions
        external_ids = study.get("external_ids") or ""
        if isinstance(external_ids, str):
            # Look for MTBLS accessions
            mtbls_matches = re.findall(r'MTBLS\d+', external_ids, re.IGNORECASE)
            accessions.extend(mtbls_matches)
            
            # Look for BioProject accessions
            bp_matches = re.findall(r'PRJNA\d+', external_ids)
            accessions.extend(bp_matches)
        
        return list(set(accessions))  # Remove duplicates
    
    def _generate_description(self, study: Dict) -> str:
        """
        Auto-generate description if not available.
        Format: "{analysis_type} metabolomics study of {species}"
        """
        # Check for existing description
        description = study.get("study_summary") or study.get("description")
        if description and len(description) > 50:
            return description
        
        # Auto-generate
        analysis_type = study.get("analysis_type") or "MS"
        species = study.get("subject_species") or study.get("species") or "unknown species"
        study_type = study.get("study_type") or "targeted"
        
        return f"{analysis_type} {study_type} metabolomics study of {species}"
    
    def _compute_keyword_relevance(self, query: str, study: Dict) -> float:
        """Compute relevance score based on keyword matching."""
        if not query:
            return 0.5
        
        query_lower = query.lower()
        
        # Filter stopwords
        stopwords = COMMON_STOPWORDS | MW_STOPWORDS
        query_words = [w.strip() for w in query_lower.split() if len(w.strip()) > 2]
        keywords = [w for w in query_words if w not in stopwords]
        
        if not keywords:
            keywords = query_words[:3]
        
        if not keywords:
            return 0.3
        
        # Get searchable text from study
        title = (study.get("title") or "").lower()
        summary = (study.get("study_summary") or study.get("description") or "").lower()
        species = (study.get("subject_species") or study.get("species") or "").lower()
        disease = (study.get("disease") or study.get("study_disease") or "").lower()
        tissue = (study.get("tissue") or study.get("sample_type") or "").lower()
        
        # Calculate matches
        total_keywords = len(keywords)
        title_matches = sum(1 for kw in keywords if kw in title)
        summary_matches = sum(1 for kw in keywords if kw in summary)
        species_matches = sum(1 for kw in keywords if kw in species)
        disease_matches = sum(1 for kw in keywords if kw in disease)
        tissue_matches = sum(1 for kw in keywords if kw in tissue)
        
        # Calculate score with weights
        title_score = (title_matches / total_keywords) * 0.35 if total_keywords > 0 else 0
        summary_score = (summary_matches / total_keywords) * 0.2 if total_keywords > 0 else 0
        species_score = (species_matches / total_keywords) * 0.15 if total_keywords > 0 else 0
        disease_score = (disease_matches / total_keywords) * 0.15 if total_keywords > 0 else 0
        tissue_score = (tissue_matches / total_keywords) * 0.15 if total_keywords > 0 else 0
        
        base_score = title_score + summary_score + species_score + disease_score + tissue_score
        
        # Boost for exact phrase match
        if query_lower in title:
            base_score = min(1.0, base_score + 0.3)
        elif query_lower in summary:
            base_score = min(1.0, base_score + 0.15)
        
        return max(0.1, min(1.0, base_score + 0.1))
    
    def _filter_matches_query(self, study: Dict, query: str) -> bool:
        """Check if study matches search query."""
        if not query:
            return True
        
        # Get all searchable text
        searchable = " ".join([
            study.get("title") or "",
            study.get("study_summary") or "",
            study.get("description") or "",
            study.get("subject_species") or "",
            study.get("species") or "",
            study.get("disease") or "",
            study.get("study_disease") or "",
            study.get("tissue") or "",
            study.get("sample_type") or "",
        ]).lower()
        
        # Extract search terms
        query_lower = query.lower()
        stopwords = COMMON_STOPWORDS | MW_STOPWORDS
        terms = [t for t in query_lower.split() if t not in stopwords and len(t) > 2]
        
        if not terms:
            terms = query_lower.split()[:3]
        
        # Require at least one term to match
        return any(term in searchable for term in terms)
    
    async def search(
        self,
        query: str,
        max_results: int = 50,
        organism: Optional[str] = None,
        **kwargs
    ) -> List[UnifiedDataset]:
        """
        Search Metabolomics Workbench for metabolomics datasets.
        
        Args:
            query: Search query
            max_results: Maximum results to return
            organism: Optional organism filter
            
        Returns:
            List of UnifiedDataset objects
        """
        try:
            logger.info("Searching MW", query=query, max_results=max_results)
            
            # Get cached study list
            started_at = asyncio.get_running_loop().time()
            study_list = await asyncio.wait_for(
                self._get_study_list_cached(),
                timeout=min(self.timeout, self.search_timeout),
            )
            
            logger.info("MW study list", count=len(study_list))
            
            # Pre-filter study list by basic query matching
            matching_studies = []
            for study in study_list:
                if self._filter_matches_query(study, query):
                    matching_studies.append(study)
                    if len(matching_studies) >= max_results * SEARCH_MULTIPLIER:
                        break
            
            logger.info("MW pre-filtered", count=len(matching_studies))
            
            # If we have few pre-filtered results, fetch details for all
            # Otherwise, limit to a reasonable number
            studies_to_fetch = matching_studies[:max_results * 2]
            
            # Fetch detailed info for matching studies in parallel batches
            study_ids = [s.get("study_id") for s in studies_to_fetch if s.get("study_id")]
            
            if study_ids:
                elapsed = asyncio.get_running_loop().time() - started_at
                remaining = max(0.1, self.search_timeout - elapsed)
                detailed_studies = await asyncio.wait_for(
                    self._fetch_batch_details(study_ids),
                    timeout=remaining,
                )
            else:
                # Use basic info if no IDs available
                detailed_studies = studies_to_fetch
            
            logger.info("MW detailed studies fetched", count=len(detailed_studies))
            
            # Convert to unified format and filter by organism
            unified = []
            for study in detailed_studies:
                try:
                    # Apply organism filter if specified
                    if organism:
                        study_organism = study.get("subject_species") or study.get("species") or ""
                        if organism.lower() not in study_organism.lower():
                            continue
                    
                    ud = self._to_unified(study, query=query)
                    if ud:
                        unified.append(ud)
                except Exception as e:
                    logger.warning("Failed to convert MW study", error=str(e))
                    continue
            
            # Sort by relevance
            unified.sort(key=lambda x: x.relevance_score, reverse=True)
            
            # Return top results
            unified = unified[:max_results]
            
            logger.info("MW search complete", results=len(unified))
            return unified
            
        except asyncio.TimeoutError:
            logger.warning(
                "MW search timed out",
                timeout_seconds=self.search_timeout,
                query=query[:100],
            )
            return []
        except httpx.HTTPStatusError as e:
            logger.error("MW API error", status=e.response.status_code, error=str(e))
            raise MetabolomicsWorkbenchError(f"MW API error: {e.response.status_code}") from e
        except httpx.RequestError as e:
            logger.error("MW request failed", error=str(e))
            raise MetabolomicsWorkbenchError(f"MW request failed: {str(e)}") from e
    
    async def get_dataset(self, accession: str) -> Optional[UnifiedDataset]:
        """Fetch a single MW study by ID."""
        try:
            study = await self._fetch_study_details(accession)
            if study:
                return self._to_unified(study)
            return None
        except asyncio.TimeoutError:
            logger.warning(
                "MW dataset fetch timed out",
                timeout_seconds=self.search_timeout,
                accession=accession,
            )
            return None
        except Exception as e:
            logger.error("Failed to fetch MW study", accession=accession, error=str(e))
            return None
    
    def _to_unified(self, study: Dict, query: Optional[str] = None) -> Optional[UnifiedDataset]:
        """Convert MW study to UnifiedDataset."""
        try:
            study_id = study.get("study_id")
            if not study_id:
                return None
            
            # Basic info
            title = study.get("title") or study.get("study_title") or f"Study {study_id}"
            description = self._generate_description(study)
            
            # Species/organism
            organism = study.get("subject_species") or study.get("species") or ""
            organisms = [organism] if organism else []
            
            # Disease
            disease = study.get("disease") or study.get("study_disease") or ""
            diseases = [disease] if disease and disease.lower() != "none" else []
            
            # Tissue/sample type
            tissue = study.get("tissue") or study.get("sample_type") or ""
            tissues = [tissue] if tissue else []
            
            # Sample count
            sample_count = 0
            try:
                sample_count = int(study.get("subject_count") or study.get("num_subjects") or 0)
            except (ValueError, TypeError):
                pass
            
            # Metabolite counts
            metabolite_count = 0
            named_metabolite_count = 0
            try:
                metabolite_count = int(study.get("num_metabolites") or study.get("metabolite_count") or 0)
                named_metabolite_count = int(study.get("num_named_metabolites") or study.get("named_metabolites") or 0)
            except (ValueError, TypeError):
                pass
            
            # Assay type detection
            assay_type = self._detect_assay_type(study)
            
            # Dates
            submission_date = study.get("submit_date") or study.get("created")
            release_date = study.get("release_date") or study.get("released")
            last_update = study.get("last_updated") or study.get("updated")
            
            # Contributors
            pi = study.get("pi") or study.get("principal_investigator")
            contributors = [pi] if pi else []
            
            # Institution
            institute = study.get("institute") or study.get("institution")
            
            # Analysis info
            analysis_type = study.get("analysis_type") or ""
            ms_type = study.get("ms_type") or ""
            ion_mode = study.get("ion_mode") or ""
            chromatography = study.get("chromatography") or ""
            study_type = study.get("study_type") or ""
            
            # Factors
            factors = []
            factor_data = study.get("factors") or study.get("study_factors")
            if isinstance(factor_data, list):
                factors = factor_data
            elif isinstance(factor_data, str) and factor_data:
                factors = [f.strip() for f in factor_data.split(",")]
            
            # Secondary accessions (for cross-source dedup with MetaboLights)
            secondary_accessions = self._extract_secondary_accessions(study)
            
            # Calculate relevance
            relevance_score = self._compute_keyword_relevance(query, study) if query else 0.5
            
            # Build source URL
            source_url = f"https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID={study_id}"
            
            # Build download links
            download_links = [
                DownloadLink(
                    url=source_url,
                    file_type=DownloadFileType.OTHER,
                    file_name=None,
                    file_size_bytes=None,
                    protocol="https",
                    description="View study and download data"
                )
            ]
            
            # Add data download link if available
            data_download_url = f"https://www.metabolomicsworkbench.org/data/DRCCData.php?Mode=StudyData&StudyID={study_id}"
            download_links.append(DownloadLink(
                url=data_download_url,
                file_type=DownloadFileType.PROCESSED,
                file_name=None,
                file_size_bytes=None,
                protocol="https",
                description="Download study data"
            ))
            
            # Build match reasons
            match_reasons = []
            if query:
                query_lower = query.lower()
                if query_lower in title.lower():
                    match_reasons.append(f"Title contains query")
                if organism and any(t in organism.lower() for t in query_lower.split()):
                    match_reasons.append(f"Species: {organism}")
                if disease and any(t in disease.lower() for t in query_lower.split()):
                    match_reasons.append(f"Disease: {disease}")
                if tissue and any(t in tissue.lower() for t in query_lower.split()):
                    match_reasons.append(f"Tissue: {tissue}")
            
            if not match_reasons:
                match_reasons.append("General match")
            
            # Build extension
            extension = MetabolomicsWorkbenchExtension(
                study_id=study_id,
                project_id=study.get("project_id"),
                institute=institute,
                analysis_type=analysis_type,
                ms_type=ms_type,
                ion_mode=ion_mode,
                chromatography=chromatography,
                metabolite_count=metabolite_count,
                named_metabolite_count=named_metabolite_count,
                study_type=study_type,
                study_status=study.get("study_status") or "public",
                subject_count=sample_count,
                factors=factors,
                metaboanalyst_link=study.get("metaboanalyst_link"),
                metabominer_link=study.get("metabominer_link")
            )
            
            # Build unified dataset
            return UnifiedDataset(
                id=f"mw:{study_id}",
                accession=study_id,
                source=DataSource.METABOLOMICS_WORKBENCH,
                source_url=source_url,
                secondary_accession=secondary_accessions,
                omics_type=OmicsType.METABOLOMICS,
                assay_types=[assay_type],
                title=title,
                description=description[:2000] if description else None,
                organism=organisms,
                sample_count=sample_count,
                sample_count_display=f"{sample_count:,}" if sample_count else "N/A",
                disease=diseases,
                tissue=tissues,
                cell_line=[],
                submission_date=submission_date,
                release_date=release_date,
                last_update=last_update,
                pubmed_ids=[],
                doi=study.get("doi"),
                citation=None,
                contributors=contributors,
                institution=institute,
                relevance_score=relevance_score,
                match_reasons=match_reasons[:5],
                quality_score=0.7,  # Community-submitted, moderate quality
                download_links=download_links,
                curation_level=CurationLevel.COMMUNITY,  # Community-submitted
                merged_sources=[],
                extensions={"mw": extension.model_dump()}
            )
            
        except Exception as e:
            logger.warning("Failed to convert MW study to unified", error=str(e))
            return None
    
    async def prewarm_cache(self) -> bool:
        """
        Pre-warm the study list cache.
        Called during startup to ensure cache is ready.
        """
        try:
            logger.info("Pre-warming MW study list cache")
            await asyncio.wait_for(
                self._get_study_list_cached(),
                timeout=self.search_timeout,
            )
            logger.info("MW cache pre-warmed successfully", count=len(self._study_list_cache))
            return True
        except Exception as e:
            logger.error("Failed to pre-warm MW cache", error=str(e))
            return False


# Singleton instance
_mw_adapter: Optional[MetabolomicsWorkbenchAdapter] = None


def get_metabolomics_workbench_adapter(batch_size: int = BATCH_SIZE) -> MetabolomicsWorkbenchAdapter:
    """Get or create Metabolomics Workbench adapter singleton."""
    global _mw_adapter
    if _mw_adapter is None:
        _mw_adapter = MetabolomicsWorkbenchAdapter(batch_size=batch_size)
    return _mw_adapter
