"""
CPTAC/PDC (Proteomic Data Commons) Adapter for BioScouter.

Searches the NCI Proteomic Data Commons for CPTAC and other cancer proteomics studies.
Uses GraphQL API for efficient querying.

API Documentation:
- GraphQL Endpoint: https://pdc.cancer.gov/graphql
- GitHub: https://github.com/esacinc/PDC-Public
"""

import asyncio
from datetime import date, datetime
from typing import Optional, List, Any, Tuple

import httpx
import structlog

from .base import build_default_limits
from bioscouter_core.models.unified import (
    UnifiedDataset,
    DataSource,
    OmicsType,
    AssayType,
    CurationLevel,
    CPTACExtension,
)

logger = structlog.get_logger(__name__)


# === CONFIGURATION ===

PDC_GRAPHQL_URL = "https://pdc.cancer.gov/graphql"
DEFAULT_TIMEOUT = 60.0  # 60 seconds for GraphQL queries
MAX_RESULTS_PER_PAGE = 100


# === GRAPHQL QUERIES ===

# Query for paginated study search
STUDY_SEARCH_QUERY = """
query getPaginatedUIStudy(
    $disease_type: String,
    $primary_site: String,
    $analytical_fraction: String,
    $experiment_type: String,
    $program_name: String,
    $offset: Int!,
    $limit: Int!
) {
    getPaginatedUIStudy(
        disease_type: $disease_type,
        primary_site: $primary_site,
        analytical_fraction: $analytical_fraction,
        experiment_type: $experiment_type,
        program_name: $program_name,
        offset: $offset,
        limit: $limit
    ) {
        total
        uiStudies {
            study_id
            pdc_study_id
            study_submitter_id
            submitter_id_name
            program_name
            project_name
            disease_type
            primary_site
            analytical_fraction
            experiment_type
            embargo_date
            cases_count
            aliquots_count
            filesCount {
                data_category
                file_type
                files_count
            }
        }
        pagination {
            count
            from
            page
            total
            pages
            size
        }
    }
}
"""

# Query for text search by study name
TEXT_SEARCH_QUERY = """
query studySearch($name: String!, $offset: Int!, $limit: Int!) {
    studySearch(name: $name, offset: $offset, limit: $limit) {
        total
        studies {
            record_type
            name
            submitter_id_name
            study_id
            study_submitter_id
            pdc_study_id
        }
        pagination {
            count
            page
            total
            pages
            size
        }
    }
}
"""

# Query for detailed study info
STUDY_DETAIL_QUERY = """
query study($pdc_study_id: String!, $acceptDUA: Boolean) {
    study(pdc_study_id: $pdc_study_id, acceptDUA: $acceptDUA) {
        study_id
        study_submitter_id
        pdc_study_id
        study_name
        study_shortname
        study_description
        program_name
        project_name
        disease_type
        primary_site
        analytical_fraction
        experiment_type
        cases_count
        aliquots_count
        embargo_date
        filesCount {
            data_category
            file_type
            files_count
        }
    }
}
"""

# Query for available filters (programs, diseases, sites, etc.)
FILTERS_QUERY = """
query uiFilters {
    uiFilters {
        program_name { filterName filterValue }
        project_name { filterName filterValue }
        disease_type { filterName filterValue }
        primary_site { filterName filterValue }
        analytical_fraction { filterName filterValue }
        experiment_type { filterName filterValue }
        acquisition_type { filterName filterValue }
        sample_type { filterName filterValue }
    }
}
"""


class CPTACSearchError(Exception):
    """Custom exception for CPTAC search errors."""
    pass


class CPTACAdapter:
    """Adapter for searching CPTAC/PDC proteomics data."""
    
    # Source identifier
    source = DataSource.CPTAC
    source_name = "CPTAC/PDC"
    
    def __init__(self):
        self._client: Optional[httpx.AsyncClient] = None
        # Filter cache
        self._filters_cache: Optional[dict] = None
        self._filters_cache_time: float = 0
        self._filters_cache_ttl: float = 3600 * 6  # 6 hours
        # Study list cache for text search
        self._studies_cache: Optional[List[dict]] = None
        self._studies_cache_time: float = 0
        self._studies_cache_ttl: float = 3600  # 1 hour
    
    async def _get_client(self) -> httpx.AsyncClient:
        """Get or create HTTP client."""
        if self._client is None or self._client.is_closed:
            self._client = httpx.AsyncClient(
                timeout=httpx.Timeout(DEFAULT_TIMEOUT),
                limits=build_default_limits(),
                headers={
                    "Content-Type": "application/json",
                    "Accept": "application/json",
                }
            )
        return self._client
    
    async def _execute_graphql(
        self, 
        query: str, 
        variables: dict,
        accept_partial: bool = True
    ) -> Tuple[Optional[dict], List[str]]:
        """
        Execute a GraphQL query against PDC.
        
        Returns:
            Tuple of (data dict or None, list of warning messages)
        """
        client = await self._get_client()
        warnings = []
        
        try:
            response = await client.post(
                PDC_GRAPHQL_URL,
                json={"query": query, "variables": variables}
            )
            response.raise_for_status()
            result = response.json()
            
            # Check for GraphQL errors
            if "errors" in result:
                error_messages = [e.get("message", str(e)) for e in result["errors"]]
                
                if accept_partial and "data" in result and result["data"]:
                    # Accept partial data with warnings
                    warnings.extend([f"GraphQL warning: {msg}" for msg in error_messages])
                    logger.warning("CPTAC GraphQL partial response", errors=error_messages)
                else:
                    # No usable data
                    raise CPTACSearchError(f"GraphQL errors: {'; '.join(error_messages)}")
            
            return result.get("data"), warnings
            
        except httpx.TimeoutException:
            raise CPTACSearchError("PDC API timeout - try a more specific query")
        except httpx.HTTPStatusError as e:
            raise CPTACSearchError(f"PDC API error: {e.response.status_code}")
        except Exception as e:
            raise CPTACSearchError(f"PDC query failed: {str(e)}")
    
    async def get_filters(self) -> dict:
        """
        Get available filter options (programs, diseases, sites, etc.).
        Results are cached for 6 hours.
        """
        now = asyncio.get_event_loop().time()
        
        if self._filters_cache and (now - self._filters_cache_time) < self._filters_cache_ttl:
            return self._filters_cache
        
        logger.info("Fetching CPTAC filter options")
        
        data, warnings = await self._execute_graphql(FILTERS_QUERY, {})
        
        if not data or "uiFilters" not in data:
            logger.warning("Failed to fetch CPTAC filters")
            return {}
        
        filters = data["uiFilters"]
        
        # Transform to more usable format
        self._filters_cache = {
            "programs": [f["filterValue"] for f in filters.get("program_name", [])],
            "diseases": [f["filterValue"] for f in filters.get("disease_type", [])],
            "sites": [f["filterValue"] for f in filters.get("primary_site", [])],
            "fractions": [f["filterValue"] for f in filters.get("analytical_fraction", [])],
            "experiment_types": [f["filterValue"] for f in filters.get("experiment_type", [])],
        }
        self._filters_cache_time = now
        
        logger.info("CPTAC filters cached", 
                   programs=len(self._filters_cache["programs"]),
                   diseases=len(self._filters_cache["diseases"]))
        
        return self._filters_cache
    
    async def _get_all_studies(self) -> List[dict]:
        """
        Fetch all studies for client-side text filtering.
        Cached for 1 hour.
        """
        now = asyncio.get_event_loop().time()
        
        if self._studies_cache and (now - self._studies_cache_time) < self._studies_cache_ttl:
            return self._studies_cache
        
        logger.info("Fetching all CPTAC studies for cache")
        
        all_studies = []
        offset = 0
        limit = 100
        
        while True:
            data, _ = await self._execute_graphql(
                STUDY_SEARCH_QUERY,
                {"offset": offset, "limit": limit}
            )
            
            if not data or "getPaginatedUIStudy" not in data:
                break
            
            result = data["getPaginatedUIStudy"]
            studies = result.get("uiStudies", [])
            
            if not studies:
                break
            
            all_studies.extend(studies)
            
            # Check if we have all
            total = result.get("total", 0)
            if len(all_studies) >= total:
                break
            
            offset += limit
            
            # Safety limit
            if offset > 1000:
                break
        
        self._studies_cache = all_studies
        self._studies_cache_time = now
        
        logger.info("CPTAC studies cached", count=len(all_studies))
        
        return all_studies
    
    def _map_experiment_type_to_assay(self, experiment_type: Optional[str]) -> List[AssayType]:
        """Map CPTAC experiment type to BioScouter AssayType."""
        if not experiment_type:
            return []
        
        exp_lower = experiment_type.lower()
        
        if "tmt" in exp_lower:
            return [AssayType.TMT]
        elif "itraq" in exp_lower:
            return [AssayType.ITRAQ]
        elif "silac" in exp_lower:
            return [AssayType.SILAC]
        elif "label" in exp_lower and "free" in exp_lower:
            return [AssayType.LABEL_FREE]
        elif "dia" in exp_lower:
            return [AssayType.DIA]
        elif "dda" in exp_lower:
            return [AssayType.DDA]
        
        return []
    
    def _check_embargo(self, embargo_date: Optional[str]) -> bool:
        """Check if a study is currently embargoed."""
        if not embargo_date:
            return False
        
        try:
            embargo = datetime.strptime(embargo_date, "%Y-%m-%d").date()
            return embargo > date.today()
        except (ValueError, TypeError):
            return False
    
    def _study_matches_query(self, study: dict, query_terms: List[str]) -> bool:
        """Check if a study matches query terms."""
        # Fields to search
        searchable = " ".join([
            study.get("submitter_id_name", "") or "",
            study.get("study_description", "") or "",
            study.get("disease_type", "") or "",
            study.get("primary_site", "") or "",
            study.get("program_name", "") or "",
            study.get("project_name", "") or "",
            study.get("analytical_fraction", "") or "",
            study.get("experiment_type", "") or "",
        ]).lower()
        
        # Check if ANY query term matches (more permissive search)
        return any(term in searchable for term in query_terms)
    
    def _study_to_dataset(self, study: dict) -> UnifiedDataset:
        """Convert a CPTAC study to UnifiedDataset."""
        pdc_id = study.get("pdc_study_id", "")
        study_name = study.get("submitter_id_name", "") or study.get("study_name", "")
        description = study.get("study_description")
        
        # If no description, generate one
        if not description:
            parts = []
            if study.get("disease_type"):
                parts.append(study["disease_type"])
            if study.get("analytical_fraction"):
                parts.append(study["analytical_fraction"].lower())
            if study.get("experiment_type"):
                parts.append(f"using {study['experiment_type']}")
            if study.get("program_name"):
                parts.append(f"from {study['program_name']}")
            description = " ".join(parts) + " proteomics study" if parts else None
        
        # Determine omics type based on analytical fraction
        fraction = study.get("analytical_fraction", "").lower()
        if "phospho" in fraction:
            omics_type = OmicsType.PROTEOMICS  # Could add specific phospho handling
        elif "glyco" in fraction:
            omics_type = OmicsType.PROTEOMICS
        else:
            omics_type = OmicsType.PROTEOMICS
        
        # Check embargo status
        embargo_date = study.get("embargo_date")
        is_embargoed = self._check_embargo(embargo_date)
        
        # Build file counts list
        file_counts = []
        for fc in study.get("filesCount", []) or []:
            if fc.get("files_count", 0) > 0:
                file_counts.append({
                    "data_category": fc.get("data_category"),
                    "file_type": fc.get("file_type"),
                    "count": fc.get("files_count", 0)
                })
        
        # Get all analytical fractions (may be multiple)
        fractions = [study.get("analytical_fraction")] if study.get("analytical_fraction") else []
        
        # Build extension
        extension = CPTACExtension(
            pdc_study_id=pdc_id,
            study_submitter_id=study.get("study_submitter_id"),
            program_name=study.get("program_name", "Unknown"),
            project_name=study.get("project_name"),
            analytical_fraction=study.get("analytical_fraction"),
            experiment_type=study.get("experiment_type"),
            embargo_date=embargo_date,
            is_embargoed=is_embargoed,
            aliquots_count=study.get("aliquots_count", 0) or 0,
            cases_count=study.get("cases_count", 0) or 0,
            file_counts=file_counts,
            analytical_fractions=fractions,
        )
        
        # Disease and tissue
        disease = [study["disease_type"]] if study.get("disease_type") else []
        tissue = [study["primary_site"]] if study.get("primary_site") else []
        
        return UnifiedDataset(
            id=f"cptac:{pdc_id}",
            accession=pdc_id,
            source=DataSource.CPTAC,
            source_url=f"https://pdc.cancer.gov/pdc/study/{pdc_id}",
            omics_type=omics_type,
            assay_types=self._map_experiment_type_to_assay(study.get("experiment_type")),
            title=study_name,
            description=description,
            organism=["Homo sapiens"],  # CPTAC is human cancer focused
            sample_count=study.get("cases_count", 0) or 0,
            sample_count_display=str(study.get("cases_count", 0) or 0),
            disease=disease,
            tissue=tissue,
            curation_level=CurationLevel.CURATED,  # CPTAC data is highly curated
            extensions={"cptac": extension.model_dump()},
        )
    
    async def search(
        self,
        query: str,
        max_results: int = 50,
        organism: Optional[str] = None,  # Accepted but not used (all CPTAC is human)
        disease_type: Optional[str] = None,
        primary_site: Optional[str] = None,
        analytical_fraction: Optional[str] = None,
        experiment_type: Optional[str] = None,
        program_name: Optional[str] = None,
        include_embargoed: bool = True,
    ) -> List[UnifiedDataset]:
        """
        Search CPTAC/PDC for proteomics studies.
        
        Args:
            query: Natural language search query
            max_results: Maximum number of results
            organism: Organism filter (accepted but ignored - all CPTAC is human)
            disease_type: Filter by cancer type (e.g., "Breast Cancer")
            primary_site: Filter by anatomical site (e.g., "Breast")
            analytical_fraction: Filter by fraction (Proteome, Phosphoproteome, etc.)
            experiment_type: Filter by experiment (TMT10, TMT11, etc.)
            program_name: Filter by program (CPTAC, TCGA, etc.)
            include_embargoed: Whether to include embargoed studies (with badge)
        
        Returns:
            List of UnifiedDataset results
        """
        logger.info("Searching CPTAC", 
                   query=query, 
                   max_results=max_results,
                   disease_type=disease_type,
                   program_name=program_name)
        
        # Extract search terms from query
        query_terms = [t.lower() for t in query.split() if len(t) > 2]
        
        # First, try filter-based search if we have specific filters
        has_filters = any([disease_type, primary_site, analytical_fraction, 
                         experiment_type, program_name])
        
        datasets = []
        
        if has_filters:
            # Use filter-based GraphQL query
            variables = {
                "offset": 0,
                "limit": min(max_results * 2, MAX_RESULTS_PER_PAGE),  # Fetch extra for filtering
            }
            if disease_type:
                variables["disease_type"] = disease_type
            if primary_site:
                variables["primary_site"] = primary_site
            if analytical_fraction:
                variables["analytical_fraction"] = analytical_fraction
            if experiment_type:
                variables["experiment_type"] = experiment_type
            if program_name:
                variables["program_name"] = program_name
            
            data, warnings = await self._execute_graphql(STUDY_SEARCH_QUERY, variables)
            
            if data and "getPaginatedUIStudy" in data:
                studies = data["getPaginatedUIStudy"].get("uiStudies", [])
                
                for study in studies:
                    # Apply text filter
                    if query_terms and not self._study_matches_query(study, query_terms):
                        continue
                    
                    dataset = self._study_to_dataset(study)
                    
                    # Skip embargoed if not included
                    if not include_embargoed and dataset.extensions.get("cptac", {}).get("is_embargoed"):
                        continue
                    
                    datasets.append(dataset)
                    
                    if len(datasets) >= max_results:
                        break
        else:
            # Text-based search: get all studies and filter locally
            all_studies = await self._get_all_studies()
            
            for study in all_studies:
                if not self._study_matches_query(study, query_terms):
                    continue
                
                dataset = self._study_to_dataset(study)
                
                # Skip embargoed if not included
                if not include_embargoed and dataset.extensions.get("cptac", {}).get("is_embargoed"):
                    continue
                
                datasets.append(dataset)
                
                if len(datasets) >= max_results:
                    break
        
        # Calculate relevance scores based on query match
        for ds in datasets:
            score = 0.0
            title_lower = (ds.title or "").lower()
            desc_lower = (ds.description or "").lower()
            
            for term in query_terms:
                if term in title_lower:
                    score += 0.3
                if term in desc_lower:
                    score += 0.1
                if ds.disease and any(term in d.lower() for d in ds.disease):
                    score += 0.2
                if ds.tissue and any(term in t.lower() for t in ds.tissue):
                    score += 0.15
            
            ds.relevance_score = min(score, 1.0)
        
        # Sort by relevance
        datasets.sort(key=lambda x: x.relevance_score, reverse=True)
        
        logger.info("CPTAC search complete", results=len(datasets))
        
        return datasets[:max_results]
    
    async def get_study_details(self, pdc_study_id: str) -> Optional[UnifiedDataset]:
        """Get detailed information for a specific study."""
        logger.info("Fetching CPTAC study details", pdc_study_id=pdc_study_id)
        
        data, warnings = await self._execute_graphql(
            STUDY_DETAIL_QUERY,
            {"pdc_study_id": pdc_study_id, "acceptDUA": True}
        )
        
        if not data or "study" not in data or not data["study"]:
            return None
        
        # The API returns a list, get the first item
        study = data["study"]
        if isinstance(study, list):
            if not study:
                return None
            study = study[0]
        
        return self._study_to_dataset(study)
    
    async def get_dataset(self, accession: str) -> Optional[UnifiedDataset]:
        """
        Fetch a single CPTAC/PDC dataset by accession.
        
        Args:
            accession: PDC study ID (e.g., "PDC000123")
            
        Returns:
            UnifiedDataset if found, None otherwise
        """
        try:
            # Normalize accession - PDC IDs look like "PDC000123"
            accession = accession.upper()
            if not accession.startswith("PDC"):
                accession = f"PDC{accession}"
            
            logger.info("Fetching CPTAC dataset", accession=accession)
            
            # Use the existing get_study_details method
            dataset = await self.get_study_details(accession)
            
            if dataset:
                logger.info("CPTAC dataset fetched", accession=accession)
            else:
                logger.warning("CPTAC dataset not found", accession=accession)
            
            return dataset
            
        except Exception as e:
            logger.error("Failed to fetch CPTAC dataset", 
                        accession=accession, error=str(e))
            return None
    
    async def close(self):
        """Close the HTTP client."""
        if self._client and not self._client.is_closed:
            await self._client.aclose()
            self._client = None


# === SINGLETON ACCESS ===

_adapter_instance: Optional[CPTACAdapter] = None


def get_cptac_adapter() -> CPTACAdapter:
    """Get singleton CPTAC adapter instance."""
    global _adapter_instance
    if _adapter_instance is None:
        _adapter_instance = CPTACAdapter()
    return _adapter_instance


# === ASYNC GENERATOR FOR STREAMING ===

async def search_cptac_streaming(
    query: str,
    max_results: int = 50,
    **kwargs
):
    """
    Stream CPTAC search results.
    Yields results as they become available.
    """
    adapter = get_cptac_adapter()
    results = await adapter.search(query, max_results, **kwargs)
    
    for result in results:
        yield result
