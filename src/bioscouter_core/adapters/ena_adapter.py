"""
ENA (European Nucleotide Archive) Adapter for BioScouter.

Searches ENA for sequencing studies including genomics, transcriptomics, and metagenomics.
Uses the ENA Portal API for efficient text search.

API Documentation:
- Portal API: https://www.ebi.ac.uk/ena/portal/api/
- Browser API: https://www.ebi.ac.uk/ena/browser/api/
"""

import asyncio
from typing import Optional, List, Any

import httpx
import structlog

from .base import build_default_limits
from bioscouter_core.models.unified import (
    UnifiedDataset,
    DataSource,
    OmicsType,
    AssayType,
    CurationLevel,
    ENAExtension,
)

logger = structlog.get_logger(__name__)


# === CONFIGURATION ===

ENA_PORTAL_API_URL = "https://www.ebi.ac.uk/ena/portal/api"
DEFAULT_TIMEOUT = 60.0  # 60 seconds as per user request
MAX_RESULTS_PER_PAGE = 100
RATE_LIMIT_DELAY = 0.02  # 50 req/sec = 20ms between requests

ENA_QUERY_STOPWORDS = {
    "and", "data", "dataset", "datasets", "human", "mouse", "rna", "seq",
    "rna-seq", "sequencing", "study", "studies", "transcriptomics",
}


# === LIBRARY STRATEGY TO OMICS/ASSAY MAPPING ===

LIBRARY_STRATEGY_TO_OMICS = {
    # Transcriptomics
    "RNA-Seq": OmicsType.TRANSCRIPTOMICS,
    "ssRNA-seq": OmicsType.TRANSCRIPTOMICS,
    "miRNA-Seq": OmicsType.TRANSCRIPTOMICS,
    "ncRNA-Seq": OmicsType.TRANSCRIPTOMICS,
    "EST": OmicsType.TRANSCRIPTOMICS,
    "FL-cDNA": OmicsType.TRANSCRIPTOMICS,
    
    # Genomics
    "WGS": OmicsType.GENOMICS,
    "WXS": OmicsType.GENOMICS,
    "WCS": OmicsType.GENOMICS,
    "CLONE": OmicsType.GENOMICS,
    "POOLCLONE": OmicsType.GENOMICS,
    "Targeted-Capture": OmicsType.GENOMICS,
    "RAD-Seq": OmicsType.GENOMICS,
    
    # Epigenomics
    "ChIP-Seq": OmicsType.EPIGENOMICS,
    "ATAC-seq": OmicsType.EPIGENOMICS,
    "DNase-Hypersensitivity": OmicsType.EPIGENOMICS,
    "Bisulfite-Seq": OmicsType.EPIGENOMICS,
    "MeDIP-Seq": OmicsType.EPIGENOMICS,
    "MRE-Seq": OmicsType.EPIGENOMICS,
    "MBD-Seq": OmicsType.EPIGENOMICS,
    "FAIRE-seq": OmicsType.EPIGENOMICS,
    "Hi-C": OmicsType.EPIGENOMICS,
    
    # Metagenomics
    "AMPLICON": OmicsType.METAGENOMICS,
    "WGA": OmicsType.METAGENOMICS,  # Often used in metagenomics
    
    # Other
    "OTHER": OmicsType.GENOMICS,
    "SYNTHETIC-LONG-READ": OmicsType.GENOMICS,
}

LIBRARY_STRATEGY_TO_ASSAY = {
    "RNA-Seq": AssayType.RNA_SEQ,
    "ssRNA-seq": AssayType.RNA_SEQ,
    "miRNA-Seq": AssayType.RNA_SEQ,
    "WGS": AssayType.WGS,
    "WXS": AssayType.WES,
    "ChIP-Seq": AssayType.CHIP_SEQ,
    "ATAC-seq": AssayType.ATAC_SEQ,
    "DNase-Hypersensitivity": AssayType.DNASE_SEQ,
    "Bisulfite-Seq": AssayType.BISULFITE_SEQ,
    "AMPLICON": AssayType.AMPLICON_16S,  # Default to 16S, may need refinement
    "Hi-C": AssayType.HI_C,
    "Targeted-Capture": AssayType.TARGETED_SEQ,
}


# === TAXONOMY MAPPING ===

COMMON_TAX_IDS = {
    "human": 9606,
    "homo sapiens": 9606,
    "mouse": 10090,
    "mus musculus": 10090,
    "rat": 10116,
    "rattus norvegicus": 10116,
    "zebrafish": 7955,
    "danio rerio": 7955,
    "fruit fly": 7227,
    "drosophila melanogaster": 7227,
    "c. elegans": 6239,
    "caenorhabditis elegans": 6239,
    "yeast": 4932,
    "saccharomyces cerevisiae": 4932,
    "arabidopsis": 3702,
    "arabidopsis thaliana": 3702,
    "chicken": 9031,
    "gallus gallus": 9031,
    "pig": 9823,
    "sus scrofa": 9823,
    "cow": 9913,
    "bos taurus": 9913,
    "dog": 9615,
    "canis familiaris": 9615,
}


class ENASearchError(Exception):
    """Custom exception for ENA search errors."""
    pass


class ENAAdapter:
    """Adapter for searching European Nucleotide Archive."""
    
    # Source identifier
    source = DataSource.ENA
    source_name = "ENA"
    
    def __init__(self):
        self._client: Optional[httpx.AsyncClient] = None
        self._last_request_time: float = 0
    
    async def _get_client(self) -> httpx.AsyncClient:
        """Get or create HTTP client."""
        if self._client is None or self._client.is_closed:
            self._client = httpx.AsyncClient(
                timeout=httpx.Timeout(DEFAULT_TIMEOUT),
                limits=build_default_limits(),
                headers={
                    "Accept": "application/json",
                }
            )
        return self._client
    
    async def _rate_limit(self):
        """Apply rate limiting (50 req/sec)."""
        now = asyncio.get_event_loop().time()
        elapsed = now - self._last_request_time
        if elapsed < RATE_LIMIT_DELAY:
            await asyncio.sleep(RATE_LIMIT_DELAY - elapsed)
        self._last_request_time = asyncio.get_event_loop().time()
    
    def _get_taxonomy_filter(self, organism: Optional[str]) -> Optional[str]:
        """
        Get taxonomy filter for ENA query.
        Uses tax_tree() for better coverage of subspecies.
        """
        if not organism:
            return None
        
        organism_lower = organism.lower().strip()
        
        # Check common names
        tax_id = COMMON_TAX_IDS.get(organism_lower)
        
        if tax_id:
            # Use tax_tree for hierarchical search (includes subspecies)
            return f"tax_tree({tax_id})"
        
        # Try partial match
        for name, tid in COMMON_TAX_IDS.items():
            if organism_lower in name or name in organism_lower:
                return f"tax_tree({tid})"
        
        # Fall back to text search on scientific_name
        return f'scientific_name="*{organism}*"'
    
    def _build_search_query(
        self,
        query: str,
        organism: Optional[str] = None,
        library_strategy: Optional[str] = None,
    ) -> str:
        """Build ENA search query string."""
        parts = []
        
        # Search significant terms across title and description. Requiring
        # every natural-language token in the title produced frequent zeros.
        if query:
            terms = []
            for raw_term in query.replace("-", " ").split():
                term = raw_term.strip().lower().replace('"', "")
                if len(term) > 2 and term not in ENA_QUERY_STOPWORDS:
                    terms.append(term)

            text_clauses = []
            for term in terms[:5]:
                text_clauses.extend([
                    f'study_title="*{term}*"',
                    f'study_description="*{term}*"',
                ])
            if text_clauses:
                parts.append(f"({' OR '.join(text_clauses)})")
        
        # Organism filter
        tax_filter = self._get_taxonomy_filter(organism)
        if tax_filter:
            parts.append(tax_filter)
        
        # Library strategy filter
        if library_strategy:
            parts.append(f'library_strategy="{library_strategy}"')
        
        # Join with AND
        if parts:
            return " AND ".join(parts)
        
        return ""

    def _infer_organism(self, query: str) -> Optional[str]:
        query_lower = query.lower()
        for common_name in sorted(COMMON_TAX_IDS, key=len, reverse=True):
            if common_name in query_lower:
                return common_name
        return None
    
    async def _search_studies(
        self,
        query: str,
        max_results: int,
        organism: Optional[str] = None,
        library_strategy: Optional[str] = None,
    ) -> List[dict]:
        """
        Search ENA for studies.
        """
        await self._rate_limit()
        
        client = await self._get_client()
        
        # Build query
        search_query = self._build_search_query(query, organism, library_strategy)
        
        # Fields to retrieve
        fields = [
            "study_accession",
            "secondary_study_accession", 
            "study_title",
            "study_description",
            "scientific_name",
            "tax_id",
            "center_name",
            "broker_name",
            "first_public",
            "last_updated",
            "study_alias",
        ]
        
        params = {
            "result": "study",
            "query": search_query if search_query else "*",
            "fields": ",".join(fields),
            "format": "json",
            "limit": min(max_results, MAX_RESULTS_PER_PAGE),
        }
        
        try:
            response = await client.get(
                f"{ENA_PORTAL_API_URL}/search",
                params=params
            )
            response.raise_for_status()
            
            # ENA returns empty array if no results
            data = response.json()
            
            if isinstance(data, list):
                return data
            elif isinstance(data, dict) and "error" in data:
                logger.warning("ENA search error", error=data.get("error"))
                return []
            
            return []
            
        except httpx.TimeoutException:
            raise ENASearchError("ENA API timeout - try a more specific query")
        except httpx.HTTPStatusError as e:
            if e.response.status_code == 204:
                # No content = no results
                return []
            raise ENASearchError(f"ENA API error: {e.response.status_code}")
        except Exception as e:
            logger.error("ENA search failed", error=str(e))
            raise ENASearchError(f"ENA query failed: {str(e)}")
    
    async def _get_study_samples(self, study_accession: str) -> dict:
        """Get sample and run counts for a study."""
        await self._rate_limit()
        
        client = await self._get_client()
        
        try:
            # Query for read_study to get sample/run info
            # Use secondary_study_accession for ERP accessions or study_accession for PRJ
            query_field = "secondary_study_accession" if study_accession.startswith("ERP") else "study_accession"
            response = await client.get(
                f"{ENA_PORTAL_API_URL}/search",
                params={
                    "result": "read_study",
                    "query": f'{query_field}="{study_accession}"',
                    "fields": "study_accession,sample_count,run_count,base_count",
                    "format": "json",
                    "limit": 1,
                }
            )
            
            if response.status_code == 204:
                return {"sample_count": 0, "run_count": 0, "base_count": None}
            
            response.raise_for_status()
            data = response.json()
            
            if data and isinstance(data, list) and len(data) > 0:
                return {
                    "sample_count": int(data[0].get("sample_count", 0) or 0),
                    "run_count": int(data[0].get("run_count", 0) or 0),
                    "base_count": int(data[0].get("base_count", 0) or 0) if data[0].get("base_count") else None,
                }
            
            return {"sample_count": 0, "run_count": 0, "base_count": None}
            
        except Exception as e:
            logger.debug("Failed to get ENA study samples", study=study_accession, error=str(e))
            return {"sample_count": 0, "run_count": 0, "base_count": None}
    
    async def _get_study_experiments(self, study_accession: str) -> dict:
        """Get experiment/library info for a study."""
        await self._rate_limit()
        
        client = await self._get_client()
        
        try:
            response = await client.get(
                f"{ENA_PORTAL_API_URL}/search",
                params={
                    "result": "read_experiment",
                    "query": f'study_accession="{study_accession}"',
                    "fields": "library_strategy,library_source,library_selection,instrument_platform",
                    "format": "json",
                    "limit": 10,
                }
            )
            
            if response.status_code == 204:
                return {}
            
            response.raise_for_status()
            data = response.json()
            
            if data and isinstance(data, list) and len(data) > 0:
                # Get most common values
                strategies = [d.get("library_strategy") for d in data if d.get("library_strategy")]
                sources = [d.get("library_source") for d in data if d.get("library_source")]
                platforms = [d.get("instrument_platform") for d in data if d.get("instrument_platform")]
                selections = [d.get("library_selection") for d in data if d.get("library_selection")]
                
                return {
                    "library_strategy": strategies[0] if strategies else None,
                    "library_source": sources[0] if sources else None,
                    "library_selection": selections[0] if selections else None,
                    "instrument_platform": platforms[0] if platforms else None,
                }
            
            return {}
            
        except Exception as e:
            logger.debug("Failed to get ENA study experiments", study=study_accession, error=str(e))
            return {}
    
    def _determine_omics_type(
        self,
        library_strategy: Optional[str],
        library_source: Optional[str],
    ) -> OmicsType:
        """Determine omics type from library strategy and source."""
        # First check library strategy
        if library_strategy:
            omics = LIBRARY_STRATEGY_TO_OMICS.get(library_strategy)
            if omics:
                return omics
        
        # Fall back to library source
        if library_source:
            source_lower = library_source.lower()
            if "transcriptomic" in source_lower:
                return OmicsType.TRANSCRIPTOMICS
            elif "metagenomic" in source_lower:
                return OmicsType.METAGENOMICS
            elif "metatranscriptomic" in source_lower:
                return OmicsType.TRANSCRIPTOMICS
            elif "genomic" in source_lower:
                return OmicsType.GENOMICS
        
        # Default to genomics
        return OmicsType.GENOMICS
    
    def _determine_assay_types(self, library_strategy: Optional[str]) -> List[AssayType]:
        """Determine assay types from library strategy."""
        if not library_strategy:
            return []
        
        assay = LIBRARY_STRATEGY_TO_ASSAY.get(library_strategy)
        return [assay] if assay else []
    
    def _study_to_dataset(
        self,
        study: dict,
        sample_info: dict,
        experiment_info: dict,
    ) -> UnifiedDataset:
        """Convert an ENA study to UnifiedDataset."""
        accession = study.get("study_accession", "")
        secondary = study.get("secondary_study_accession")
        
        # Determine omics type and assays
        library_strategy = experiment_info.get("library_strategy")
        library_source = experiment_info.get("library_source")
        omics_type = self._determine_omics_type(library_strategy, library_source)
        assay_types = self._determine_assay_types(library_strategy)
        
        # Build extension
        extension = ENAExtension(
            study_accession=accession,
            secondary_accession=secondary,
            center_name=study.get("center_name"),
            broker_name=study.get("broker_name"),
            tax_id=int(study.get("tax_id")) if study.get("tax_id") else None,
            library_strategy=library_strategy,
            library_source=library_source,
            library_selection=experiment_info.get("library_selection"),
            instrument_platform=experiment_info.get("instrument_platform"),
            sample_count=sample_info.get("sample_count", 0),
            run_count=sample_info.get("run_count", 0),
            base_count=sample_info.get("base_count"),
        )
        
        # Build secondary accessions list
        secondary_accessions = []
        if secondary:
            secondary_accessions.append(secondary)
        
        # Organism
        organism = []
        if study.get("scientific_name"):
            organism = [study["scientific_name"]]
        
        return UnifiedDataset(
            id=f"ena:{accession}",
            accession=accession,
            source=DataSource.ENA,
            source_url=f"https://www.ebi.ac.uk/ena/browser/view/{accession}",
            secondary_accession=secondary_accessions,
            omics_type=omics_type,
            assay_types=assay_types,
            title=study.get("study_title", accession),
            description=study.get("study_description"),
            organism=organism,
            sample_count=sample_info.get("sample_count", 0),
            sample_count_display=str(sample_info.get("sample_count", 0)),
            submission_date=study.get("first_public"),
            last_update=study.get("last_updated"),
            institution=study.get("center_name"),
            curation_level=CurationLevel.COMMUNITY,  # ENA is community-submitted
            extensions={"ena": extension.model_dump()},
        )
    
    async def search(
        self,
        query: str,
        max_results: int = 50,
        organism: Optional[str] = None,
        library_strategy: Optional[str] = None,
        fetch_details: bool = True,
    ) -> List[UnifiedDataset]:
        """
        Search ENA for sequencing studies.
        
        Args:
            query: Natural language search query
            max_results: Maximum number of results
            organism: Filter by organism (uses tax_tree for hierarchical search)
            library_strategy: Filter by library strategy (RNA-Seq, WGS, etc.)
            fetch_details: Whether to fetch sample/run counts (slower but richer)
        
        Returns:
            List of UnifiedDataset results
        """
        logger.info("Searching ENA", 
                   query=query, 
                   max_results=max_results,
                   organism=organism,
                   library_strategy=library_strategy)
        
        try:
            organism = organism or self._infer_organism(query)
            # Search for studies
            studies = await self._search_studies(query, max_results, organism, library_strategy)
            
            logger.info("ENA study search complete", count=len(studies))
            
            if not studies:
                return []
            
            datasets = []
            
            # Fetch details in parallel (batched to respect rate limits)
            batch_size = 10
            for i in range(0, len(studies), batch_size):
                batch = studies[i:i + batch_size]
                
                if fetch_details:
                    # Fetch sample info and experiment info in parallel
                    tasks = []
                    for study in batch:
                        accession = study.get("study_accession", "")
                        tasks.append(self._get_study_samples(accession))
                        tasks.append(self._get_study_experiments(accession))
                    
                    results = await asyncio.gather(*tasks, return_exceptions=True)
                    
                    # Process results
                    for j, study in enumerate(batch):
                        sample_info = results[j * 2] if not isinstance(results[j * 2], Exception) else {}
                        experiment_info = results[j * 2 + 1] if not isinstance(results[j * 2 + 1], Exception) else {}
                        
                        dataset = self._study_to_dataset(study, sample_info, experiment_info)
                        datasets.append(dataset)
                else:
                    # Basic conversion without details
                    for study in batch:
                        dataset = self._study_to_dataset(study, {}, {})
                        datasets.append(dataset)
                
                if len(datasets) >= max_results:
                    break
            
            # Calculate relevance scores
            query_terms = [t.lower() for t in query.split() if len(t) > 2]
            
            for ds in datasets:
                score = 0.0
                title_lower = (ds.title or "").lower()
                desc_lower = (ds.description or "").lower()
                
                for term in query_terms:
                    if term in title_lower:
                        score += 0.3
                    if term in desc_lower:
                        score += 0.1
                    if ds.organism and any(term in o.lower() for o in ds.organism):
                        score += 0.2
                
                ds.relevance_score = min(score, 1.0)
            
            # Sort by relevance
            datasets.sort(key=lambda x: x.relevance_score, reverse=True)
            
            logger.info("ENA search complete", results=len(datasets))
            
            return datasets[:max_results]
            
        except ENASearchError:
            raise
        except Exception as e:
            logger.error("ENA search failed", error=str(e))
            raise ENASearchError(f"ENA search failed: {str(e)}")
    
    async def get_dataset(self, accession: str) -> Optional[UnifiedDataset]:
        """
        Fetch a single ENA study by accession.
        
        Args:
            accession: ENA accession (e.g., "PRJNA123456", "ERP123456")
            
        Returns:
            UnifiedDataset if found, None otherwise
        """
        try:
            # Normalize accession
            accession = accession.upper()
            
            logger.info("Fetching ENA dataset", accession=accession)
            
            await self._rate_limit()
            client = await self._get_client()
            
            # Determine query field based on accession format
            if accession.startswith("ERP"):
                query_field = "secondary_study_accession"
            elif accession.startswith("PRJ"):
                query_field = "study_accession"
            else:
                query_field = "study_accession"
            
            # Fields to retrieve
            fields = [
                "study_accession",
                "secondary_study_accession", 
                "study_title",
                "study_description",
                "scientific_name",
                "tax_id",
                "center_name",
                "broker_name",
                "first_public",
                "last_updated",
                "study_alias",
            ]
            
            # Direct query for the specific accession
            response = await client.get(
                f"{ENA_PORTAL_API_URL}/search",
                params={
                    "result": "study",
                    "query": f'{query_field}="{accession}"',
                    "fields": ",".join(fields),
                    "format": "json",
                    "limit": 1,
                }
            )
            
            if response.status_code == 204:
                # No content = no results, try secondary field
                if query_field == "study_accession":
                    response = await client.get(
                        f"{ENA_PORTAL_API_URL}/search",
                        params={
                            "result": "study",
                            "query": f'secondary_study_accession="{accession}"',
                            "fields": ",".join(fields),
                            "format": "json",
                            "limit": 1,
                        }
                    )
            
            if response.status_code == 204:
                logger.warning("ENA dataset not found", accession=accession)
                return None
            
            response.raise_for_status()
            studies = response.json()
            
            if not studies or not isinstance(studies, list) or len(studies) == 0:
                logger.warning("ENA dataset not found", accession=accession)
                return None
            
            # Get additional details
            study = studies[0]
            study_acc = study.get("study_accession", accession)
            sample_info = await self._get_study_samples(study_acc)
            experiment_info = await self._get_study_experiments(study_acc)
            
            # Convert to unified dataset
            dataset = self._study_to_dataset(study, sample_info=sample_info, experiment_info=experiment_info)
            
            logger.info("ENA dataset fetched", accession=accession)
            return dataset
            
        except Exception as e:
            logger.error("Failed to fetch ENA dataset", 
                        accession=accession, error=str(e))
            return None
    
    async def close(self):
        """Close the HTTP client."""
        if self._client and not self._client.is_closed:
            await self._client.aclose()
            self._client = None


# === SINGLETON ACCESS ===

_adapter_instance: Optional[ENAAdapter] = None


def get_ena_adapter() -> ENAAdapter:
    """Get singleton ENA adapter instance."""
    global _adapter_instance
    if _adapter_instance is None:
        _adapter_instance = ENAAdapter()
    return _adapter_instance


# === ASYNC GENERATOR FOR STREAMING ===

async def search_ena_streaming(
    query: str,
    max_results: int = 50,
    **kwargs
):
    """
    Stream ENA search results.
    Yields results as they become available.
    """
    adapter = get_ena_adapter()
    results = await adapter.search(query, max_results, **kwargs)
    
    for result in results:
        yield result
