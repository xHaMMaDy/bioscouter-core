"""
MassIVE (Mass Spectrometry Interactive Virtual Environment) Adapter for BioScouter.

Searches MassIVE for proteomics and GNPS metabolomics datasets.
Includes ProteomeXchange cross-references for deduplication with PRIDE.

API Documentation:
- Dataset JSON: https://massive.ucsd.edu/ProteoSAFe/datasets_json.jsp
- Dataset Details: https://massive.ucsd.edu/ProteoSAFe/dataset_summary.jsp
- ReDU API: https://redu.ucsd.edu/
"""

import asyncio
import re
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
    MassIVEExtension,
)

logger = structlog.get_logger(__name__)


# === CONFIGURATION ===

MASSIVE_BASE_URL = "https://massive.ucsd.edu/ProteoSAFe"
GNPS_BASE_URL = "https://gnps.ucsd.edu/ProteoSAFe"
REDU_API_URL = "https://redu.ucsd.edu"
DEFAULT_TIMEOUT = 8.0
DEFAULT_SEARCH_TIMEOUT = 10.0
DEFAULT_DATASET_RECORD_CAP = 2000
MAX_RESULTS_PER_PAGE = 100
RATE_LIMIT_DELAY = 0.5  # Be conservative with MassIVE


class MassIVESearchError(Exception):
    """Custom exception for MassIVE search errors."""
    pass


class MassIVEAdapter:
    """Adapter for searching MassIVE proteomics/metabolomics datasets."""
    
    # Source identifier
    source = DataSource.MASSIVE
    source_name = "MassIVE"
    
    def __init__(
        self,
        timeout: float = DEFAULT_TIMEOUT,
        search_timeout: float = DEFAULT_SEARCH_TIMEOUT,
        dataset_record_cap: int = DEFAULT_DATASET_RECORD_CAP,
    ):
        self.timeout = timeout
        self.search_timeout = search_timeout
        self.dataset_record_cap = dataset_record_cap
        self._client: Optional[httpx.AsyncClient] = None
        self._last_request_time: float = 0
        # Dataset list cache
        self._datasets_cache: Optional[List[dict]] = None
        self._datasets_cache_time: float = 0
        self._datasets_cache_ttl: float = 3600  # 1 hour
    
    async def _get_client(self) -> httpx.AsyncClient:
        """Get or create HTTP client."""
        if self._client is None or self._client.is_closed:
            self._client = httpx.AsyncClient(
                timeout=httpx.Timeout(self.timeout),
                limits=build_default_limits(),
                headers={
                    "Accept": "application/json",
                    "User-Agent": "BioScouter/1.0 (Multi-omics Search Platform)",
                },
                follow_redirects=True,
            )
        return self._client
    
    async def _rate_limit(self):
        """Apply rate limiting."""
        now = asyncio.get_event_loop().time()
        elapsed = now - self._last_request_time
        if elapsed < RATE_LIMIT_DELAY:
            await asyncio.sleep(RATE_LIMIT_DELAY - elapsed)
        self._last_request_time = asyncio.get_event_loop().time()
    
    async def _fetch_all_datasets(self) -> List[dict]:
        """
        Fetch all datasets from MassIVE.
        Cached for 1 hour.
        """
        now = asyncio.get_event_loop().time()
        
        if self._datasets_cache and (now - self._datasets_cache_time) < self._datasets_cache_ttl:
            return self._datasets_cache
        
        logger.info("Fetching MassIVE dataset list")
        
        await self._rate_limit()
        client = await self._get_client()
        
        try:
            # Fetch JSON dataset list
            response = await client.get(
                f"{MASSIVE_BASE_URL}/datasets_json.jsp",
                params={"pageSize": self.dataset_record_cap, "offset": 0},
            )
            response.raise_for_status()
            
            # Handle potential encoding issues in MassIVE response
            try:
                data = response.json()
            except Exception:
                # Try decoding with error handling
                content = response.content.decode('utf-8', errors='replace')
                import json
                data = json.loads(content)
            
            # Extract datasets from response
            if isinstance(data, dict) and "datasets" in data:
                datasets = data["datasets"]
            elif isinstance(data, list):
                datasets = data
            else:
                datasets = []
            
            self._datasets_cache = datasets
            self._datasets_cache_time = now
            
            logger.info("MassIVE datasets cached", count=len(datasets))
            
            return datasets
            
        except Exception as e:
            logger.error("Failed to fetch MassIVE datasets", error=str(e))
            raise MassIVESearchError(f"Failed to fetch MassIVE datasets: {str(e)}")
    
    async def _fetch_dataset_details(self, task_id: str) -> Optional[dict]:
        """Fetch detailed information for a specific dataset."""
        await self._rate_limit()
        client = await self._get_client()
        
        try:
            # Try to get dataset info from result.xml or JSON endpoint
            response = await client.get(
                f"{MASSIVE_BASE_URL}/result_json.jsp",
                params={"task": task_id, "view": "view_all_datasets"},
            )
            
            if response.status_code == 200:
                try:
                    data = response.json()
                    if isinstance(data, dict):
                        return data
                except:
                    pass
            
            return None
            
        except Exception as e:
            logger.debug("Failed to fetch MassIVE dataset details", 
                        task_id=task_id, error=str(e))
            return None
    
    def _extract_msv_id(self, dataset: dict) -> Optional[str]:
        """Extract MSV accession from dataset."""
        # Check various fields where MSV ID might be stored
        for field in ["dataset", "massive_id", "msv_id", "accession"]:
            value = dataset.get(field, "")
            if isinstance(value, str) and value.startswith("MSV"):
                return value
        
        # Try to extract from task ID
        task = dataset.get("task", "")
        if task:
            # MSV ID is sometimes encoded in the task
            msv_match = re.search(r"MSV\d+", str(task))
            if msv_match:
                return msv_match.group()
        
        return None
    
    def _extract_px_accession(self, dataset: dict) -> Optional[str]:
        """Extract ProteomeXchange accession from dataset."""
        # Check description and title for PXD IDs
        text = " ".join([
            str(dataset.get("title", "")),
            str(dataset.get("description", "")),
            str(dataset.get("comment", "")),
        ])
        
        px_match = re.search(r"PXD\d{6,}", text)
        if px_match:
            return px_match.group()
        
        # Check specific field
        px_id = dataset.get("px_accession") or dataset.get("proteomexchange")
        if px_id and px_id.startswith("PXD"):
            return px_id
        
        return None
    
    def _is_gnps_dataset(self, dataset: dict) -> bool:
        """Check if dataset is from GNPS (metabolomics)."""
        title = str(dataset.get("title", "")).lower()
        desc = str(dataset.get("description", "")).lower()
        
        gnps_keywords = ["gnps", "molecular networking", "fbmn", "metabolomics", 
                        "natural product", "mass spec network"]
        
        return any(kw in title or kw in desc for kw in gnps_keywords)
    
    def _is_reanalysis(self, dataset: dict) -> bool:
        """Check if dataset is a reanalysis."""
        title = str(dataset.get("title", "")).lower()
        desc = str(dataset.get("description", "")).lower()
        
        reanalysis_keywords = ["reanalysis", "re-analysis", "benchmark", "test data",
                              "tutorial", "demo", "example"]
        
        return any(kw in title or kw in desc for kw in reanalysis_keywords)
    
    def _extract_species(self, dataset: dict) -> List[str]:
        """Extract species from dataset metadata."""
        species = []
        
        # Check species field
        if dataset.get("species"):
            sp = dataset["species"]
            if isinstance(sp, list):
                species.extend(sp)
            elif isinstance(sp, str):
                species.append(sp)
        
        # Check title/description for common species
        text = " ".join([
            str(dataset.get("title", "")),
            str(dataset.get("description", "")),
        ]).lower()
        
        species_map = {
            "human": "Homo sapiens",
            "homo sapiens": "Homo sapiens",
            "mouse": "Mus musculus",
            "mus musculus": "Mus musculus",
            "rat": "Rattus norvegicus",
            "rattus": "Rattus norvegicus",
            "yeast": "Saccharomyces cerevisiae",
            "saccharomyces": "Saccharomyces cerevisiae",
            "e. coli": "Escherichia coli",
            "escherichia": "Escherichia coli",
            "arabidopsis": "Arabidopsis thaliana",
        }
        
        for keyword, scientific_name in species_map.items():
            if keyword in text and scientific_name not in species:
                species.append(scientific_name)
                break  # Take first match
        
        return species if species else ["Unknown"]
    
    def _extract_instruments(self, dataset: dict) -> List[str]:
        """Extract MS instruments from dataset."""
        instruments = []
        
        # Check instrument field
        if dataset.get("instrument"):
            inst = dataset["instrument"]
            if isinstance(inst, list):
                instruments.extend(inst)
            elif isinstance(inst, str):
                instruments.append(inst)
        
        # Check description for common instruments
        text = str(dataset.get("description", "")).lower()
        
        instrument_keywords = [
            "orbitrap", "q-tof", "qtof", "triple quad", "qq", "tof",
            "thermo", "bruker", "agilent", "waters", "sciex", "shimadzu"
        ]
        
        for kw in instrument_keywords:
            if kw in text and kw.title() not in [i.lower() for i in instruments]:
                instruments.append(kw.title())
        
        return instruments
    
    def _extract_modifications(self, dataset: dict) -> List[str]:
        """Extract PTMs from dataset."""
        mods = []
        
        text = " ".join([
            str(dataset.get("title", "")),
            str(dataset.get("description", "")),
        ]).lower()
        
        ptm_keywords = {
            "phospho": "Phosphorylation",
            "acetyl": "Acetylation",
            "methyl": "Methylation",
            "ubiquitin": "Ubiquitination",
            "glyco": "Glycosylation",
            "sumo": "SUMOylation",
            "nitro": "Nitrosylation",
        }
        
        for keyword, name in ptm_keywords.items():
            if keyword in text:
                mods.append(name)
        
        return mods
    
    def _determine_omics_type(self, dataset: dict) -> OmicsType:
        """Determine omics type from dataset."""
        if self._is_gnps_dataset(dataset):
            return OmicsType.METABOLOMICS
        return OmicsType.PROTEOMICS
    
    def _determine_assay_types(self, dataset: dict) -> List[AssayType]:
        """Determine assay types from dataset."""
        text = str(dataset.get("description", "")).lower()
        assays = []
        
        if "tmt" in text:
            assays.append(AssayType.TMT)
        if "itraq" in text:
            assays.append(AssayType.ITRAQ)
        if "silac" in text:
            assays.append(AssayType.SILAC)
        if "label-free" in text or "label free" in text:
            assays.append(AssayType.LABEL_FREE)
        if "dia" in text or "swath" in text:
            assays.append(AssayType.DIA)
        if "dda" in text:
            assays.append(AssayType.DDA)
        if "lc-ms" in text or "lcms" in text:
            assays.append(AssayType.LC_MS)
        
        return assays if assays else [AssayType.LABEL_FREE]  # Default
    
    def _dataset_matches_query(self, dataset: dict, query_terms: List[str]) -> bool:
        """Check if dataset matches search terms."""
        searchable = " ".join([
            str(dataset.get("title", "")),
            str(dataset.get("description", "")),
            str(dataset.get("species", "")),
            str(dataset.get("comment", "")),
        ]).lower()
        
        # Check if ANY term matches (more lenient)
        matches = sum(1 for term in query_terms if term in searchable)
        return matches >= max(1, len(query_terms) // 2)
    
    def _dataset_to_unified(self, dataset: dict) -> UnifiedDataset:
        """Convert MassIVE dataset to UnifiedDataset."""
        msv_id = self._extract_msv_id(dataset) or f"MSV{dataset.get('task', 'unknown')}"
        px_accession = self._extract_px_accession(dataset)
        is_gnps = self._is_gnps_dataset(dataset)
        is_reanalysis = self._is_reanalysis(dataset)
        
        title = dataset.get("title", msv_id)
        description = dataset.get("description")
        
        # Get metadata
        species = self._extract_species(dataset)
        instruments = self._extract_instruments(dataset)
        modifications = self._extract_modifications(dataset)
        omics_type = self._determine_omics_type(dataset)
        assay_types = self._determine_assay_types(dataset)
        
        # File info
        file_count = int(dataset.get("files_count", 0) or 0)
        total_size = None
        if dataset.get("size"):
            try:
                total_size = int(dataset["size"])
            except:
                pass
        
        # Check for GNPS visualization
        gnps_viz_url = None
        molecular_networking = False
        if is_gnps:
            gnps_task = dataset.get("gnps_task") or dataset.get("task")
            if gnps_task:
                gnps_viz_url = f"https://gnps.ucsd.edu/ProteoSAFe/result.jsp?task={gnps_task}&view=network_displayer"
                molecular_networking = True
        
        # Build extension
        extension = MassIVEExtension(
            msv_id=msv_id,
            px_accession=px_accession,
            gnps_task_id=dataset.get("gnps_task"),
            submitter=dataset.get("user"),
            file_count=file_count,
            total_size_bytes=total_size,
            instruments=instruments,
            species=species,
            modifications=modifications,
            is_gnps_dataset=is_gnps,
            is_reanalysis=is_reanalysis,
            has_redu_metadata=False,  # Would need separate check
            molecular_networking_available=molecular_networking,
            gnps_visualization_url=gnps_viz_url,
        )
        
        # Secondary accessions
        secondary = []
        if px_accession:
            secondary.append(px_accession)
        
        # Source URL
        source_url = f"https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?accession={msv_id}"
        if is_gnps:
            source_url = f"https://gnps.ucsd.edu/ProteoSAFe/dataset.jsp?accession={msv_id}"
        
        # Submission date
        submission_date = dataset.get("create_time") or dataset.get("created")
        
        return UnifiedDataset(
            id=f"massive:{msv_id}",
            accession=msv_id,
            source=DataSource.MASSIVE,
            source_url=source_url,
            secondary_accession=secondary,
            omics_type=omics_type,
            assay_types=assay_types,
            title=title,
            description=description,
            organism=species,
            sample_count=file_count,  # Use file count as proxy
            sample_count_display=f"{file_count} files",
            submission_date=submission_date,
            contributors=[dataset.get("user")] if dataset.get("user") else [],
            curation_level=CurationLevel.COMMUNITY,
            extensions={"massive": extension.model_dump()},
        )
    
    async def search(
        self,
        query: str,
        max_results: int = 50,
        organism: Optional[str] = None,  # Accepted for interface compatibility
        include_gnps: bool = True,
        include_reanalysis: bool = True,
        fetch_details: bool = True,
    ) -> List[UnifiedDataset]:
        """
        Search MassIVE for proteomics/metabolomics datasets.
        
        Args:
            query: Natural language search query
            max_results: Maximum number of results
            organism: Organism filter (used in text matching)
            include_gnps: Include GNPS metabolomics datasets
            include_reanalysis: Include reanalysis/benchmark datasets
            fetch_details: Fetch additional details for each dataset
        
        Returns:
            List of UnifiedDataset results
        """
        logger.info("Searching MassIVE", 
                   query=query, 
                   max_results=max_results,
                   include_gnps=include_gnps,
                   include_reanalysis=include_reanalysis)
        
        # Extract query terms
        query_terms = [t.lower() for t in query.split() if len(t) > 2]
        
        try:
            # Fetch all datasets (cached)
            all_datasets = await asyncio.wait_for(
                self._fetch_all_datasets(),
                timeout=self.search_timeout,
            )
            
            logger.info("MassIVE filtering datasets", total=len(all_datasets))
            
            # Filter datasets
            matching = []
            for dataset in all_datasets:
                # Apply filters
                is_gnps = self._is_gnps_dataset(dataset)
                is_reanalysis = self._is_reanalysis(dataset)
                
                if is_gnps and not include_gnps:
                    continue
                if is_reanalysis and not include_reanalysis:
                    continue
                
                # Match query
                if not self._dataset_matches_query(dataset, query_terms):
                    continue
                
                matching.append(dataset)
                
                if len(matching) >= max_results * 2:  # Get extra for scoring
                    break
            
            logger.info("MassIVE pre-filtered", count=len(matching))
            
            # Convert to unified datasets
            datasets = []
            for dataset in matching:
                try:
                    unified = self._dataset_to_unified(dataset)
                    datasets.append(unified)
                except Exception as e:
                    logger.debug("Failed to convert MassIVE dataset", error=str(e))
                    continue
            
            # Fetch additional details if requested
            if fetch_details and datasets:
                # Batch fetch details
                batch_size = 5
                for i in range(0, min(len(datasets), max_results), batch_size):
                    batch = datasets[i:i + batch_size]
                    tasks = []
                    for ds in batch:
                        msv_id = ds.accession
                        # Find original dataset
                        orig = next((d for d in matching if self._extract_msv_id(d) == msv_id), None)
                        if orig and orig.get("task"):
                            tasks.append(self._fetch_dataset_details(orig["task"]))
                        else:
                            tasks.append(asyncio.sleep(0))  # Placeholder
                    
                    # Fetch in parallel
                    await asyncio.wait_for(
                        asyncio.gather(*tasks, return_exceptions=True),
                        timeout=min(2.0, self.search_timeout),
                    )
            
            # Calculate relevance scores
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
            
            logger.info("MassIVE search complete", results=len(datasets))
            
            return datasets[:max_results]
            
        except asyncio.TimeoutError:
            logger.warning(
                "MassIVE search timed out",
                timeout_seconds=self.search_timeout,
                query=query[:100],
            )
            return []
        except MassIVESearchError:
            raise
        except Exception as e:
            logger.error("MassIVE search failed", error=str(e))
            return []  # Graceful failure
    
    async def get_dataset(self, accession: str) -> Optional[UnifiedDataset]:
        """
        Fetch a single MassIVE dataset by accession.
        
        Args:
            accession: MSV accession (e.g., "MSV000099441")
            
        Returns:
            UnifiedDataset if found, None otherwise
        """
        try:
            # Normalize accession
            if not accession.upper().startswith("MSV"):
                accession = f"MSV{accession}"
            accession = accession.upper()
            
            logger.info("Fetching MassIVE dataset", accession=accession)
            
            # Fetch all datasets (cached) and find the one we want
            all_datasets = await self._fetch_all_datasets()
            
            # Find matching dataset
            matching_dataset = None
            for dataset in all_datasets:
                msv_id = self._extract_msv_id(dataset)
                if msv_id and msv_id.upper() == accession:
                    matching_dataset = dataset
                    break
            
            if not matching_dataset:
                logger.warning("MassIVE dataset not found", accession=accession)
                return None
            
            # Convert to unified dataset
            unified = self._dataset_to_unified(matching_dataset)
            
            logger.info("MassIVE dataset fetched", accession=accession)
            return unified
            
        except Exception as e:
            logger.error("Failed to fetch MassIVE dataset", 
                        accession=accession, error=str(e))
            return None
    
    async def close(self):
        """Close the HTTP client."""
        if self._client and not self._client.is_closed:
            await self._client.aclose()
            self._client = None


# === SINGLETON ACCESS ===

_adapter_instance: Optional[MassIVEAdapter] = None


def get_massive_adapter() -> MassIVEAdapter:
    """Get singleton MassIVE adapter instance."""
    global _adapter_instance
    if _adapter_instance is None:
        _adapter_instance = MassIVEAdapter()
    return _adapter_instance


# === ASYNC GENERATOR FOR STREAMING ===

async def search_massive_streaming(
    query: str,
    max_results: int = 50,
    **kwargs
):
    """
    Stream MassIVE search results.
    Yields results as they become available.
    """
    adapter = get_massive_adapter()
    results = await adapter.search(query, max_results, **kwargs)
    
    for result in results:
        yield result
