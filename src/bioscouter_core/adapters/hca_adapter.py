"""
Human Cell Atlas (HCA) Data Portal Adapter
Searches HCA Data Portal for single-cell datasets via the Azul API.
API Docs: https://service.azul.data.humancellatlas.org/
"""

import asyncio
import time
from typing import List, Optional, Any, Set, Dict
from datetime import datetime, timezone

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
    CurationLevel,
    UnifiedDataset,
    HCAExtension,
    DownloadLink,
    DownloadFileType,
)

logger = structlog.get_logger(__name__)

# HCA Azul API base URL
HCA_API_BASE = "https://service.azul.data.humancellatlas.org"

# Fallback used only when catalog discovery is unavailable.
DEFAULT_CATALOG = "dcp59"

# HCA-specific stopwords
HCA_STOPWORDS: Set[str] = frozenset({
    'single', 'cell', 'cells', 'scrna', 'scrnaseq', 'rna', 'seq',
    'sequencing', 'atlas', 'human', 'hca'
})

# Cache duration for projects list (1 hour as per plan)
PROJECTS_CACHE_TTL = 3600  # 1 hour


class HCAError(Exception):
    """HCA API error."""
    pass


class HCAAdapter(BaseSourceAdapter):
    """
    Adapter for Human Cell Atlas Data Portal.
    Searches for single-cell datasets via the Azul API.
    Expert-curated datasets with standardized annotations.
    """
    
    def __init__(self, timeout: float = 30.0, catalog: Optional[str] = None):
        super().__init__()
        self.timeout = timeout
        self.catalog = catalog or DEFAULT_CATALOG
        self._catalog_explicit = catalog is not None
        self._client: Optional[httpx.AsyncClient] = None
        
        # Cache for projects list
        self._projects_cache: Optional[List[dict]] = None
        self._projects_cache_time: float = 0
        
        # Catalog version tracking for admin notifications
        self._catalog_version: Optional[str] = self.catalog if self._catalog_explicit else None
        self._catalog_version_check_time: float = 0
    
    @property
    def source(self) -> DataSource:
        return DataSource.HCA
    
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
                    "User-Agent": "BioScouter/1.0",
                }
            )
        return self._client
    
    async def close(self):
        """Close HTTP client."""
        if self._client and not self._client.is_closed:
            await self._client.aclose()
            self._client = None
    
    async def _fetch_with_retry(self, url: str, params: Optional[Dict] = None) -> httpx.Response:
        """Fetch URL with retry logic. Delegates to the shared base helper."""
        client = await self._get_client()
        return await super()._fetch_with_retry(client, url, params=params)

    async def _resolve_catalog(self) -> str:
        """Discover Azul's current default catalog once per adapter instance."""
        if self._catalog_version:
            return self._catalog_version

        try:
            response = await self._fetch_with_retry(f"{HCA_API_BASE}/index/catalogs")
            response.raise_for_status()
            discovered = response.json().get("default_catalog")
            if discovered:
                self.catalog = discovered
                self._catalog_version = discovered
                self._catalog_version_check_time = time.time()
                logger.info("Discovered HCA default catalog", catalog=discovered)
                return discovered
        except Exception as e:
            logger.warning(
                "HCA catalog discovery failed; using fallback",
                catalog=self.catalog,
                error=str(e),
            )

        self._catalog_version = self.catalog
        return self.catalog
    
    async def _get_projects_cached(self, max_results: int = 500) -> List[dict]:
        """
        Get projects with caching to avoid repeated API calls.
        Uses on-demand pagination with max_results limit.
        """
        current_time = time.time()
        
        # Check if we have a valid cached copy
        if (self._projects_cache is not None and 
            current_time - self._projects_cache_time < PROJECTS_CACHE_TTL):
            logger.debug("Using cached HCA projects", count=len(self._projects_cache))
            return self._projects_cache
        
        logger.info("Fetching HCA projects", catalog=self.catalog)
        catalog = await self._resolve_catalog()
        
        # Fetch projects with pagination
        # Note: HCA Azul API may reject large page sizes (>25) with 400
        all_projects = []
        url = f"{HCA_API_BASE}/index/projects"
        params = {
            "catalog": catalog,
            "size": min(25, max_results),  # Page size (API limit ~25)
        }
        
        while len(all_projects) < max_results:
            response = await self._fetch_with_retry(url, params)
            response.raise_for_status()
            
            data = response.json()
            hits = data.get("hits", [])
            
            if not hits:
                break
            
            all_projects.extend(hits)
            
            # Check for next page
            pagination = data.get("pagination", {})
            search_after = pagination.get("search_after")
            
            if not search_after or len(all_projects) >= max_results:
                break
            
            # Use search_after for pagination
            params["search_after"] = search_after
        
        # Cache the result
        self._projects_cache = all_projects[:max_results]
        self._projects_cache_time = current_time
        
        logger.info("HCA projects cached", count=len(self._projects_cache))
        return self._projects_cache
    
    def _parse_organism(self, organism: Optional[str]) -> Optional[str]:
        """Convert organism name to HCA format."""
        if not organism:
            return None
        
        organism_lower = organism.lower()
        organism_mapping = {
            "human": "Homo sapiens",
            "homo sapiens": "Homo sapiens",
            "mouse": "Mus musculus",
            "mus musculus": "Mus musculus",
        }
        
        return organism_mapping.get(organism_lower, organism)
    
    def _detect_assay_type(self, methods: List[str]) -> List[AssayType]:
        """Detect assay types from HCA library construction methods."""
        assay_types = []
        
        for method in methods:
            method_lower = method.lower() if method else ""
            
            if "10x" in method_lower or "chromium" in method_lower:
                assay_types.append(AssayType.SCRNA_SEQ)
            elif "smart-seq" in method_lower:
                assay_types.append(AssayType.SMART_SEQ)
            elif "cite-seq" in method_lower:
                assay_types.append(AssayType.CITE_SEQ)
            elif "drop-seq" in method_lower:
                assay_types.append(AssayType.SCRNA_SEQ)
            elif "multiome" in method_lower:
                assay_types.append(AssayType.MULTIOME)
            elif "atac" in method_lower:
                assay_types.append(AssayType.ATAC_SEQ)
            elif "sn" in method_lower and "rna" in method_lower:
                assay_types.append(AssayType.SNRNA_SEQ)
            elif "spatial" in method_lower:
                assay_types.append(AssayType.SPATIAL_TRANSCRIPTOMICS)
        
        # Default to scRNA-seq if nothing detected
        if not assay_types:
            assay_types.append(AssayType.SCRNA_SEQ)
        
        return list(set(assay_types))  # Remove duplicates
    
    def _extract_secondary_accessions(self, project: dict) -> List[str]:
        """Extract GEO/SRA/ENA accessions from project for cross-source deduplication."""
        accessions = []
        
        # Check various fields where external accessions may be stored
        protocols = project.get("protocols", [])
        for protocol in protocols:
            if isinstance(protocol, dict):
                # Check for linked accessions in protocol metadata
                supplementary_links = protocol.get("supplementary_links", [])
                for link in supplementary_links:
                    if isinstance(link, str):
                        if "ncbi.nlm.nih.gov/geo" in link:
                            # Extract GSE ID from URL
                            import re
                            match = re.search(r'GSE\d+', link)
                            if match:
                                accessions.append(match.group())
                        elif "ncbi.nlm.nih.gov/sra" in link or "SRP" in link or "SRR" in link:
                            match = re.search(r'(SRP|SRR|SRX)\d+', link)
                            if match:
                                accessions.append(match.group())
        
        # Check project supplementary links
        sources = project.get("sources", [])
        for source in sources:
            if isinstance(source, dict):
                source_id = source.get("sourceId", "")
                if source_id.startswith("GSE"):
                    accessions.append(source_id)
                elif source_id.startswith(("SRP", "SRR", "SRX")):
                    accessions.append(source_id)
        
        # Also check contributorMatrices and supplementaryLinks
        contributor_matrices = project.get("contributorMatrices", {})
        if isinstance(contributor_matrices, dict):
            for matrix_info in contributor_matrices.values():
                if isinstance(matrix_info, dict):
                    for link in matrix_info.get("links", []):
                        if isinstance(link, str):
                            if "geo" in link.lower():
                                import re
                                match = re.search(r'GSE\d+', link)
                                if match:
                                    accessions.append(match.group())
        
        return list(set(accessions))  # Remove duplicates
    
    def _get_cell_count(self, project: dict) -> int:
        """
        Get cell count with fallback chain:
        effectiveCellCount → cellCount → matrixCellCount → 0
        """
        cell_suspensions = project.get("cellSuspensions", [])
        
        total_cells = 0
        for cs in cell_suspensions:
            if isinstance(cs, dict):
                # Try fallback chain
                count = (
                    cs.get("effectiveCellCount") or 
                    cs.get("cellCount") or 
                    cs.get("totalCells") or 
                    0
                )
                if isinstance(count, (int, float)):
                    total_cells += int(count)
        
        # Also check project-level cell count
        if total_cells == 0:
            project_cell_count = project.get("cellCount") or project.get("totalCells") or 0
            if isinstance(project_cell_count, (int, float)):
                total_cells = int(project_cell_count)
        
        return total_cells
    
    async def search(
        self,
        query: str,
        max_results: int = 50,
        organism: Optional[str] = None,
        **kwargs
    ) -> List[UnifiedDataset]:
        """
        Search HCA for single-cell datasets.
        
        Args:
            query: Search query
            max_results: Maximum results to return
            organism: Optional organism filter
            
        Returns:
            List of UnifiedDataset objects
        """
        try:
            logger.info("Searching HCA", query=query, max_results=max_results)
            
            # Get projects with caching
            projects = await self._get_projects_cached(max_results=max_results * 3)
            
            logger.info("HCA response", project_count=len(projects))
            
            # Parse organism filter
            target_organism = self._parse_organism(organism)
            
            # Extract search terms
            query_lower = query.lower()
            search_terms = [
                term for term in query_lower.split()
                if term not in COMMON_STOPWORDS and term not in HCA_STOPWORDS and len(term) > 2
            ]
            
            # Filter and convert projects
            unified = []
            for project_hit in projects:
                try:
                    project = project_hit.get("projects", [{}])[0] if project_hit.get("projects") else {}
                    if not project:
                        continue
                    
                    ud = self._to_unified(
                        project_hit,
                        query=query,
                        search_terms=search_terms,
                        target_organism=target_organism
                    )
                    if ud:
                        unified.append(ud)
                except Exception as e:
                    logger.warning("Failed to parse HCA project", error=str(e))
                    continue
            
            # Sort by relevance score
            unified.sort(key=lambda x: x.relevance_score, reverse=True)
            
            # Return top results
            unified = unified[:max_results]
            
            logger.info("HCA search complete", results=len(unified))
            return unified
            
        except httpx.HTTPStatusError as e:
            logger.error("HCA API error", status=e.response.status_code, error=str(e))
            raise HCAError(f"HCA API error: {e.response.status_code}") from e
        except httpx.RequestError as e:
            logger.error("HCA request failed", error=str(e))
            raise HCAError(f"HCA request failed: {str(e)}") from e
    
    async def get_dataset(self, accession: str) -> Optional[UnifiedDataset]:
        """Fetch a single HCA project by ID."""
        try:
            url = f"{HCA_API_BASE}/index/projects/{accession}"
            params = {"catalog": await self._resolve_catalog()}
            
            response = await self._fetch_with_retry(url, params)
            
            if response.status_code == 404:
                return None
            
            response.raise_for_status()
            project_hit = response.json()
            
            return self._to_unified(project_hit)
            
        except httpx.HTTPStatusError as e:
            if e.response.status_code == 404:
                return None
            logger.error("HCA API error", accession=accession, error=str(e))
            return None
        except Exception as e:
            logger.error("Failed to fetch HCA project", accession=accession, error=str(e))
            return None
    
    def _to_unified(
        self, 
        project_hit: dict, 
        query: Optional[str] = None,
        search_terms: Optional[List[str]] = None,
        target_organism: Optional[str] = None
    ) -> Optional[UnifiedDataset]:
        """Convert an HCA project to UnifiedDataset."""
        try:
            # Extract project data from hit
            projects = project_hit.get("projects", [])
            if not projects:
                return None
            
            project = projects[0]
            project_id = project.get("projectId") or ""
            if not project_id:
                return None
            
            # Basic info
            title = project.get("projectTitle") or project.get("projectShortname") or ""
            description = project.get("projectDescription") or ""
            short_name = project.get("projectShortname") or ""
            
            # Extract organisms from specimens/donors
            organisms = set()
            specimens = project_hit.get("specimens", [])
            for spec in specimens:
                if isinstance(spec, dict):
                    genus_species = spec.get("genusSpecies", [])
                    if isinstance(genus_species, list):
                        organisms.update(genus_species)
                    elif isinstance(genus_species, str):
                        organisms.add(genus_species)
            
            # Also check donors
            donors = project_hit.get("donorOrganisms", [])
            for donor in donors:
                if isinstance(donor, dict):
                    genus_species = donor.get("genusSpecies", [])
                    if isinstance(genus_species, list):
                        organisms.update(genus_species)
                    elif isinstance(genus_species, str):
                        organisms.add(genus_species)
            
            # Filter empty organisms
            organisms = {o for o in organisms if o}
            
            # Organism filter
            if target_organism and organisms:
                if not any(target_organism.lower() in o.lower() for o in organisms):
                    return None
            
            # Extract tissues/organs
            tissues = set()
            organs = set()
            organ_parts = set()
            
            for spec in specimens:
                if isinstance(spec, dict):
                    organ_list = spec.get("organ", [])
                    if isinstance(organ_list, list):
                        organs.update(organ_list)
                    elif isinstance(organ_list, str):
                        organs.add(organ_list)
                    
                    organ_part_list = spec.get("organPart", [])
                    if isinstance(organ_part_list, list):
                        organ_parts.update(organ_part_list)
                    elif isinstance(organ_part_list, str):
                        organ_parts.add(organ_part_list)
            
            # Filter out None values
            organs = {o for o in organs if o}
            organ_parts = {op for op in organ_parts if op}
            
            tissues = organs | organ_parts
            tissues = {t for t in tissues if t}
            
            # Extract diseases
            diseases = set()
            for spec in specimens:
                if isinstance(spec, dict):
                    disease_list = spec.get("disease", [])
                    if isinstance(disease_list, list):
                        diseases.update(disease_list)
                    elif isinstance(disease_list, str):
                        diseases.add(disease_list)
            
            diseases = {d for d in diseases if d and d.lower() != "normal"}
            
            # Extract library construction methods for assay detection
            library_methods = []
            cell_suspensions = project_hit.get("cellSuspensions", [])
            protocols = project_hit.get("protocols", [])
            
            for protocol in protocols:
                if isinstance(protocol, dict):
                    library_construction = protocol.get("libraryConstructionApproach", [])
                    if isinstance(library_construction, list):
                        library_methods.extend(library_construction)
                    elif isinstance(library_construction, str):
                        library_methods.append(library_construction)
            
            # Also check cellSuspensions
            for cs in cell_suspensions:
                if isinstance(cs, dict):
                    selected_cell_type = cs.get("selectedCellType", [])
                    # Cell suspension doesn't typically have library method, check protocols
            
            assay_types = self._detect_assay_type(library_methods)
            
            # Get cell count with fallback chain
            cell_count = self._get_cell_count(project_hit)
            
            # Sample count (donors)
            donor_count = 0
            for donor in donors:
                if isinstance(donor, dict):
                    donor_count += donor.get("donorCount", 1)
            
            sample_count = donor_count if donor_count > 0 else len(specimens)
            
            # Extract secondary accessions for cross-source deduplication
            secondary_accessions = self._extract_secondary_accessions(project_hit)
            
            # Also extract from project metadata
            geo_accessions = []
            sra_accessions = []
            ena_accessions = []
            
            accessions_data = project.get("accessions", [])
            for acc in accessions_data:
                if isinstance(acc, dict):
                    acc_value = acc.get("accession", "")
                    if acc_value.startswith("GSE"):
                        geo_accessions.append(acc_value)
                        if acc_value not in secondary_accessions:
                            secondary_accessions.append(acc_value)
                    elif acc_value.startswith(("SRP", "SRR", "SRX")):
                        sra_accessions.append(acc_value)
                        if acc_value not in secondary_accessions:
                            secondary_accessions.append(acc_value)
                    elif acc_value.startswith(("ERP", "ERR", "ERX")):
                        ena_accessions.append(acc_value)
                        if acc_value not in secondary_accessions:
                            secondary_accessions.append(acc_value)
            
            # Extract publications
            pubmed_ids = []
            doi = None
            publications = project.get("publications", [])
            for pub in publications:
                if isinstance(pub, dict):
                    pmid = pub.get("publicationId") or pub.get("pmid")
                    if pmid:
                        pubmed_ids.append(str(pmid))
                    pub_doi = pub.get("doi")
                    if pub_doi and not doi:
                        doi = pub_doi
            
            # Extract contributors
            contributors = []
            contrib_data = project.get("contributors", [])
            for contrib in contrib_data:
                if isinstance(contrib, dict):
                    name = contrib.get("contactName") or contrib.get("name")
                    if name:
                        contributors.append(name)
            
            # Extract institution
            institution = None
            if contrib_data:
                first_contrib = contrib_data[0] if isinstance(contrib_data[0], dict) else {}
                institution = first_contrib.get("institution") or first_contrib.get("laboratory")
            
            # Dates
            submission_date = project.get("submissionDate")
            last_update = project.get("updateDate")
            
            # Calculate relevance score
            relevance_score = self._calculate_relevance(
                title=title,
                description=description,
                organisms=list(organisms),
                tissues=list(tissues),
                diseases=list(diseases),
                query=query,
                search_terms=search_terms or []
            )
            
            # Build download links - point to HCA Data Browser
            download_links = []
            browser_url = f"https://data.humancellatlas.org/explore/projects/{project_id}"
            download_links.append(DownloadLink(
                url=browser_url,
                file_type=DownloadFileType.OTHER,
                file_name=None,
                file_size_bytes=None,
                protocol="https",
                description="View and download in HCA Data Browser"
            ))
            
            # Also add direct matrix download link if available
            matrices = project_hit.get("matrices", {})
            if matrices:
                for matrix_key, matrix_info in matrices.items():
                    if isinstance(matrix_info, dict):
                        matrix_url = matrix_info.get("url")
                        if matrix_url:
                            download_links.append(DownloadLink(
                                url=matrix_url,
                                file_type=DownloadFileType.MATRIX,
                                file_name=matrix_info.get("name"),
                                file_size_bytes=matrix_info.get("size"),
                                protocol="https",
                                description="Expression matrix"
                            ))
            
            # Development stages
            development_stages = []
            for donor in donors:
                if isinstance(donor, dict):
                    dev_stage = donor.get("developmentStage", [])
                    if isinstance(dev_stage, list):
                        development_stages.extend(dev_stage)
                    elif isinstance(dev_stage, str):
                        development_stages.append(dev_stage)
            
            # Nucleic acid source
            nucleic_acid_sources = []
            for cs in cell_suspensions:
                if isinstance(cs, dict):
                    nas = cs.get("nucleicAcidSource", [])
                    if isinstance(nas, list):
                        nucleic_acid_sources.extend(nas)
                    elif isinstance(nas, str):
                        nucleic_acid_sources.append(nas)
            
            # File formats
            file_formats = []
            files = project_hit.get("files", [])
            for file_info in files:
                if isinstance(file_info, dict):
                    file_format = file_info.get("format") or file_info.get("fileFormat")
                    if file_format and file_format not in file_formats:
                        file_formats.append(file_format)
            
            # Match reasons
            match_reasons = []
            if search_terms:
                title_lower = title.lower()
                desc_lower = description.lower()
                for term in search_terms:
                    if term in title_lower:
                        match_reasons.append(f"Title contains '{term}'")
                    elif term in desc_lower:
                        match_reasons.append(f"Description contains '{term}'")
                
                for org in organisms:
                    if any(term in org.lower() for term in search_terms):
                        match_reasons.append(f"Organism: {org}")
                
                for tissue in tissues:
                    if any(term in tissue.lower() for term in search_terms):
                        match_reasons.append(f"Tissue: {tissue}")
                
                for disease in diseases:
                    if any(term in disease.lower() for term in search_terms):
                        match_reasons.append(f"Disease: {disease}")
            
            if not match_reasons:
                match_reasons.append("General match")
            
            # Build extension
            extension = HCAExtension(
                project_id=project_id,
                project_short_name=short_name,
                cell_count=cell_count,
                effective_cell_count=cell_count,  # Same for now
                organ=list(organs),
                organ_part=list(organ_parts),
                development_stage=list(set(development_stages)),
                library_construction_method=list(set(library_methods)),
                nucleic_acid_source=list(set(nucleic_acid_sources)),
                file_formats=list(set(file_formats)),
                protocols=[],  # Could extract protocol types
                donor_count=donor_count,
                specimen_count=len(specimens),
                geo_accessions=geo_accessions,
                sra_accessions=sra_accessions,
                ena_accessions=ena_accessions,
                catalog=self.catalog
            )
            
            # Build unified dataset
            return UnifiedDataset(
                id=f"hca:{project_id}",
                accession=project_id,
                source=DataSource.HCA,
                source_url=browser_url,
                secondary_accession=secondary_accessions,
                omics_type=OmicsType.SINGLE_CELL,
                assay_types=assay_types,
                title=title,
                description=description[:2000] if description else None,  # Truncate long descriptions
                organism=list(organisms),
                sample_count=sample_count,
                sample_count_display=f"{sample_count:,}" if sample_count else "N/A",
                disease=list(diseases),
                tissue=list(tissues),
                cell_line=[],
                submission_date=submission_date,
                release_date=None,
                last_update=last_update,
                pubmed_ids=pubmed_ids,
                doi=doi,
                citation=None,
                contributors=contributors[:5],  # Limit contributors
                institution=institution,
                relevance_score=relevance_score,
                match_reasons=match_reasons[:5],  # Limit reasons
                quality_score=0.9,  # HCA is curated, high quality
                download_links=download_links,
                curation_level=CurationLevel.CURATED,  # HCA is expert-curated
                merged_sources=[],
                extensions={"hca": extension.model_dump()}
            )
            
        except Exception as e:
            logger.warning("Failed to convert HCA project to unified", error=str(e))
            return None
    
    def _calculate_relevance(
        self,
        title: str,
        description: str,
        organisms: List[str],
        tissues: List[str],
        diseases: List[str],
        query: Optional[str],
        search_terms: List[str]
    ) -> float:
        """Calculate relevance score for a project."""
        if not search_terms:
            return 0.5  # Default score if no query
        
        score = 0.0
        
        # Check title (highest weight)
        title_lower = title.lower()
        for term in search_terms:
            if term in title_lower:
                score += 0.3
        
        # Check description
        desc_lower = description.lower() if description else ""
        for term in search_terms:
            if term in desc_lower:
                score += 0.15
        
        # Check organisms
        for org in organisms:
            org_lower = org.lower()
            for term in search_terms:
                if term in org_lower:
                    score += 0.15
        
        # Check tissues
        for tissue in tissues:
            tissue_lower = tissue.lower()
            for term in search_terms:
                if term in tissue_lower:
                    score += 0.1
        
        # Check diseases
        for disease in diseases:
            disease_lower = disease.lower()
            for term in search_terms:
                if term in disease_lower:
                    score += 0.1
        
        # Cap at 1.0
        return min(score, 1.0)
    
    async def prewarm_cache(self) -> bool:
        """
        Pre-warm the projects cache.
        Called during startup to ensure cache is ready.
        """
        try:
            logger.info("Pre-warming HCA projects cache")
            await self._get_projects_cached(max_results=500)
            logger.info("HCA cache pre-warmed successfully", count=len(self._projects_cache or []))
            return True
        except Exception as e:
            logger.error("Failed to pre-warm HCA cache", error=str(e))
            return False


# Singleton instance
_hca_adapter: Optional[HCAAdapter] = None


def get_hca_adapter(catalog: Optional[str] = None) -> HCAAdapter:
    """Get or create HCA adapter singleton."""
    global _hca_adapter
    if _hca_adapter is None:
        _hca_adapter = HCAAdapter(catalog=catalog)
    return _hca_adapter
