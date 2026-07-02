"""
PRIDE Archive Adapter
Searches the PRIDE (Proteomics Identifications Database) for proteomics datasets.
API Docs: https://www.ebi.ac.uk/pride/ws/archive/v2/
"""

import asyncio
from typing import List, Optional, Any, Set
from urllib.parse import urlencode

import httpx
import structlog

from .base import (
    BaseSourceAdapter,
    COMMON_STOPWORDS,
    build_default_limits,
)
from ..models.unified import (
    DataSource,
    OmicsType,
    AssayType,
    UnifiedDataset,
    PRIDEExtension,
    DownloadLink,
    DownloadFileType,
)

logger = structlog.get_logger(__name__)

# PRIDE API base URL (configurable via environment in future)
PRIDE_API_BASE = "https://www.ebi.ac.uk/pride/ws/archive/v2"

# PRIDE-specific stopwords to add to common set
PRIDE_STOPWORDS: Set[str] = frozenset({
    'proteomics', 'proteomic', 'protein', 'proteins', 'proteome'
})


class PRIDEError(Exception):
    """PRIDE API error."""
    pass


class PRIDEAdapter(BaseSourceAdapter):
    """
    Adapter for PRIDE Archive (proteomics database).
    Uses PRIDE Archive REST API v2.
    """
    
    def __init__(self, timeout: float = 30.0):
        super().__init__()
        self.timeout = timeout
        self._client: Optional[httpx.AsyncClient] = None
    
    @property
    def source(self) -> DataSource:
        return DataSource.PRIDE
    
    @property
    def supported_omics(self) -> List[OmicsType]:
        return [OmicsType.PROTEOMICS]
    
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
    
    def _simplify_query(self, query: str) -> str:
        """
        Simplify a natural language query for PRIDE search.
        PRIDE works better with fewer, more specific terms, but we preserve
        important qualifying terms that distinguish between different types of data.
        """
        # Only remove truly generic filler words that don't add search value
        # Preserve protein-related terms since they help distinguish proteomics queries
        generic_stopwords = {
            'data', 'dataset', 'datasets', 'analysis', 'study', 'experiment', 
            'experiments', 'research', 'large', 'high', 'quality', 'based', 
            'using', 'from', 'with', 'the', 'and', 'for', 'that', 'this',
            'are', 'was', 'were', 'been', 'being', 'have', 'has', 'had',
            'show', 'shows', 'showing', 'find', 'finds', 'finding', 'look',
            'search', 'searching', 'get', 'want', 'need', 'looking'
        }
        
        # Terms that are valuable qualifiers for proteomics searches
        valuable_terms = {
            'protein', 'proteins', 'proteome', 'proteomics', 'proteomic',
            'phospho', 'phosphorylation', 'phosphoproteome', 'phosphoproteomics',
            'glyco', 'glycosylation', 'glycoproteome', 'glycoproteomics',
            'ubiquitin', 'ubiquitination', 'sumo', 'sumoylation',
            'acetyl', 'acetylation', 'methyl', 'methylation',
            'expression', 'interaction', 'binding', 'signaling', 'pathway',
            'quantitative', 'quantitation', 'label-free', 'tmt', 'itraq', 'silac',
            'secretome', 'membrane', 'nuclear', 'cytoplasmic', 'mitochondrial',
            'plasma', 'serum', 'tissue', 'cell', 'cells', 'line', 'lines'
        }
        
        # Tokenize and filter
        words = query.lower().split()
        filtered = []
        
        for word in words:
            # Remove punctuation for matching
            clean_word = word.strip('.,;:!?()[]{}"\'-')
            if len(clean_word) <= 2:
                continue
            if clean_word in generic_stopwords:
                continue
            filtered.append(clean_word)
        
        # Categorize terms by priority
        priority_terms = []
        valuable_qualifier_terms = []
        other_terms = []
        
        # High priority: disease, organism, tissue-specific terms
        high_priority_keywords = {
            'cancer', 'tumor', 'tumour', 'disease', 'syndrome', 'disorder', 
            'carcinoma', 'leukemia', 'lymphoma', 'melanoma', 'adenoma',
            'alzheimer', 'parkinson', 'diabetes', 'obesity', 'inflammation',
            'infection', 'viral', 'bacterial', 'covid', 'sars',
            'human', 'mouse', 'rat', 'zebrafish', 'drosophila', 'yeast',
            'brain', 'liver', 'heart', 'kidney', 'lung', 'breast', 'colon',
            'blood', 'bone', 'skin', 'muscle', 'nerve', 'pancreas'
        }
        
        for word in filtered:
            if any(pk in word for pk in high_priority_keywords):
                priority_terms.append(word)
            elif word in valuable_terms or any(vt in word for vt in valuable_terms):
                valuable_qualifier_terms.append(word)
            else:
                other_terms.append(word)
        
        # Build result: prioritize disease/organism, then valuable qualifiers, then others
        # Allow up to 6 terms for more specific searches
        result_terms = []
        result_terms.extend(priority_terms[:3])
        result_terms.extend(valuable_qualifier_terms[:2])
        remaining_slots = 6 - len(result_terms)
        if remaining_slots > 0:
            result_terms.extend(other_terms[:remaining_slots])
        
        result = ' '.join(result_terms) if result_terms else query.split()[0]
        logger.info("Simplified PRIDE query", original=query, simplified=result, 
                   priority=priority_terms, qualifiers=valuable_qualifier_terms)
        return result
    
    async def _fetch_with_retry(self, url: str) -> httpx.Response:
        """Fetch URL with retry logic.

        Thin wrapper over :meth:`BaseSourceAdapter._fetch_with_retry` so
        existing call sites keep their signature. The shared helper handles
        timeout / network / 5xx retries with exponential backoff.
        """
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
        Search PRIDE Archive for proteomics datasets.
        
        Args:
            query: Search query (keywords, disease, etc.)
            max_results: Maximum results to return
            organism: Optional organism filter (e.g., "Homo sapiens")
            
        Returns:
            List of UnifiedDataset objects
            
        Raises:
            PRIDEError: On API or network failures after retries
        """
        try:
            # Simplify query for PRIDE - extract key terms
            # PRIDE works better with simpler queries
            simplified_query = self._simplify_query(query)
            
            # Build PRIDE search URL
            params = {
                "keyword": simplified_query,
                "pageSize": min(max_results, 100),  # PRIDE max is 100 per page
                "page": 0,
                "sortDirection": "DESC",
                "sortFields": "submissionDate",
            }
            
            # Add organism filter if provided
            if organism:
                params["speciesFilter"] = organism
            
            url = f"{PRIDE_API_BASE}/search/projects?{urlencode(params)}"
            
            logger.info("Searching PRIDE", query=query, max_results=max_results)
            
            response = await self._fetch_with_retry(url)
            response.raise_for_status()
            
            data = response.json()
            
            # Parse projects from response - handle different API response formats
            # PRIDE API may return: 
            # 1. {"_embedded": {"compactprojects": [...]}} 
            # 2. A direct list of projects
            # 3. {"content": [...]} or similar
            if isinstance(data, list):
                projects = data
            elif isinstance(data, dict):
                projects = (
                    data.get("_embedded", {}).get("compactprojects") or
                    data.get("_embedded", {}).get("projects") or
                    data.get("content") or
                    data.get("projects") or
                    []
                )
            else:
                projects = []
            
            logger.info("PRIDE response parsed", project_count=len(projects), response_type=type(data).__name__)
            
            unified = []
            for project in projects[:max_results]:
                try:
                    ud = self._to_unified(project, query=query)
                    if ud:
                        unified.append(ud)
                except Exception as e:
                    logger.warning("Failed to parse PRIDE project", error=str(e))
                    continue
            
            # Sort by relevance score (descending)
            unified.sort(key=lambda x: x.relevance_score, reverse=True)
            
            logger.info("PRIDE search complete", results=len(unified))
            return unified
            
        except httpx.HTTPStatusError as e:
            logger.error("PRIDE API error", status=e.response.status_code, error=str(e))
            raise PRIDEError(f"PRIDE API error: {e.response.status_code}") from e
        except httpx.RequestError as e:
            logger.error("PRIDE request failed", error=str(e))
            raise PRIDEError(f"PRIDE request failed: {str(e)}") from e
    
    async def get_dataset(self, accession: str) -> Optional[UnifiedDataset]:
        """
        Fetch a single PRIDE project by accession with file details.
        
        Args:
            accession: PXD accession (e.g., "PXD012345")
            
        Returns:
            UnifiedDataset if found
        """
        try:
            # Normalize accession
            if not accession.upper().startswith("PXD"):
                accession = f"PXD{accession}"
            
            url = f"{PRIDE_API_BASE}/projects/{accession}"
            
            client = await self._get_client()
            response = await client.get(url)
            
            if response.status_code == 404:
                return None
            
            response.raise_for_status()
            project = response.json()
            
            # Fetch individual files for this project
            file_list = await self._fetch_project_files(accession)
            
            return self._to_unified_full(project, file_list=file_list)
            
        except httpx.HTTPStatusError as e:
            if e.response.status_code == 404:
                return None
            logger.error("PRIDE API error", accession=accession, error=str(e))
            return None
        except Exception as e:
            logger.error("Failed to fetch PRIDE project", accession=accession, error=str(e))
            return None
    
    async def _fetch_project_files(self, accession: str) -> List[dict]:
        """Fetch file list for a PRIDE project."""
        try:
            url = f"{PRIDE_API_BASE}/projects/{accession}/files"
            client = await self._get_client()
            response = await client.get(url)
            
            if response.status_code != 200:
                return []
            
            data = response.json()
            return data if isinstance(data, list) else []
        except Exception as e:
            logger.debug("Failed to fetch PRIDE files", accession=accession, error=str(e))
            return []
    
    def _to_unified(self, project: dict, query: Optional[str] = None) -> Optional[UnifiedDataset]:
        """
        Convert a compact PRIDE project to UnifiedDataset.
        Used for search results.
        
        Args:
            project: PRIDE project dict from API
            query: Optional search query for relevance scoring
        """
        accession = project.get("accession", "")
        if not accession:
            return None
        
        # Extract assay types from experiment types
        assay_types = self._detect_assay_types(project)
        
        # Helper to extract names from list of strings or dicts
        def extract_names(items):
            if not items:
                return []
            result = []
            for item in items:
                if isinstance(item, str) and item:
                    result.append(item)
                elif isinstance(item, dict):
                    name = item.get("name") or item.get("value") or item.get("accession")
                    if name:
                        result.append(str(name))
            return result
        
        # Helper to extract a single name from string or dict
        def extract_single_name(item):
            if isinstance(item, str):
                return item
            elif isinstance(item, dict):
                return item.get("name") or item.get("value") or item.get("accession") or ""
            return ""
        
        # Parse organisms
        organisms = project.get("organisms", []) or []
        organism_list = extract_names(organisms) if isinstance(organisms, list) else [str(organisms)] if organisms else []
        
        # Parse diseases
        diseases = project.get("diseases", []) or []
        disease_list = extract_names(diseases) if isinstance(diseases, list) else [str(diseases)] if diseases else []
        
        # Parse tissues/cell lines
        tissues = project.get("tissues", []) or []
        tissue_list = extract_names(tissues) if isinstance(tissues, list) else [str(tissues)] if tissues else []
        
        # Get quantification method (handle both string and dict)
        quant_methods = project.get("quantificationMethods", [])
        quant_method = None
        if quant_methods and isinstance(quant_methods, list) and quant_methods[0]:
            quant_method = extract_single_name(quant_methods[0])
        
        # Build PRIDE extension
        pride_ext = PRIDEExtension(
            pxd_id=accession,
            instruments=self._normalize_list(project.get("instruments", [])),
            modifications=self._normalize_list(project.get("ptms", [])),
            quantification_method=quant_method,
            experiment_types=self._normalize_list(project.get("experimentTypes", [])),
        )
        
        # Extract title and description for relevance scoring
        title = project.get("title", "")
        description = project.get("projectDescription") or project.get("description", "")
        
        # Compute relevance score if query provided (uses base class method with PRIDE stopwords)
        relevance = self._compute_keyword_relevance(query, title, description, PRIDE_STOPWORDS) if query else 0.0
        
        # Generate download links
        download_links = self._generate_download_links(accession, project)
        
        return UnifiedDataset(
            # Identifiers
            id=self.build_unified_id(accession),
            accession=accession,
            source=self.source,
            source_url=f"https://www.ebi.ac.uk/pride/archive/projects/{accession}",
            
            # Classification
            omics_type=OmicsType.PROTEOMICS,
            assay_types=assay_types,
            
            # Basic metadata
            title=title,
            description=description,
            organism=organism_list,
            
            # Quantitative - attempt to use reported counts, otherwise will be enriched later
            sample_count=int(project.get("numAssays") or project.get("numberOfSamples") or project.get("numberOfExperiments") or 0),
            sample_count_display=str(project.get("numAssays") or project.get("numberOfSamples") or project.get("numberOfExperiments") or "N/A"),
            
            # Study context
            disease=disease_list,
            tissue=tissue_list,
            
            # Dates
            submission_date=project.get("submissionDate"),
            release_date=project.get("publicationDate"),
            
            # Relevance
            relevance_score=relevance,
            
            # Publications
            pubmed_ids=[str(p) for p in project.get("references", []) if p] if project.get("references") else [],
            
            # Download links
            download_links=download_links,
            
            # Extensions
            extensions={"pride": pride_ext.model_dump()},
        )
    
    def _to_unified_full(self, project: dict, file_list: Optional[List[dict]] = None) -> Optional[UnifiedDataset]:
        """
        Convert a full PRIDE project to UnifiedDataset.
        Used for single project fetch.
        
        Args:
            project: PRIDE project dict from API
            file_list: Optional list of files from the files API
        """
        # Start with compact conversion
        unified = self._to_unified(project)
        if not unified:
            return None
        
        # Enhance with additional data from full project
        
        # Get sample count from assays or experimental design
        num_assays = len(project.get("assays", []))
        if num_assays > 0:
            unified.sample_count = num_assays
            unified.sample_count_display = str(num_assays)
        
        # Get file count and size
        if "pride" in unified.extensions:
            pride_ext_dict = unified.extensions["pride"]
            pride_ext_dict["file_count"] = len(file_list) if file_list else len(project.get("files", []))
            
            # Calculate total size from file_list if available
            if file_list:
                total_size = sum(f.get("fileSizeBytes", 0) for f in file_list)
            else:
                total_size = sum(f.get("fileSize", 0) for f in project.get("files", []))
            pride_ext_dict["total_size_mb"] = round(total_size / (1024 * 1024), 2)
            
            # Software used
            pride_ext_dict["software"] = project.get("softwares", [])
        
        # Regenerate download links with individual files if available
        if file_list:
            accession = unified.accession
            unified.download_links = self._generate_download_links(accession, project, file_list=file_list)
        
        # Get contributors - handle both list of dicts and list of strings
        submitters = project.get("submitters", []) or []
        if submitters:
            contributors = []
            for s in submitters:
                if isinstance(s, dict):
                    name = s.get('name', '')
                    if name:
                        affiliation = s.get('affiliation', '')
                        contributors.append(f"{name} ({affiliation})" if affiliation else name)
                elif isinstance(s, str) and s:
                    contributors.append(s)
            unified.contributors = contributors
        
        return unified
    
    def _detect_assay_types(self, project: dict) -> List[AssayType]:
        """
        Detect proteomics assay types from PRIDE project metadata.
        """
        assay_types = []
        
        # Check experiment types - handle both list of strings and list of dicts
        exp_types = project.get("experimentTypes", []) or []
        quant_methods = project.get("quantificationMethods", []) or []
        
        # Safely extract text from items (could be strings or dicts)
        def extract_text(items):
            texts = []
            for item in items:
                if isinstance(item, str):
                    texts.append(item)
                elif isinstance(item, dict):
                    texts.append(item.get("name", "") or item.get("value", "") or str(item))
                else:
                    texts.append(str(item))
            return " ".join(texts)
        
        text = extract_text(exp_types + quant_methods).lower()
        
        # Also check title and description
        title = project.get('title', '') or ''
        desc = project.get('projectDescription', '') or project.get('description', '') or ''
        text += f" {title} {desc}".lower()
        
        # TMT
        if "tmt" in text or "tandem mass tag" in text:
            assay_types.append(AssayType.TMT)
        
        # iTRAQ
        if "itraq" in text:
            assay_types.append(AssayType.ITRAQ)
        
        # SILAC
        if "silac" in text:
            assay_types.append(AssayType.SILAC)
        
        # Label-free
        if "label-free" in text or "label free" in text or "lfq" in text:
            assay_types.append(AssayType.LABEL_FREE)
        
        # DIA
        if "dia" in text or "data independent" in text or "swath" in text:
            assay_types.append(AssayType.DIA)
        
        # DDA (default for most proteomics)
        if "dda" in text or "data dependent" in text:
            assay_types.append(AssayType.DDA)
        
        # Default to DDA if nothing detected
        if not assay_types:
            assay_types.append(AssayType.DDA)
        
        return assay_types
    
    def _generate_download_links(self, accession: str, project: dict, file_list: Optional[List[dict]] = None) -> List[DownloadLink]:
        """
        Generate download links for a PRIDE project.
        
        PRIDE FTP structure:
        - ftp.pride.ebi.ac.uk/pride/data/archive/{year}/{month}/{PXD}/
        
        Args:
            accession: PXD accession (e.g., "PXD012345")
            project: Project dict from API (may contain submission date)
            file_list: Optional list of individual files from files API
            
        Returns:
            List of DownloadLink objects
        """
        links = []
        
        # If we have individual files, add them first
        if file_list:
            for file_info in file_list:
                file_name = file_info.get("fileName", "")
                file_size = file_info.get("fileSizeBytes")
                file_category = file_info.get("fileCategory", {})
                category_value = file_category.get("value", "") if isinstance(file_category, dict) else ""
                
                # Get FTP URL from publicFileLocations
                ftp_url = None
                public_locations = file_info.get("publicFileLocations", [])
                for loc in public_locations:
                    if isinstance(loc, dict) and "FTP" in loc.get("name", ""):
                        ftp_url = loc.get("value")
                        break
                
                if not ftp_url or not file_name:
                    continue
                
                # Determine file type based on category
                if category_value.upper() == "RAW":
                    dl_file_type = DownloadFileType.RAW
                elif category_value.upper() in ("SEARCH", "RESULT", "PEAK"):
                    dl_file_type = DownloadFileType.OTHER
                else:
                    dl_file_type = DownloadFileType.OTHER
                
                # Build description
                desc = f"{file_name}"
                if category_value:
                    desc += f" ({category_value})"
                
                links.append(DownloadLink(
                    url=ftp_url,
                    file_type=dl_file_type,
                    file_name=file_name,
                    file_size_bytes=file_size,
                    protocol="ftp",
                    description=desc,
                ))
        
        # Add utility links after individual files
        
        # Main archive page - always available
        links.append(DownloadLink(
            url=f"https://www.ebi.ac.uk/pride/archive/projects/{accession}",
            file_type=DownloadFileType.OTHER,
            protocol="https",
            description="PRIDE Archive project page (file list)",
        ))
        
        # ProteomeXchange viewer link
        links.append(DownloadLink(
            url=f"http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID={accession}",
            file_type=DownloadFileType.METADATA,
            protocol="http",
            description="ProteomeXchange dataset info",
        ))
        
        # FTP archive link - requires date extraction for exact path
        submission_date = project.get("submissionDate") or project.get("publicationDate")
        if submission_date:
            # Parse date like "2023-05-15" to extract year and month
            try:
                parts = submission_date.split("-")
                if len(parts) >= 2:
                    year, month = parts[0], parts[1]
                    ftp_url = f"ftp://ftp.pride.ebi.ac.uk/pride/data/archive/{year}/{month}/{accession}/"
                    links.append(DownloadLink(
                        url=ftp_url,
                        file_type=DownloadFileType.RAW,
                        protocol="ftp",
                        description="FTP archive (all raw files)",
                    ))
            except Exception:
                pass
        
        # Only add API endpoint if we don't have individual files (for search results)
        if not file_list:
            links.append(DownloadLink(
                url=f"https://www.ebi.ac.uk/pride/ws/archive/v2/files/byProject?accession={accession}",
                file_type=DownloadFileType.METADATA,
                protocol="https",
                description="PRIDE API - file list (JSON)",
                needs_refresh=True,  # Dynamic API response
            ))
        
        return links


# Factory function
_pride_adapter: Optional[PRIDEAdapter] = None


def get_pride_adapter() -> PRIDEAdapter:
    """Get or create PRIDE adapter instance."""
    global _pride_adapter
    if _pride_adapter is None:
        _pride_adapter = PRIDEAdapter()
    return _pride_adapter
