"""
NCBI GEO Search Service
Handles all interactions with NCBI E-utilities API via Biopython.
"""

import asyncio
import re
from typing import Optional
from xml.etree import ElementTree

import structlog
from Bio import Entrez
from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception_type

from bioscouter_core.core.config import get_settings
from bioscouter_core.models.schemas import Dataset, DownloadLink, DownloadFileType

logger = structlog.get_logger(__name__)

# Constants for sample count validation
MIN_VALID_SAMPLE_COUNT = 2
MAX_VALID_SAMPLE_COUNT = 10000


class NCBIError(Exception):
    """Base exception for NCBI-related errors."""
    pass


class NCBIRateLimitError(NCBIError):
    """Raised when NCBI rate limit is exceeded."""
    pass


class NCBIService:
    """Service for searching and retrieving data from NCBI GEO."""
    
    def __init__(self):
        self.settings = get_settings()
        
        # Configure Entrez
        Entrez.email = self.settings.ncbi_email
        Entrez.tool = self.settings.ncbi_tool_name
        if self.settings.ncbi_api_key:
            Entrez.api_key = self.settings.ncbi_api_key
    
    def _sanitize_query(self, query: str) -> str:
        """
        Sanitize user query to prevent injection attacks.
        Removes potentially dangerous characters and limits length.
        """
        # Remove characters that could affect NCBI query parsing
        sanitized = re.sub(r'["\[\]\(\)]', '', query)
        # Limit length to prevent abuse
        sanitized = sanitized[:500]
        return sanitized.strip()
    
    def _extract_gse_ids(self, query: str) -> list[str]:
        """
        Extract GSE accession IDs from query if present.
        Returns list of GSE IDs found (e.g., ['GSE123456', 'GSE789']).
        """
        # Match GSE followed by digits, case-insensitive
        pattern = r'\b(GSE\d+)\b'
        matches = re.findall(pattern, query, re.IGNORECASE)
        # Normalize to uppercase
        return [m.upper() for m in matches]
    
    @retry(
        stop=stop_after_attempt(3),
        wait=wait_exponential(multiplier=1, min=2, max=10),
        retry=retry_if_exception_type((IOError, TimeoutError)),
    )
    async def search(
        self,
        query: str,
        organism: Optional[str] = None,
        max_results: int = 100,
        soft_organism_filter: bool = True,
    ) -> list[str]:
        """
        Search NCBI GEO for datasets matching the query.
        Returns a list of GSE IDs.
        
        Args:
            query: Search query
            organism: Optional organism filter
            max_results: Maximum results to return
            soft_organism_filter: If True and organism is specified, also search
                                  without organism and merge results (improves recall)
        """
        # Check for exact GSE ID(s) in query - bypass keyword search
        gse_ids = self._extract_gse_ids(query)
        if gse_ids:
            logger.info("Detected GSE accession IDs in query, using direct lookup", gse_ids=gse_ids)
            # Return the GSE IDs directly - they'll be fetched by accession
            return gse_ids[:max_results]
        
        # Sanitize and validate input
        query = self._sanitize_query(query)
        if not query:
            raise NCBIError("Empty or invalid query after sanitization")
        
        # Build base search query (without organism)
        base_terms = [query, '"gse"[Entry Type]']
        base_query = " AND ".join(base_terms)
        
        # Build organism-filtered query if organism specified
        organism_query = None
        if organism:
            organism = self._sanitize_query(organism)
            if organism:
                organism_terms = base_terms + [f'"{organism}"[Organism]']
                organism_query = " AND ".join(organism_terms)
        
        logger.info("Searching NCBI GEO", 
                    base_query=base_query, 
                    organism_query=organism_query,
                    max_results=max_results,
                    soft_organism_filter=soft_organism_filter)
        
        try:
            loop = asyncio.get_event_loop()
            all_ids = []
            seen_ids = set()
            
            # Helper to run a single paginated search
            def _sync_search_paginated(search_query: str, target_count: int) -> list[str]:
                import socket
                old_timeout = socket.getdefaulttimeout()
                socket.setdefaulttimeout(30)
                collected_ids = []
                batch_size = min(100, target_count)  # NCBI recommends max 100 per request
                retstart = 0
                
                try:
                    while len(collected_ids) < target_count:
                        handle = Entrez.esearch(
                            db="gds",
                            term=search_query,
                            retmax=batch_size,
                            retstart=retstart,
                            usehistory="y",
                            sort="relevance",
                        )
                        results = Entrez.read(handle)
                        handle.close()
                        
                        batch_ids = results.get("IdList", [])
                        if not batch_ids:
                            break  # No more results
                        
                        collected_ids.extend(batch_ids)
                        retstart += len(batch_ids)
                        
                        # Check if we've exhausted results
                        total_count = int(results.get("Count", 0))
                        if retstart >= total_count:
                            break
                        
                        # Safety: limit pagination to avoid excessive API calls
                        if retstart >= 500:
                            logger.info("Reached pagination limit (500 IDs)")
                            break
                    
                    return collected_ids[:target_count]
                finally:
                    socket.setdefaulttimeout(old_timeout)
            
            # Strategy: If organism specified and soft filter enabled,
            # first search WITH organism, then WITHOUT to catch any missed results
            if organism_query and soft_organism_filter:
                # Primary search: with organism filter (prioritized)
                primary_ids = await loop.run_in_executor(
                    None, _sync_search_paginated, organism_query, max_results
                )
                for id_ in primary_ids:
                    if id_ not in seen_ids:
                        all_ids.append(id_)
                        seen_ids.add(id_)
                
                # Secondary search: without organism (fill remaining slots)
                remaining = max_results - len(all_ids)
                if remaining > 0:
                    secondary_ids = await loop.run_in_executor(
                        None, _sync_search_paginated, base_query, remaining + 20  # fetch extra for dedup
                    )
                    for id_ in secondary_ids:
                        if id_ not in seen_ids and len(all_ids) < max_results:
                            all_ids.append(id_)
                            seen_ids.add(id_)
                    
                    logger.info("Soft organism filter added extra results", 
                               primary=len(primary_ids), 
                               added_from_secondary=len(all_ids) - len(primary_ids))
            elif organism_query:
                # Hard organism filter (original behavior)
                all_ids = await loop.run_in_executor(
                    None, _sync_search_paginated, organism_query, max_results
                )
            else:
                # No organism filter
                all_ids = await loop.run_in_executor(
                    None, _sync_search_paginated, base_query, max_results
                )
            
            id_list = all_ids
            logger.info("NCBI search complete", total_found=len(id_list))
            
            return id_list
            
        except Exception as e:
            logger.error("NCBI search failed", error=str(e))
            raise NCBIError(f"Search failed: {str(e)}") from e
    
    @retry(
        stop=stop_after_attempt(3),
        wait=wait_exponential(multiplier=1, min=2, max=10),
        retry=retry_if_exception_type((IOError, TimeoutError)),
    )
    async def fetch_summaries(self, id_list: list[str]) -> list[Dataset]:
        """
        Fetch summary information for a list of GEO dataset IDs.
        Returns parsed Dataset objects.
        
        Handles both:
        - Numeric internal IDs from esearch (e.g., "200012345")
        - GSE accession strings from direct ID detection (e.g., "GSE12345")
        """
        if not id_list:
            return []
        
        # Separate GSE accessions from numeric IDs
        gse_accessions = [id_ for id_ in id_list if id_.upper().startswith("GSE")]
        numeric_ids = [id_ for id_ in id_list if not id_.upper().startswith("GSE")]
        
        datasets = []
        loop = asyncio.get_event_loop()
        
        # Handle GSE accessions via direct lookup
        if gse_accessions:
            logger.info("Fetching GSE accessions directly", count=len(gse_accessions))
            for gse_id in gse_accessions:
                try:
                    dataset = await self.fetch_by_gse_id(gse_id)
                    if dataset:
                        datasets.append(dataset)
                except Exception as e:
                    logger.warning("Failed to fetch GSE by accession", gse_id=gse_id, error=str(e))
        
        # Handle numeric IDs via esummary
        if numeric_ids:
            logger.info("Fetching NCBI summaries by internal ID", count=len(numeric_ids))
            
            try:
                def _sync_fetch():
                    handle = Entrez.esummary(
                        db="gds",
                        id=",".join(numeric_ids),
                        retmode="xml",
                    )
                    results = Entrez.read(handle)
                    handle.close()
                    return results
                
                results = await loop.run_in_executor(None, _sync_fetch)
                
                for item in results:
                    dataset = self._parse_summary(item)
                    if dataset:
                        datasets.append(dataset)
                        
            except Exception as e:
                logger.error("NCBI fetch failed", error=str(e))
                raise NCBIError(f"Fetch failed: {str(e)}") from e
        
        logger.info("Parsed summaries", count=len(datasets))
        return datasets
    
    def _parse_summary(self, item: dict) -> Optional[Dataset]:
        """Parse a single esummary result into a Dataset object."""
        try:
            # Extract GSE accession
            accession = item.get("Accession", "")
            if not accession.startswith("GSE"):
                # Skip non-series entries (platforms, samples, etc.)
                return None
            
            # Extract fields
            title = item.get("title", "")
            summary = item.get("summary", "")
            organism = item.get("taxon", "Unknown")
            platform = item.get("GPL", "")
            dataset_type = item.get("gdsType", "")
            submission_date = item.get("PDAT", "")
            
            # Extract sample count - this can be tricky
            sample_count = self._extract_sample_count(item, summary)
            
            # Generate download links
            download_links = self._generate_download_links(accession)
            
            return Dataset(
                gse_id=accession,
                title=title,
                summary=summary[:500] + "..." if len(summary) > 500 else summary,
                organism=organism,
                sample_count=sample_count,
                platform=platform,
                dataset_type=dataset_type,
                submission_date=submission_date,
                download_links=download_links,
            )
            
        except Exception as e:
            logger.warning("Failed to parse summary", error=str(e))
            return None
    
    def _generate_download_links(self, gse_id: str, sra_ids: list[str] = None) -> list[DownloadLink]:
        """
        Generate download links for a GEO dataset.
        
        GEO FTP structure:
        - ftp.ncbi.nlm.nih.gov/geo/series/{GSEnnn}nnn/{GSE}/
          - soft/ - SOFT formatted family file
          - matrix/ - Series matrix file (expression matrix)
          - suppl/ - Supplementary files (raw data, processed)
        """
        links = []
        
        # Extract GSE prefix for FTP path (e.g., "GSE123" becomes "GSE123nnn")
        gse_prefix = gse_id[:len(gse_id) - 3] if len(gse_id) > 6 else gse_id[:3]
        ftp_base = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{gse_prefix}nnn/{gse_id}"
        
        # SOFT file - Contains all metadata and data in SOFT format
        links.append(DownloadLink(
            url=f"{ftp_base}/soft/{gse_id}_family.soft.gz",
            file_type=DownloadFileType.SOFT,
            file_name=f"{gse_id}_family.soft.gz",
            protocol="https",
            description="SOFT formatted family file (metadata + data)",
        ))
        
        # Series Matrix file - Tab-delimited expression matrix
        # Note: Some datasets have multiple matrix files, but most have just one
        links.append(DownloadLink(
            url=f"{ftp_base}/matrix/{gse_id}_series_matrix.txt.gz",
            file_type=DownloadFileType.MATRIX,
            file_name=f"{gse_id}_series_matrix.txt.gz",
            protocol="https",
            description="Series matrix file (expression values)",
        ))
        
        # Supplementary files directory
        links.append(DownloadLink(
            url=f"{ftp_base}/suppl/",
            file_type=DownloadFileType.SUPPLEMENTARY,
            protocol="https",
            description="Supplementary files (raw/processed data)",
        ))
        
        # SRA raw data link if available
        if sra_ids and len(sra_ids) > 0:
            links.append(DownloadLink(
                url=f"https://www.ncbi.nlm.nih.gov/Traces/study/?acc={gse_id}",
                file_type=DownloadFileType.RAW,
                file_name="SRA Run Selector",
                protocol="https",
                description="SRA raw sequencing data (FASTQ)",
            ))
        
        return links
    
    def _extract_sample_count(self, item: dict, summary: str) -> Optional[int]:
        """
        Extract sample count from various sources.
        NCBI doesn't always have a clean sample count field.
        """
        # Try direct field first
        n_samples = item.get("n_samples")
        if n_samples and isinstance(n_samples, (int, str)):
            try:
                return int(n_samples)
            except (ValueError, TypeError):
                pass
        
        # Try Samples field (list of sample IDs)
        samples = item.get("Samples", [])
        if samples and isinstance(samples, list):
            return len(samples)
        
        # Try parsing from summary text
        # Look for patterns like "n=50", "50 samples", "50 subjects"
        patterns = [
            r'n\s*=\s*(\d+)',
            r'(\d+)\s*samples?',
            r'(\d+)\s*subjects?',
            r'(\d+)\s*patients?',
            r'(\d+)\s*specimens?',
            r'total\s+of\s+(\d+)',
        ]
        
        for pattern in patterns:
            match = re.search(pattern, summary, re.IGNORECASE)
            if match:
                try:
                    count = int(match.group(1))
                    if 1 < count < 10000:  # Sanity check
                        return count
                except (ValueError, IndexError):
                    continue
        
        return None
    
    async def search_and_fetch(
        self,
        query: str,
        organism: Optional[str] = None,
        max_results: int = 100,
        soft_organism_filter: bool = True,
    ) -> list[Dataset]:
        """
        Convenience method to search and fetch in one call.
        
        Args:
            query: Search query
            organism: Optional organism filter
            max_results: Maximum results (default increased to 100 for better recall)
            soft_organism_filter: If True, also search without organism to improve recall
        """
        id_list = await self.search(
            query=query,
            organism=organism,
            max_results=max_results,
            soft_organism_filter=soft_organism_filter,
        )
        if not id_list:
            return []
        
        return await self.fetch_summaries(id_list)

    @retry(
        stop=stop_after_attempt(3),
        wait=wait_exponential(multiplier=1, min=2, max=10),
        retry=retry_if_exception_type((IOError, TimeoutError)),
    )
    async def fetch_by_gse_id(self, gse_id: str) -> Optional[Dataset]:
        """
        Fetch a single dataset by its GSE accession number.
        This first searches for the GSE ID to get the internal ID,
        then fetches the full summary.
        """
        # Ensure we have the GSE prefix
        if not gse_id.startswith("GSE"):
            gse_id = f"GSE{gse_id}"
        
        logger.info("Fetching dataset by GSE ID", gse_id=gse_id)
        
        try:
            # Search for the specific GSE ID
            handle = Entrez.esearch(
                db="gds",
                term=f'{gse_id}[Accession]',
                retmax=1,
            )
            results = Entrez.read(handle)
            handle.close()
            
            id_list = results.get("IdList", [])
            if not id_list:
                logger.warning("GSE ID not found in NCBI", gse_id=gse_id)
                return None
            
            # Fetch the summary using the internal ID
            datasets = await self.fetch_summaries(id_list)
            return datasets[0] if datasets else None
            
        except Exception as e:
            logger.error("Failed to fetch by GSE ID", gse_id=gse_id, error=str(e))
            raise NCBIError(f"Failed to fetch {gse_id}: {str(e)}") from e


# Singleton instance
_ncbi_service: Optional[NCBIService] = None


def get_ncbi_service() -> NCBIService:
    """Get or create NCBI service instance."""
    global _ncbi_service
    if _ncbi_service is None:
        _ncbi_service = NCBIService()
    return _ncbi_service
