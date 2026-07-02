"""
Citation Generation Service

Generates properly formatted citations for datasets using standard citation formats.
No LLM required - uses programmatic formatting based on citation style guidelines.
"""

import re
from datetime import datetime
from typing import Optional, Dict, Any
from enum import Enum

import structlog

logger = structlog.get_logger(__name__)


class CitationFormat(str, Enum):
    """Supported citation formats."""
    APA = "apa"
    MLA = "mla"
    CHICAGO = "chicago"
    BIBTEX = "bibtex"
    RIS = "ris"
    ENDNOTE = "endnote"


class CitationService:
    """
    Service for generating properly formatted citations.
    
    Supports GEO, PRIDE, and other omics database datasets.
    Uses standard citation format rules - no LLM needed.
    """
    
    def __init__(self):
        self._format_handlers = {
            CitationFormat.APA: self._format_apa,
            CitationFormat.MLA: self._format_mla,
            CitationFormat.CHICAGO: self._format_chicago,
            CitationFormat.BIBTEX: self._format_bibtex,
            CitationFormat.RIS: self._format_ris,
            CitationFormat.ENDNOTE: self._format_endnote,
        }
    
    def generate_citation(
        self,
        accession: str,
        title: str,
        source: str,  # 'geo', 'pride', etc.
        format: str = "apa",
        organism: Optional[str] = None,
        submission_date: Optional[str] = None,
        authors: Optional[list[str]] = None,
        pubmed_ids: Optional[list[str]] = None,
        description: Optional[str] = None,
    ) -> str:
        """
        Generate a formatted citation for a dataset.
        
        Args:
            accession: Dataset accession (e.g., GSE123456, PXD012345)
            title: Dataset title
            source: Data source (geo, pride, etc.)
            format: Citation format (apa, mla, chicago, bibtex, ris, endnote)
            organism: Organism studied
            submission_date: Date submitted (various formats accepted)
            authors: List of author names (if available)
            pubmed_ids: Associated PubMed IDs
            description: Dataset description/abstract
            
        Returns:
            Formatted citation string
        """
        # Normalize format
        try:
            citation_format = CitationFormat(format.lower())
        except ValueError:
            citation_format = CitationFormat.APA
            logger.warning(f"Unknown citation format '{format}', falling back to APA")
        
        # Parse and normalize the date
        year = self._extract_year(submission_date)
        
        # Get source-specific info
        source_info = self._get_source_info(source, accession)
        
        # Build citation data
        citation_data = {
            "accession": accession,
            "title": self._clean_title(title),
            "year": year,
            "source": source.upper(),
            "source_name": source_info["name"],
            "url": source_info["url"],
            "publisher": source_info["publisher"],
            "authors": authors or ["Dataset Contributors"],
            "organism": organism,
            "pubmed_ids": pubmed_ids or [],
            "description": description,
            "access_date": datetime.now().strftime("%Y-%m-%d"),
        }
        
        # Get the appropriate formatter
        formatter = self._format_handlers.get(citation_format, self._format_apa)
        
        return formatter(citation_data)
    
    def _extract_year(self, date_str: Optional[str]) -> str:
        """Extract year from various date formats."""
        if not date_str:
            return "n.d."  # "no date" in citation style
        
        # Try to find a 4-digit year
        match = re.search(r'(\d{4})', str(date_str))
        if match:
            return match.group(1)
        
        return "n.d."
    
    def _clean_title(self, title: str) -> str:
        """Clean and normalize the title."""
        if not title:
            return "Untitled Dataset"
        
        # Remove extra whitespace
        title = " ".join(title.split())
        
        # Ensure title ends without period (we'll add punctuation as needed)
        title = title.rstrip(".")
        
        return title
    
    def _get_source_info(self, source: str, accession: str) -> Dict[str, str]:
        """Get source-specific metadata for citation."""
        source_lower = source.lower()
        
        sources = {
            "geo": {
                "name": "Gene Expression Omnibus",
                "publisher": "NCBI",
                "url": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}",
            },
            "pride": {
                "name": "PRIDE Archive",
                "publisher": "EMBL-EBI",
                "url": f"https://www.ebi.ac.uk/pride/archive/projects/{accession}",
            },
            "arrayexpress": {
                "name": "ArrayExpress",
                "publisher": "EMBL-EBI",
                "url": f"https://www.ebi.ac.uk/arrayexpress/experiments/{accession}",
            },
            "metabolights": {
                "name": "MetaboLights",
                "publisher": "EMBL-EBI",
                "url": f"https://www.ebi.ac.uk/metabolights/{accession}",
            },
            "encode": {
                "name": "ENCODE Project",
                "publisher": "ENCODE Consortium",
                "url": f"https://www.encodeproject.org/experiments/{accession}",
            },
            "sra": {
                "name": "Sequence Read Archive",
                "publisher": "NCBI",
                "url": f"https://www.ncbi.nlm.nih.gov/sra/{accession}",
            },
        }
        
        return sources.get(source_lower, {
            "name": source.upper(),
            "publisher": "Unknown Publisher",
            "url": f"https://example.com/{accession}",
        })
    
    def _format_author_list(self, authors: list[str], max_authors: int = 7) -> str:
        """Format author list according to citation style."""
        if not authors:
            return "Dataset Contributors"
        
        if len(authors) == 1:
            return authors[0]
        elif len(authors) == 2:
            return f"{authors[0]} & {authors[1]}"
        elif len(authors) <= max_authors:
            return ", ".join(authors[:-1]) + f", & {authors[-1]}"
        else:
            return ", ".join(authors[:6]) + ", ... & " + authors[-1]
    
    def _format_apa(self, data: Dict[str, Any]) -> str:
        """
        Format citation in APA 7th edition style.
        
        Format: Author. (Year). Title. Publisher. URL
        """
        authors = self._format_author_list(data["authors"])
        year = data["year"]
        title = data["title"]
        publisher = f"{data['source_name']}, {data['publisher']}"
        url = data["url"]
        
        return f"{authors}. ({year}). {title} [{data['source']} Accession: {data['accession']}]. {publisher}. {url}"
    
    def _format_mla(self, data: Dict[str, Any]) -> str:
        """
        Format citation in MLA 9th edition style.
        
        Format: "Title." Publisher, Year, URL.
        """
        title = data["title"]
        publisher = data["source_name"]
        year = data["year"]
        url = data["url"]
        access_date = datetime.strptime(data["access_date"], "%Y-%m-%d").strftime("%d %b. %Y")
        
        return f'"{title}." {publisher}, {data["publisher"]}, {year}, {url}. Accessed {access_date}.'
    
    def _format_chicago(self, data: Dict[str, Any]) -> str:
        """
        Format citation in Chicago Manual of Style (17th ed).
        
        Format: "Title." Publisher. Year. URL.
        """
        title = data["title"]
        publisher = f"{data['source_name']}, {data['publisher']}"
        year = data["year"]
        url = data["url"]
        
        return f'"{title}." {publisher}. {year}. {url}.'
    
    def _format_bibtex(self, data: Dict[str, Any]) -> str:
        """
        Format citation in BibTeX format.
        
        For use with LaTeX documents and reference managers.
        """
        # Create a valid BibTeX key
        key = f"{data['source'].lower()}_{data['accession'].lower()}"
        key = re.sub(r'[^a-z0-9_]', '', key)
        
        # Format authors for BibTeX (Last, First and Last, First)
        authors_str = " and ".join(data["authors"]) if data["authors"] else "Dataset Contributors"
        
        bibtex = f"""@misc{{{key},
  author = {{{authors_str}}},
  title = {{{data['title']}}},
  year = {{{data['year']}}},
  howpublished = {{{data['source_name']}}},
  publisher = {{{data['publisher']}}},
  note = {{{data['source']} Accession: {data['accession']}}},
  url = {{{data['url']}}},
  urldate = {{{data['access_date']}}}
}}"""
        return bibtex
    
    def _format_ris(self, data: Dict[str, Any]) -> str:
        """
        Format citation in RIS format.
        
        Compatible with EndNote, Mendeley, Zotero, etc.
        """
        lines = [
            "TY  - DATA",  # Type: Dataset
            f"TI  - {data['title']}",
            f"PY  - {data['year']}",
            f"PB  - {data['publisher']}",
            f"DB  - {data['source_name']}",
            f"AN  - {data['accession']}",
            f"UR  - {data['url']}",
            f"Y2  - {data['access_date']}",
        ]
        
        # Add authors
        for author in data["authors"]:
            lines.append(f"AU  - {author}")
        
        # Add organism if available
        if data.get("organism"):
            lines.append(f"KW  - {data['organism']}")
        
        # Add abstract if available
        if data.get("description"):
            lines.append(f"AB  - {data['description'][:500]}")
        
        lines.append("ER  - ")  # End of record
        
        return "\n".join(lines)
    
    def _format_endnote(self, data: Dict[str, Any]) -> str:
        """
        Format citation in EndNote XML format.
        
        Can be imported directly into EndNote.
        """
        authors_xml = "\n".join([
            f"        <author>{author}</author>"
            for author in data["authors"]
        ])
        
        xml = f"""<?xml version="1.0" encoding="UTF-8"?>
<xml>
  <records>
    <record>
      <ref-type name="Dataset">59</ref-type>
      <contributors>
        <authors>
{authors_xml}
        </authors>
      </contributors>
      <titles>
        <title>{data['title']}</title>
        <secondary-title>{data['source_name']}</secondary-title>
      </titles>
      <dates>
        <year>{data['year']}</year>
      </dates>
      <publisher>{data['publisher']}</publisher>
      <urls>
        <related-urls>
          <url>{data['url']}</url>
        </related-urls>
      </urls>
      <accession-num>{data['accession']}</accession-num>
      <database>{data['source']}</database>
    </record>
  </records>
</xml>"""
        return xml


# Singleton instance
_citation_service: Optional[CitationService] = None


def get_citation_service() -> CitationService:
    """Get or create citation service instance."""
    global _citation_service
    if _citation_service is None:
        _citation_service = CitationService()
    return _citation_service
