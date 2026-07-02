"""
Faceted Search Service

Calculates facet counts from search results.
Supports progressive updates and "would-be" counts.
"""

from typing import List, Dict, Optional
from collections import Counter
from datetime import datetime, timezone

import structlog

from bioscouter_core.models.unified import (
    UnifiedDataset,
    DataSource,
    OmicsType,
    SearchFacets,
    FacetValue,
    SOURCE_REGISTRY,
)

logger = structlog.get_logger(__name__)


# Icon and color mappings for facets
SOURCE_ICONS = {
    DataSource.GEO: ("🧬", "bg-blue-900/30 text-blue-400"),
    DataSource.PRIDE: ("🔬", "bg-purple-900/30 text-purple-400"),
    DataSource.ENCODE: ("🧪", "bg-green-900/30 text-green-400"),
    DataSource.METABOLIGHTS: ("⚗️", "bg-amber-900/30 text-amber-400"),
    DataSource.METABOLOMICS_WORKBENCH: ("🧪", "bg-orange-900/30 text-orange-400"),
    DataSource.GDC: ("🎗️", "bg-red-900/30 text-red-400"),
    DataSource.CELLXGENE: ("🔮", "bg-indigo-900/30 text-indigo-400"),
    DataSource.HCA: ("🧫", "bg-violet-900/30 text-violet-400"),
    DataSource.MGNIFY: ("🦠", "bg-teal-900/30 text-teal-400"),
    DataSource.MGRAST: ("🌍", "bg-cyan-900/30 text-cyan-400"),
    DataSource.CPTAC: ("🎗️", "bg-pink-900/30 text-pink-400"),
    DataSource.ENA: ("🧬", "bg-emerald-900/30 text-emerald-400"),
    DataSource.MASSIVE: ("🔬", "bg-fuchsia-900/30 text-fuchsia-400"),
    DataSource.GTEX: ("🧬", "bg-sky-900/30 text-sky-400"),
}

OMICS_COLORS = {
    OmicsType.TRANSCRIPTOMICS: ("📊", "bg-blue-500"),
    OmicsType.PROTEOMICS: ("🔬", "bg-purple-500"),
    OmicsType.METABOLOMICS: ("⚗️", "bg-amber-500"),
    OmicsType.EPIGENOMICS: ("🧪", "bg-green-500"),
    OmicsType.GENOMICS: ("🧬", "bg-cyan-500"),
    OmicsType.SINGLE_CELL: ("🔮", "bg-indigo-500"),
    OmicsType.METAGENOMICS: ("🦠", "bg-teal-500"),
    OmicsType.MULTI_OMICS: ("📈", "bg-pink-500"),
}


def calculate_facets(
    datasets: List[UnifiedDataset],
    sources_searched: List[DataSource],
    max_organisms: int = 20,
    calculate_would_be: bool = True,
) -> SearchFacets:
    """
    Calculate facet counts from a list of datasets.
    
    Args:
        datasets: List of datasets to analyze
        sources_searched: All sources that were searched (for zero-count display)
        max_organisms: Maximum organism values to return
        calculate_would_be: Whether to calculate "would-be" intersection counts
        
    Returns:
        SearchFacets with counts for each facet dimension
    """
    # Count by source
    source_counts: Counter = Counter()
    omics_counts: Counter = Counter()
    organism_counts: Counter = Counter()
    tissue_counts: Counter = Counter()
    disease_counts: Counter = Counter()
    
    for dataset in datasets:
        # Source count
        source_counts[dataset.source] += 1
        
        # Omics type count
        omics_counts[dataset.omics_type] += 1
        
        # Organism count (normalize to first if multiple)
        for org in dataset.organism[:2]:  # Limit to first 2 organisms per dataset
            normalized = _normalize_organism_name(org)
            organism_counts[normalized] += 1
        
        # Tissue count
        for tissue in dataset.tissue[:2]:
            tissue_counts[tissue.lower().strip()] += 1
        
        # Disease count
        for disease in dataset.disease[:2]:
            disease_counts[disease.lower().strip()] += 1
    
    # Build facet values for sources
    source_facets = []
    for source in sources_searched:
        count = source_counts.get(source, 0)
        source_info = SOURCE_REGISTRY.get(source)
        icon, color = SOURCE_ICONS.get(source, ("📁", "bg-gray-500"))
        
        source_facets.append(FacetValue(
            value=source.value if hasattr(source, 'value') else str(source),
            display_name=source_info.name if source_info else str(source),
            count=count,
            icon=icon,
            color=color,
        ))
    
    # Sort sources by count (descending), then alphabetically
    source_facets.sort(key=lambda x: (-x.count, x.display_name))
    
    # Build facet values for omics types
    omics_facets = []
    for omics_type in OmicsType:
        count = omics_counts.get(omics_type, 0)
        icon, color = OMICS_COLORS.get(omics_type, ("📊", "bg-gray-500"))
        
        omics_facets.append(FacetValue(
            value=omics_type.value if hasattr(omics_type, 'value') else str(omics_type),
            display_name=_format_omics_name(omics_type),
            count=count,
            icon=icon,
            color=color,
        ))
    
    # Sort omics by count (descending)
    omics_facets.sort(key=lambda x: (-x.count, x.display_name))
    
    # Build facet values for organisms (top N)
    organism_facets = []
    for org, count in organism_counts.most_common(max_organisms):
        organism_facets.append(FacetValue(
            value=org,
            display_name=_format_organism_name(org),
            count=count,
            icon="🧬",
        ))
    
    # Build tissue facets (top 10)
    tissue_facets = []
    for tissue, count in tissue_counts.most_common(10):
        tissue_facets.append(FacetValue(
            value=tissue,
            display_name=tissue.title(),
            count=count,
        ))
    
    # Build disease facets (top 10)
    disease_facets = []
    for disease, count in disease_counts.most_common(10):
        disease_facets.append(FacetValue(
            value=disease,
            display_name=disease.title(),
            count=count,
        ))
    
    return SearchFacets(
        sources=source_facets,
        omics_types=omics_facets,
        organisms=organism_facets,
        tissues=tissue_facets,
        diseases=disease_facets,
        last_updated=datetime.now(timezone.utc).isoformat(),
        is_complete=True,
    )


def calculate_progressive_facets(
    current_facets: Optional[SearchFacets],
    new_datasets: List[UnifiedDataset],
    source: DataSource,
    sources_searched: List[DataSource],
) -> SearchFacets:
    """
    Update facet counts progressively as new source results arrive.
    
    Args:
        current_facets: Existing facets (or None for first update)
        new_datasets: New datasets from one source
        source: The source that just returned
        sources_searched: All sources being searched
        
    Returns:
        Updated SearchFacets
    """
    if current_facets is None:
        # First update - calculate from scratch
        facets = calculate_facets(new_datasets, sources_searched)
        facets.is_complete = False
        return facets
    
    # Merge counts from new datasets
    # Create a copy and update
    source_counts = {fv.value: fv.count for fv in current_facets.sources}
    omics_counts = {fv.value: fv.count for fv in current_facets.omics_types}
    organism_counts = {fv.value: fv.count for fv in current_facets.organisms}
    
    for dataset in new_datasets:
        # Update source count
        src_val = dataset.source.value if hasattr(dataset.source, 'value') else str(dataset.source)
        source_counts[src_val] = source_counts.get(src_val, 0) + 1
        
        # Update omics count
        omics_val = dataset.omics_type.value if hasattr(dataset.omics_type, 'value') else str(dataset.omics_type)
        omics_counts[omics_val] = omics_counts.get(omics_val, 0) + 1
        
        # Update organism count
        for org in dataset.organism[:2]:
            normalized = _normalize_organism_name(org)
            organism_counts[normalized] = organism_counts.get(normalized, 0) + 1
    
    # Rebuild facet lists with updated counts
    updated_sources = []
    for fv in current_facets.sources:
        updated_sources.append(FacetValue(
            value=fv.value,
            display_name=fv.display_name,
            count=source_counts.get(fv.value, 0),
            icon=fv.icon,
            color=fv.color,
        ))
    updated_sources.sort(key=lambda x: (-x.count, x.display_name))
    
    updated_omics = []
    for fv in current_facets.omics_types:
        updated_omics.append(FacetValue(
            value=fv.value,
            display_name=fv.display_name,
            count=omics_counts.get(fv.value, 0),
            icon=fv.icon,
            color=fv.color,
        ))
    updated_omics.sort(key=lambda x: (-x.count, x.display_name))
    
    # For organisms, merge and re-sort
    all_organisms = {}
    for fv in current_facets.organisms:
        all_organisms[fv.value] = fv.count
    for org, count in organism_counts.items():
        all_organisms[org] = count
    
    updated_organisms = []
    for org, count in sorted(all_organisms.items(), key=lambda x: -x[1])[:20]:
        updated_organisms.append(FacetValue(
            value=org,
            display_name=_format_organism_name(org),
            count=count,
            icon="🧬",
        ))
    
    return SearchFacets(
        sources=updated_sources,
        omics_types=updated_omics,
        organisms=updated_organisms,
        tissues=current_facets.tissues,
        diseases=current_facets.diseases,
        last_updated=datetime.now(timezone.utc).isoformat(),
        is_complete=False,
    )


def calculate_would_be_counts(
    datasets: List[UnifiedDataset],
    current_filters: Dict[str, List[str]],
) -> Dict[str, Dict[str, int]]:
    """
    Calculate "would-be" counts for facet values.
    
    This shows how many results would exist if a filter value were added,
    like Amazon's product filtering.
    
    Args:
        datasets: All datasets (before filtering)
        current_filters: Currently applied filters {facet_name: [values]}
        
    Returns:
        Dict mapping facet_name -> value -> would_be_count
    """
    would_be: Dict[str, Dict[str, int]] = {
        "sources": {},
        "omics_types": {},
        "organisms": {},
    }
    
    for dataset in datasets:
        # Check if dataset passes current filters (excluding each facet in turn)
        passes_without_source = _passes_filters(dataset, current_filters, exclude="sources")
        passes_without_omics = _passes_filters(dataset, current_filters, exclude="omics_types")
        passes_without_organism = _passes_filters(dataset, current_filters, exclude="organisms")
        
        # Count for source
        if passes_without_source:
            src = dataset.source.value if hasattr(dataset.source, 'value') else str(dataset.source)
            would_be["sources"][src] = would_be["sources"].get(src, 0) + 1
        
        # Count for omics type
        if passes_without_omics:
            omics = dataset.omics_type.value if hasattr(dataset.omics_type, 'value') else str(dataset.omics_type)
            would_be["omics_types"][omics] = would_be["omics_types"].get(omics, 0) + 1
        
        # Count for organism
        if passes_without_organism:
            for org in dataset.organism[:2]:
                normalized = _normalize_organism_name(org)
                would_be["organisms"][normalized] = would_be["organisms"].get(normalized, 0) + 1
    
    return would_be


def _passes_filters(
    dataset: UnifiedDataset,
    filters: Dict[str, List[str]],
    exclude: Optional[str] = None,
) -> bool:
    """Check if a dataset passes all filters (optionally excluding one facet)."""
    for facet_name, values in filters.items():
        if facet_name == exclude or not values:
            continue
            
        if facet_name == "sources":
            src = dataset.source.value if hasattr(dataset.source, 'value') else str(dataset.source)
            if src not in values:
                return False
        elif facet_name == "omics_types":
            omics = dataset.omics_type.value if hasattr(dataset.omics_type, 'value') else str(dataset.omics_type)
            if omics not in values:
                return False
        elif facet_name == "organisms":
            dataset_orgs = [_normalize_organism_name(o) for o in dataset.organism]
            if not any(o in values for o in dataset_orgs):
                return False
    
    return True


def _normalize_organism_name(org: str) -> str:
    """Normalize organism name for consistent counting."""
    if not org:
        return "Unknown"
    
    # Basic normalization
    normalized = org.strip().lower()
    
    # Common mappings
    mappings = {
        "homo sapiens": "Homo sapiens",
        "human": "Homo sapiens",
        "mus musculus": "Mus musculus",
        "mouse": "Mus musculus",
        "rattus norvegicus": "Rattus norvegicus",
        "rat": "Rattus norvegicus",
        "danio rerio": "Danio rerio",
        "zebrafish": "Danio rerio",
        "drosophila melanogaster": "Drosophila melanogaster",
        "fruit fly": "Drosophila melanogaster",
        "caenorhabditis elegans": "Caenorhabditis elegans",
        "c. elegans": "Caenorhabditis elegans",
        "saccharomyces cerevisiae": "Saccharomyces cerevisiae",
        "yeast": "Saccharomyces cerevisiae",
        "arabidopsis thaliana": "Arabidopsis thaliana",
    }
    
    if normalized in mappings:
        return mappings[normalized]
    
    # Title case the original
    return org.strip().title() if org else "Unknown"


def _format_omics_name(omics_type: OmicsType) -> str:
    """Format omics type for display."""
    name_map = {
        OmicsType.TRANSCRIPTOMICS: "Transcriptomics",
        OmicsType.PROTEOMICS: "Proteomics",
        OmicsType.METABOLOMICS: "Metabolomics",
        OmicsType.EPIGENOMICS: "Epigenomics",
        OmicsType.GENOMICS: "Genomics",
        OmicsType.SINGLE_CELL: "Single-cell",
        OmicsType.METAGENOMICS: "Metagenomics",
        OmicsType.MULTI_OMICS: "Multi-omics",
    }
    return name_map.get(omics_type, str(omics_type))


def _format_organism_name(org: str) -> str:
    """Format organism name for display."""
    if not org:
        return "Unknown"
    
    # Already formatted
    if org[0].isupper():
        return org
    
    return org.title()
