"""
Unified Data Models for Multi-Omics Discovery
Supports datasets from GEO, PRIDE, ENCODE, MetaboLights, GDC, and more.
"""

from datetime import date, datetime
from enum import Enum
from typing import Optional, List, Any, Literal

from pydantic import BaseModel, ConfigDict, Field


# === ENUMS ===

class OmicsType(str, Enum):
    """Primary omics classification"""
    TRANSCRIPTOMICS = "transcriptomics"
    PROTEOMICS = "proteomics"
    METABOLOMICS = "metabolomics"
    EPIGENOMICS = "epigenomics"       # ChIP-seq, ATAC-seq, methylation
    GENOMICS = "genomics"              # WGS, WES
    SINGLE_CELL = "single_cell"
    METAGENOMICS = "metagenomics"
    MULTI_OMICS = "multi_omics"


class DataSource(str, Enum):
    """Supported data sources"""
    # Transcriptomics
    GEO = "geo"
    SRA = "sra"
    ARRAYEXPRESS = "arrayexpress"
    
    # Proteomics
    PRIDE = "pride"
    MASSIVE = "massive"
    PROTEOMEXCHANGE = "proteomexchange"
    
    # Metabolomics
    METABOLIGHTS = "metabolights"
    METABOLOMICS_WORKBENCH = "metabolomics_workbench"
    GNPS = "gnps"
    
    # Epigenomics
    ENCODE = "encode"
    ROADMAP = "roadmap"
    
    # Genomics
    ENA = "ena"
    EGA = "ega"
    
    # Single-cell
    HCA = "hca"
    CELLXGENE = "cellxgene"
    
    # Multi-omics / Cancer
    GDC = "gdc"
    CPTAC = "cptac"
    
    # Expression/Reference
    GTEX = "gtex"
    
    # Metagenomics
    MGNIFY = "mgnify"
    MGRAST = "mgrast"


class CurationLevel(str, Enum):
    """Dataset curation level."""
    CURATED = "curated"        # Expert-curated with standardized annotations
    COMMUNITY = "community"    # Community-submitted with some validation
    AUTO = "auto"              # Auto-indexed with minimal curation


class AssayType(str, Enum):
    """Common assay types across omics"""
    # Transcriptomics
    RNA_SEQ = "rna-seq"
    MRNA_SEQ = "mrna-seq"
    MICROARRAY = "microarray"
    SCRNA_SEQ = "scrna-seq"
    SNRNA_SEQ = "snrna-seq"
    SPATIAL_TRANSCRIPTOMICS = "spatial-transcriptomics"
    
    # Proteomics
    TMT = "tmt"
    ITRAQ = "itraq"
    SILAC = "silac"
    LABEL_FREE = "label-free"
    DIA = "dia"
    DDA = "dda"
    
    # Metabolomics
    LC_MS = "lc-ms"
    GC_MS = "gc-ms"
    CE_MS = "ce-ms"
    NMR = "nmr"
    
    # Epigenomics
    CHIP_SEQ = "chip-seq"
    ATAC_SEQ = "atac-seq"
    DNASE_SEQ = "dnase-seq"
    BISULFITE_SEQ = "bisulfite-seq"
    RRBS = "rrbs"
    WGBS = "wgbs"
    CUT_AND_RUN = "cut-and-run"
    CUT_AND_TAG = "cut-and-tag"
    MINT_CHIP = "mint-chip"
    HI_C = "hi-c"
    
    # Genomics
    WGS = "wgs"
    WES = "wes"
    TARGETED_SEQ = "targeted-seq"
    AMPLICON_SEQ = "amplicon-seq"
    
    # Multi-omics / Single-cell
    CITE_SEQ = "cite-seq"
    MULTIOME = "multiome"
    SMART_SEQ = "smart-seq"
    
    # Metagenomics
    AMPLICON_16S = "amplicon-16s"
    AMPLICON_18S = "amplicon-18s"
    AMPLICON_ITS = "amplicon-its"
    SHOTGUN_METAGENOMICS = "shotgun-metagenomics"
    METATRANSCRIPTOMICS = "metatranscriptomics"
    METAPROTEOMICS = "metaproteomics"
    
    # Other
    OTHER = "other"
    UNKNOWN = "unknown"


# === SOURCE METADATA ===

# === DOWNLOAD LINK ===

class DownloadFileType(str, Enum):
    """Types of downloadable files"""
    SOFT = "soft"              # GEO SOFT format
    MATRIX = "matrix"          # Expression matrix
    SUPPLEMENTARY = "supplementary"  # Supplementary files
    RAW = "raw"                # Raw/SRA data
    PROCESSED = "processed"    # Processed data files
    METADATA = "metadata"      # Metadata files
    OTHER = "other"


class DownloadLink(BaseModel):
    """A download link for a dataset file."""
    url: str = Field(..., description="Download URL (FTP or HTTPS)")
    file_type: DownloadFileType = Field(..., description="Type of file")
    file_name: Optional[str] = Field(None, description="File name if known")
    file_size_bytes: Optional[int] = Field(None, description="File size in bytes (fetched async on demand)")
    protocol: Literal["ftp", "https", "http"] = Field(default="https", description="URL protocol")
    description: Optional[str] = Field(None, description="Human-readable description")
    needs_refresh: bool = Field(default=False, description="Whether URL may expire and need refresh")
    
    model_config = ConfigDict(use_enum_values=True)


# === SOURCE METADATA ===

class SourceInfo(BaseModel):
    """Information about a data source"""
    source: DataSource
    name: str
    description: str
    url: str
    supported_omics: List[OmicsType]
    has_api: bool = True
    rate_limit: Optional[str] = None
    icon: str = "🧬"


# Source registry with metadata
SOURCE_REGISTRY: dict[DataSource, SourceInfo] = {
    DataSource.GEO: SourceInfo(
        source=DataSource.GEO,
        name="NCBI GEO",
        description="Gene Expression Omnibus - transcriptomics and more",
        url="https://www.ncbi.nlm.nih.gov/geo/",
        supported_omics=[OmicsType.TRANSCRIPTOMICS, OmicsType.EPIGENOMICS, OmicsType.SINGLE_CELL],
        rate_limit="3/sec (10 with API key)",
        icon="🧬"
    ),
    DataSource.PRIDE: SourceInfo(
        source=DataSource.PRIDE,
        name="PRIDE Archive",
        description="Proteomics Identifications Database",
        url="https://www.ebi.ac.uk/pride/",
        supported_omics=[OmicsType.PROTEOMICS],
        rate_limit="No strict limit",
        icon="🔬"
    ),
    DataSource.ENCODE: SourceInfo(
        source=DataSource.ENCODE,
        name="ENCODE",
        description="Encyclopedia of DNA Elements",
        url="https://www.encodeproject.org/",
        supported_omics=[OmicsType.EPIGENOMICS, OmicsType.TRANSCRIPTOMICS],
        rate_limit="Generous",
        icon="🧪"
    ),
    DataSource.METABOLIGHTS: SourceInfo(
        source=DataSource.METABOLIGHTS,
        name="MetaboLights",
        description="Metabolomics experiments and derived information",
        url="https://www.ebi.ac.uk/metabolights/",
        supported_omics=[OmicsType.METABOLOMICS],
        rate_limit="Generous",
        icon="⚗️"
    ),
    DataSource.GDC: SourceInfo(
        source=DataSource.GDC,
        name="GDC Data Portal",
        description="Genomic Data Commons - TCGA, CPTAC and more",
        url="https://portal.gdc.cancer.gov/",
        supported_omics=[OmicsType.MULTI_OMICS, OmicsType.GENOMICS, OmicsType.TRANSCRIPTOMICS],
        rate_limit="Generous",
        icon="🎗️"
    ),
    DataSource.CELLXGENE: SourceInfo(
        source=DataSource.CELLXGENE,
        name="CellxGene",
        description="Single-cell data portal by CZI",
        url="https://cellxgene.cziscience.com/",
        supported_omics=[OmicsType.SINGLE_CELL],
        rate_limit="Generous",
        icon="🔮"
    ),
    DataSource.MGNIFY: SourceInfo(
        source=DataSource.MGNIFY,
        name="MGnify",
        description="EBI Metagenomics - microbiome analysis platform",
        url="https://www.ebi.ac.uk/metagenomics/",
        supported_omics=[OmicsType.METAGENOMICS, OmicsType.TRANSCRIPTOMICS],
        rate_limit="5/sec",
        icon="🦠"
    ),
    DataSource.MGRAST: SourceInfo(
        source=DataSource.MGRAST,
        name="MG-RAST",
        description="Metagenomics RAST - metagenomic analysis server",
        url="https://www.mg-rast.org/",
        supported_omics=[OmicsType.METAGENOMICS, OmicsType.TRANSCRIPTOMICS],
        rate_limit="2/sec",
        icon="🌍"
    ),
    DataSource.HCA: SourceInfo(
        source=DataSource.HCA,
        name="Human Cell Atlas",
        description="Human Cell Atlas Data Portal - single-cell reference data",
        url="https://data.humancellatlas.org/",
        supported_omics=[OmicsType.SINGLE_CELL],
        rate_limit="Generous",
        icon="🧫"
    ),
    DataSource.METABOLOMICS_WORKBENCH: SourceInfo(
        source=DataSource.METABOLOMICS_WORKBENCH,
        name="Metabolomics Workbench",
        description="NIH Metabolomics Workbench - metabolomics data repository",
        url="https://www.metabolomicsworkbench.org/",
        supported_omics=[OmicsType.METABOLOMICS],
        rate_limit="10/sec",
        icon="🧪"
    ),
    DataSource.CPTAC: SourceInfo(
        source=DataSource.CPTAC,
        name="PDC/CPTAC",
        description="Proteomic Data Commons - CPTAC and NCI cancer proteomics",
        url="https://pdc.cancer.gov/",
        supported_omics=[OmicsType.PROTEOMICS, OmicsType.MULTI_OMICS],
        rate_limit="Generous (GraphQL)",
        icon="🎗️"
    ),
    DataSource.ENA: SourceInfo(
        source=DataSource.ENA,
        name="ENA",
        description="European Nucleotide Archive - sequencing data",
        url="https://www.ebi.ac.uk/ena/",
        supported_omics=[OmicsType.GENOMICS, OmicsType.TRANSCRIPTOMICS, OmicsType.METAGENOMICS],
        rate_limit="50/sec",
        icon="🧬"
    ),
    DataSource.MASSIVE: SourceInfo(
        source=DataSource.MASSIVE,
        name="MassIVE",
        description="Mass Spectrometry Interactive Virtual Environment",
        url="https://massive.ucsd.edu/",
        supported_omics=[OmicsType.PROTEOMICS, OmicsType.METABOLOMICS],
        rate_limit="Generous",
        icon="🔬"
    ),
    DataSource.GTEX: SourceInfo(
        source=DataSource.GTEX,
        name="GTEx",
        description="Genotype-Tissue Expression - human tissue expression atlas",
        url="https://gtexportal.org/",
        supported_omics=[OmicsType.TRANSCRIPTOMICS],
        rate_limit="Generous",
        icon="🧬"
    ),
}


# === UNIFIED DATASET ===

class UnifiedDataset(BaseModel):
    """
    Universal dataset representation that works across all omics types.
    This is the standardized format returned to the frontend.
    """
    
    # === IDENTIFIERS ===
    id: str = Field(..., description="Unique ID in format: {source}:{accession}")
    accession: str = Field(..., description="Native accession (GSE123456, PXD012345, etc.)")
    source: DataSource = Field(..., description="Origin database")
    source_url: str = Field(..., description="Direct link to dataset page")
    secondary_accession: List[str] = Field(default_factory=list, description="Related accessions (SRA, ENA, BioProject, etc.)")
    
    # === CLASSIFICATION ===
    omics_type: OmicsType = Field(..., description="Primary omics category")
    assay_types: List[AssayType] = Field(default_factory=list, description="Experimental techniques used")
    
    # === BASIC METADATA ===
    title: str = Field(..., description="Dataset title")
    description: Optional[str] = Field(None, description="Abstract or summary")
    organism: List[str] = Field(default_factory=list, description="Species studied")
    
    # === QUANTITATIVE METRICS ===
    sample_count: int = Field(0, description="Number of samples/runs")
    sample_count_display: str = Field("0", description="Formatted sample count")
    
    # === STUDY CONTEXT ===
    disease: List[str] = Field(default_factory=list, description="Disease/condition studied")
    tissue: List[str] = Field(default_factory=list, description="Tissue/cell type")
    cell_line: List[str] = Field(default_factory=list, description="Cell lines used")
    
    # === DATES ===
    submission_date: Optional[str] = Field(None, description="When submitted")
    release_date: Optional[str] = Field(None, description="When made public")
    last_update: Optional[str] = Field(None, description="Last modification")
    
    # === PUBLICATIONS ===
    pubmed_ids: List[str] = Field(default_factory=list)
    doi: Optional[str] = Field(None)
    citation: Optional[str] = Field(None)
    
    # === CONTRIBUTORS ===
    contributors: List[str] = Field(default_factory=list, description="Authors/submitters")
    institution: Optional[str] = Field(None)
    
    # === AI/SEARCH ===
    relevance_score: float = Field(0.0, ge=0.0, le=1.0, description="Computed relevance score")
    match_reasons: List[str] = Field(default_factory=list, description="Why this matched")
    analysis_tags: List[str] = Field(default_factory=list, description="Suggested analyses this dataset is suited for (from research question)")
    
    # === METADATA READINESS ===
    quality_score: Optional[float] = Field(
        None,
        ge=0.0,
        le=1.0,
        description="Metadata-readiness heuristic (legacy field name)",
    )
    
    # === DOWNLOAD LINKS ===
    download_links: List[DownloadLink] = Field(default_factory=list, description="Direct download links for dataset files")
    
    # === CURATION & MERGING ===
    curation_level: Optional[CurationLevel] = Field(None, description="Dataset curation level")
    merged_sources: List[DataSource] = Field(default_factory=list, description="Sources merged into this result (for cross-source dedup)")
    
    # === SOURCE-SPECIFIC EXTENSIONS ===
    extensions: dict = Field(default_factory=dict, description="Source-specific metadata")
    
    model_config = ConfigDict(use_enum_values=True)
    
    def model_post_init(self, __context) -> None:
        """Set computed fields after initialization."""
        if self.sample_count > 0 and self.sample_count_display == "0":
            self.sample_count_display = str(self.sample_count)


# === SOURCE-SPECIFIC EXTENSIONS ===

class GEOExtension(BaseModel):
    """GEO-specific metadata"""
    gse_id: str
    platform: Optional[str] = None
    platform_id: Optional[str] = None
    series_type: Optional[str] = None
    supplementary_files: List[str] = Field(default_factory=list)
    sra_ids: List[str] = Field(default_factory=list)


class PRIDEExtension(BaseModel):
    """PRIDE-specific metadata"""
    pxd_id: str
    instruments: List[str] = Field(default_factory=list)
    modifications: List[str] = Field(default_factory=list)  # PTMs
    quantification_method: Optional[str] = None
    software: List[str] = Field(default_factory=list)
    file_count: int = 0
    total_size_mb: float = 0.0
    experiment_types: List[str] = Field(default_factory=list)


class ENCODEExtension(BaseModel):
    """ENCODE-specific metadata"""
    experiment_id: str
    target: Optional[str] = None  # For ChIP-seq (e.g., H3K27ac)
    biosample_term: Optional[str] = None
    biosample_type: Optional[str] = None
    lab: Optional[str] = None
    award: Optional[str] = None
    status: str = "released"
    replicates: int = 0


class MetaboLightsExtension(BaseModel):
    """MetaboLights-specific metadata"""
    mtbls_id: str
    metabolites_identified: int = 0
    study_design: Optional[str] = None
    analytical_platform: Optional[str] = None
    study_factors: List[str] = Field(default_factory=list)


class GDCExtension(BaseModel):
    """GDC (TCGA/CPTAC) specific metadata"""
    project_id: str
    program: str  # TCGA, CPTAC, etc.
    primary_site: Optional[str] = None
    data_categories: List[str] = Field(default_factory=list)
    experimental_strategies: List[str] = Field(default_factory=list)
    case_count: int = 0
    file_count: int = 0


class CellxGeneExtension(BaseModel):
    """CellxGene Census specific metadata"""
    collection_id: str
    collection_name: Optional[str] = None
    tissue: Optional[str] = None
    tissue_ontology_term_id: Optional[str] = None
    disease: Optional[str] = None
    disease_ontology_term_id: Optional[str] = None
    assay: Optional[str] = None
    assay_ontology_term_id: Optional[str] = None
    cell_type_count: int = 0
    cell_count: int = 0
    is_primary_data: bool = True
    publisher_metadata: Optional[dict] = None


class MGnifyExtension(BaseModel):
    """MGnify (EBI Metagenomics) specific metadata"""
    study_id: str
    secondary_accession: Optional[str] = None  # ENA/SRA accession for deduplication
    biomes: List[str] = Field(default_factory=list)  # Environmental biome lineages
    pipeline_version: Optional[str] = None
    sample_count: int = 0
    analysis_count: int = 0
    experiment_type: Optional[str] = None  # amplicon, metagenomic, metatranscriptomic
    top_phyla: List[str] = Field(default_factory=list)  # Top taxonomic phyla found
    functional_categories: List[str] = Field(default_factory=list)  # GO terms, KEGG, etc.
    total_size_bytes: Optional[int] = None  # Total download size
    is_metatranscriptomics: bool = False  # Flag for metatranscriptomics studies


class MGRASTExtension(BaseModel):
    """MG-RAST specific metadata"""
    mgm_id: Optional[str] = None  # Metagenome ID (mgm4xxxxx)
    project_id: str  # Project ID (mgp#####)
    secondary_accession: Optional[str] = None  # ENA/SRA accession for deduplication
    biome: Optional[str] = None  # Environmental biome
    feature: Optional[str] = None  # Environmental feature
    material: Optional[str] = None  # Environmental material
    env_package: Optional[str] = None  # Environment package type
    sequence_type: Optional[str] = None  # amplicon, shotgun, etc.
    bp_count: Optional[int] = None  # Total base pairs
    sequence_count: Optional[int] = None  # Number of sequences
    top_phyla: List[str] = Field(default_factory=list)  # Top taxonomic phyla
    functional_categories: List[str] = Field(default_factory=list)  # Functional annotations
    total_size_bytes: Optional[int] = None  # Total download size
    is_metatranscriptomics: bool = False  # Flag for metatranscriptomics


class HCAExtension(BaseModel):
    """Human Cell Atlas (HCA) Data Portal specific metadata"""
    project_id: str  # HCA project UUID
    project_short_name: Optional[str] = None  # Short project name
    cell_count: int = 0  # Total cell count
    effective_cell_count: Optional[int] = None  # Matrix cell count (processed)
    organ: List[str] = Field(default_factory=list)  # Organs/tissues studied
    organ_part: List[str] = Field(default_factory=list)  # Organ parts
    development_stage: List[str] = Field(default_factory=list)  # Development stages
    library_construction_method: List[str] = Field(default_factory=list)  # e.g., 10x 3' v3
    nucleic_acid_source: List[str] = Field(default_factory=list)  # single cell, single nucleus
    file_formats: List[str] = Field(default_factory=list)  # Available file formats
    protocols: List[str] = Field(default_factory=list)  # Protocol types
    donor_count: int = 0  # Number of donors
    specimen_count: int = 0  # Number of specimens
    geo_accessions: List[str] = Field(default_factory=list)  # Linked GEO accessions
    sra_accessions: List[str] = Field(default_factory=list)  # Linked SRA accessions
    ena_accessions: List[str] = Field(default_factory=list)  # Linked ENA accessions
    catalog: str = "dcp56"  # HCA catalog version


class MetabolomicsWorkbenchExtension(BaseModel):
    """Metabolomics Workbench specific metadata"""
    study_id: str  # MW study ID (STxxxxxx)
    project_id: Optional[str] = None  # Project ID if available
    institute: Optional[str] = None  # Submitting institute
    analysis_type: Optional[str] = None  # MS, NMR, etc.
    ms_type: Optional[str] = None  # MS instrument type
    ion_mode: Optional[str] = None  # Positive/negative
    chromatography: Optional[str] = None  # LC, GC, etc.
    metabolite_count: int = 0  # Number of metabolites identified
    named_metabolite_count: int = 0  # Named metabolites
    study_type: Optional[str] = None  # Targeted/untargeted
    study_status: str = "public"  # Study status
    subject_count: int = 0  # Number of subjects
    factors: List[str] = Field(default_factory=list)  # Study factors
    metaboanalyst_link: Optional[str] = None  # Link to MetaboAnalyst if available
    metabominer_link: Optional[str] = None  # Link to MetaboMiner if available


class CPTACExtension(BaseModel):
    """CPTAC/PDC (Proteomic Data Commons) specific metadata"""
    pdc_study_id: str  # PDC study ID (PDC000xxx)
    study_submitter_id: Optional[str] = None  # Original submitter ID
    program_name: str  # CPTAC, TCGA, etc.
    project_name: Optional[str] = None  # Project within program
    analytical_fraction: Optional[str] = None  # Proteome, Phosphoproteome, Glycoproteome
    experiment_type: Optional[str] = None  # TMT10, TMT11, iTRAQ4, Label Free
    acquisition_type: Optional[str] = None  # DDA, DIA
    embargo_date: Optional[str] = None  # Data embargo date
    is_embargoed: bool = False  # Whether data is currently embargoed
    aliquots_count: int = 0  # Number of aliquots/samples
    cases_count: int = 0  # Number of cases/patients
    file_counts: List[dict] = Field(default_factory=list)  # {data_category, file_type, count}
    analytical_fractions: List[str] = Field(default_factory=list)  # All fractions in study


class ENAExtension(BaseModel):
    """ENA (European Nucleotide Archive) specific metadata"""
    study_accession: str  # Primary accession (PRJEBxxxxx or ERPxxxxxx)
    secondary_accession: Optional[str] = None  # SRA accession if exists (SRPxxxxxx)
    center_name: Optional[str] = None  # Submitting center
    broker_name: Optional[str] = None  # Broker if submitted via third party
    tax_id: Optional[int] = None  # NCBI Taxonomy ID
    library_strategy: Optional[str] = None  # RNA-Seq, WGS, AMPLICON, etc.
    library_source: Optional[str] = None  # GENOMIC, TRANSCRIPTOMIC, METAGENOMIC
    library_selection: Optional[str] = None  # RANDOM, PCR, etc.
    instrument_platform: Optional[str] = None  # ILLUMINA, OXFORD_NANOPORE, etc.
    sample_count: int = 0  # Number of samples
    run_count: int = 0  # Number of sequencing runs
    base_count: Optional[int] = None  # Total base pairs
    # Cross-references to other EBI resources
    arrayexpress_accession: Optional[str] = None  # ArrayExpress link
    biosamples_accession: Optional[str] = None  # BioSamples link
    eva_accession: Optional[str] = None  # European Variation Archive link


class MassIVEExtension(BaseModel):
    """MassIVE specific metadata"""
    msv_id: str  # MassIVE ID (MSVxxxxxxxxx)
    px_accession: Optional[str] = None  # ProteomeXchange ID (PXDxxxxxx)
    gnps_task_id: Optional[str] = None  # GNPS workflow task ID if applicable
    submitter: Optional[str] = None  # Submitter username
    file_count: int = 0  # Number of files
    total_size_bytes: Optional[int] = None  # Total dataset size
    instruments: List[str] = Field(default_factory=list)  # MS instruments used
    species: List[str] = Field(default_factory=list)  # Organisms studied
    modifications: List[str] = Field(default_factory=list)  # PTMs if applicable
    is_gnps_dataset: bool = False  # True if submitted via GNPS (metabolomics)
    is_reanalysis: bool = False  # True if this is a reanalysis of existing data
    has_redu_metadata: bool = False  # True if ReDU annotations exist
    molecular_networking_available: bool = False  # True if FBMN results exist
    gnps_visualization_url: Optional[str] = None  # Link to GNPS network visualization


# === SEARCH REQUEST/RESPONSE ===

class UnifiedSearchQuery(BaseModel):
    """Multi-omics search query input."""
    query: str = Field(..., min_length=3, max_length=500, description="Natural language search query")
    provider: str = Field(default="openai", description="AI provider to use")
    ranking_mode: Optional[Literal["keyword", "embedding", "llm", "hybrid"]] = Field(
        default=None,
        description="Optional per-request ranking override",
    )
    include_vector: Optional[bool] = Field(
        default=None,
        description="Optional per-request semantic-index retrieval override",
    )
    min_relevance_score: Optional[float] = Field(
        default=None,
        ge=0.0,
        le=1.0,
        description="Optional per-request relevance threshold override",
    )
    auto_index_results: bool = Field(
        default=True,
        description="Whether returned datasets may be added to the semantic index",
    )
    concept_expansion: bool = Field(
        default=False,
        description="Use deterministic curated biomedical concept expansion for this request",
    )
    
    # Filtering
    omics_types: Optional[List[OmicsType]] = Field(None, description="Filter by omics types (auto-detected if not specified)")
    sources: Optional[List[DataSource]] = Field(None, description="Specific sources to search")
    min_samples: Optional[int] = Field(None, ge=1, description="Minimum sample count")
    organism: Optional[str] = Field(None, description="Filter by organism")
    
    # Pagination - None means use the package runtime default.
    max_results: Optional[int] = Field(default=None, ge=1, le=200, description="Maximum results per source")
    sort_by: Optional[str] = Field(default=None, description="Sort order: relevance, date_desc, date_asc")


class UnifiedSearchStatus(BaseModel):
    """Status update during multi-omics search streaming."""
    stage: str = Field(..., description="Current stage of search")
    message: str = Field(..., description="Human-readable status message")
    progress: float = Field(default=0.0, ge=0.0, le=1.0, description="Progress percentage")
    current_source: Optional[DataSource] = None
    sources_completed: List[DataSource] = Field(default_factory=list)
    total_found: int = Field(default=0, description="Total results found so far")
    
    # Progressive facets (updated as sources return)
    facets: Optional["SearchFacets"] = Field(None, description="Current facet counts")


class FacetValue(BaseModel):
    """A single value within a facet, with count."""
    value: str = Field(..., description="The facet value (e.g., 'geo', 'proteomics', 'Homo sapiens')")
    display_name: str = Field(..., description="Human-readable display name")
    count: int = Field(0, description="Number of results with this value")
    would_be_count: Optional[int] = Field(None, description="Count if this filter were added (Amazon-style)")
    icon: Optional[str] = Field(None, description="Optional icon emoji for UI")
    color: Optional[str] = Field(None, description="Optional color class for UI")


class SearchFacets(BaseModel):
    """Facet data for search results - enables dynamic filtering."""
    
    # Source facet (GEO, PRIDE, etc.)
    sources: List[FacetValue] = Field(default_factory=list, description="Results by data source")
    
    # Omics type facet (transcriptomics, proteomics, etc.)
    omics_types: List[FacetValue] = Field(default_factory=list, description="Results by omics type")
    
    # Organism facet (Homo sapiens, Mus musculus, etc.)
    organisms: List[FacetValue] = Field(default_factory=list, description="Results by organism")
    
    # Optional additional facets
    tissues: List[FacetValue] = Field(default_factory=list, description="Results by tissue type")
    diseases: List[FacetValue] = Field(default_factory=list, description="Results by disease")
    
    # Metadata about facets
    last_updated: Optional[str] = Field(None, description="When facets were last calculated")
    is_complete: bool = Field(False, description="Whether all sources have returned")
    
    model_config = ConfigDict(use_enum_values=True)


class UnifiedSearchResponse(BaseModel):
    """Complete multi-omics search response."""
    query: str = Field(..., description="Original query")
    detected_omics_types: List[OmicsType] = Field(default_factory=list, description="AI-detected omics types")
    sources_searched: List[DataSource] = Field(default_factory=list, description="Sources that were searched")
    
    datasets: List[UnifiedDataset] = Field(default_factory=list, description="Matching datasets")
    
    total_results: int = Field(default=0, description="Total number of results")
    results_by_source: dict = Field(default_factory=dict, description="Count by source")
    results_by_omics: dict = Field(default_factory=dict, description="Count by omics type")
    
    # Faceted search data
    facets: Optional[SearchFacets] = Field(None, description="Facet counts for filtering UI")
    
    execution_time_ms: float = Field(default=0.0, description="Total execution time")
    ai_summary: Optional[str] = Field(None, description="AI-generated summary of results")
    
    model_config = ConfigDict(use_enum_values=True)


# === ORGANISM SYNONYM MAPPING ===

# Maps common organism names to their scientific names and synonyms
ORGANISM_SYNONYMS: dict[str, List[str]] = {
    # Human
    "homo sapiens": ["human", "homo sapiens", "h. sapiens", "h sapiens", "humans"],
    "human": ["homo sapiens", "human", "h. sapiens", "humans"],
    
    # Mouse
    "mus musculus": ["mouse", "mus musculus", "m. musculus", "mice", "murine"],
    "mouse": ["mus musculus", "mouse", "m. musculus", "mice", "murine"],
    
    # Rat
    "rattus norvegicus": ["rat", "rattus norvegicus", "r. norvegicus", "rats", "norway rat"],
    "rat": ["rattus norvegicus", "rat", "r. norvegicus", "rats"],
    
    # Zebrafish
    "danio rerio": ["zebrafish", "danio rerio", "d. rerio", "zebra fish"],
    "zebrafish": ["danio rerio", "zebrafish", "d. rerio", "zebra fish"],
    
    # Fruit fly
    "drosophila melanogaster": ["fruit fly", "drosophila", "drosophila melanogaster", "d. melanogaster", "fly"],
    "drosophila": ["drosophila melanogaster", "fruit fly", "drosophila", "d. melanogaster", "fly"],
    "fruit fly": ["drosophila melanogaster", "drosophila", "fruit fly", "fly"],
    
    # Worm
    "caenorhabditis elegans": ["c. elegans", "caenorhabditis elegans", "worm", "c elegans", "nematode"],
    "c. elegans": ["caenorhabditis elegans", "c. elegans", "worm", "c elegans", "nematode"],
    
    # Yeast
    "saccharomyces cerevisiae": ["yeast", "saccharomyces cerevisiae", "s. cerevisiae", "s cerevisiae", "budding yeast"],
    "yeast": ["saccharomyces cerevisiae", "yeast", "s. cerevisiae", "budding yeast"],
    
    # Arabidopsis
    "arabidopsis thaliana": ["arabidopsis", "arabidopsis thaliana", "a. thaliana", "thale cress"],
    "arabidopsis": ["arabidopsis thaliana", "arabidopsis", "a. thaliana", "thale cress"],
    
    # Pig
    "sus scrofa": ["pig", "sus scrofa", "s. scrofa", "swine", "porcine"],
    "pig": ["sus scrofa", "pig", "swine", "porcine"],
    
    # Dog
    "canis familiaris": ["dog", "canis familiaris", "c. familiaris", "canine"],
    "dog": ["canis familiaris", "dog", "canine", "canis lupus familiaris"],
    
    # Chicken
    "gallus gallus": ["chicken", "gallus gallus", "g. gallus"],
    "chicken": ["gallus gallus", "chicken", "g. gallus"],
    
    # Cow
    "bos taurus": ["cow", "bos taurus", "b. taurus", "cattle", "bovine"],
    "cow": ["bos taurus", "cow", "cattle", "bovine"],
}


def normalize_organism(organism: str) -> List[str]:
    """
    Get all synonyms for an organism name.
    Returns a list of all equivalent names to search for.
    """
    if not organism:
        return []
    
    organism_lower = organism.lower().strip()
    
    # Check if we have synonyms for this organism
    if organism_lower in ORGANISM_SYNONYMS:
        return ORGANISM_SYNONYMS[organism_lower]
    
    # Try partial matching
    for key, synonyms in ORGANISM_SYNONYMS.items():
        if organism_lower in key or key in organism_lower:
            return synonyms
        for syn in synonyms:
            if organism_lower in syn or syn in organism_lower:
                return synonyms
    
    # Return original if no match
    return [organism_lower]


def organisms_match(query_organism: str, dataset_organism: str) -> bool:
    """
    Check if a query organism matches a dataset organism.
    Uses synonym mapping for flexible matching.
    """
    if not query_organism or not dataset_organism:
        return True  # No filter = match all
    
    query_lower = query_organism.lower().strip()
    dataset_lower = dataset_organism.lower().strip()
    
    # Direct match
    if query_lower in dataset_lower or dataset_lower in query_lower:
        return True
    
    # Get synonyms for query organism
    query_synonyms = normalize_organism(query_lower)
    
    # Check if any synonym matches
    for synonym in query_synonyms:
        if synonym in dataset_lower or dataset_lower in synonym:
            return True
    
    return False


# === OMICS TYPE DETECTION ===

# Keywords for detecting omics types from natural language
OMICS_KEYWORDS: dict[OmicsType, List[str]] = {
    OmicsType.TRANSCRIPTOMICS: [
        "rna-seq", "rnaseq", "rna seq", "transcriptome", "transcriptomic",
        "gene expression", "mrna", "expression profiling", "microarray",
        "differential expression", "deseq", "edger"
    ],
    OmicsType.PROTEOMICS: [
        "proteom", "protein", "mass spec", "ms/ms", "tmt", "itraq", "silac",
        "label-free", "label free", "dia", "dda", "peptide", "lc-ms",
        "quantitative proteomics", "phosphoproteom", "ubiquitin"
    ],
    OmicsType.METABOLOMICS: [
        "metabolom", "metabolite", "metabono", "lipidom", "lipid",
        "small molecule", "nmr", "gc-ms", "lc-ms metabol", "untargeted metabol",
        "targeted metabol", "flux analysis"
    ],
    OmicsType.EPIGENOMICS: [
        "chip-seq", "chipseq", "atac-seq", "atacseq", "methylation", "bisulfite",
        "histone", "h3k", "dnase", "chromatin", "epigenom", "cut&run",
        "cut&tag", "rrbs", "wgbs", "enhancer", "promoter accessibility"
    ],
    OmicsType.GENOMICS: [
        "wgs", "wes", "whole genome", "whole exome", "variant", "mutation",
        "snp", "indel", "cnv", "structural variant", "genome sequencing",
        "dna sequencing", "targeted sequencing"
    ],
    OmicsType.SINGLE_CELL: [
        "single-cell", "single cell", "scrna", "snrna", "10x genomics",
        "dropseq", "drop-seq", "smart-seq", "cell atlas", "trajectory",
        "pseudotime", "clustering cells", "cite-seq", "multiome"
    ],
    OmicsType.METAGENOMICS: [
        "metagenom", "microbiome", "16s", "shotgun metagen", "amplicon",
        "gut microb", "oral microb", "skin microb", "environmental sample"
    ],
    OmicsType.MULTI_OMICS: [
        "multi-omics", "multiomics", "integrated omics", "tcga", "cptac",
        "proteogenomic", "multi-modal", "pan-cancer"
    ]
}


def detect_omics_types(query: str) -> List[OmicsType]:
    """
    Detect omics types from a natural language query.
    Returns list of likely omics types based on keyword matching.
    """
    query_lower = query.lower()
    detected = []
    
    for omics_type, keywords in OMICS_KEYWORDS.items():
        for keyword in keywords:
            if keyword in query_lower:
                if omics_type not in detected:
                    detected.append(omics_type)
                break
    
    # Default to transcriptomics if nothing detected (most common)
    if not detected:
        detected = [OmicsType.TRANSCRIPTOMICS]
    
    return detected


# === MAPPING OMICS TO SOURCES ===

OMICS_TO_SOURCES: dict[OmicsType, List[DataSource]] = {
    OmicsType.TRANSCRIPTOMICS: [DataSource.GEO, DataSource.ENA, DataSource.SRA, DataSource.ARRAYEXPRESS, DataSource.MGNIFY, DataSource.MGRAST],
    OmicsType.PROTEOMICS: [DataSource.PRIDE, DataSource.CPTAC, DataSource.MASSIVE, DataSource.GEO, DataSource.PROTEOMEXCHANGE],
    OmicsType.METABOLOMICS: [DataSource.METABOLIGHTS, DataSource.METABOLOMICS_WORKBENCH, DataSource.MASSIVE, DataSource.GNPS],
    OmicsType.EPIGENOMICS: [DataSource.ENCODE, DataSource.GEO, DataSource.ROADMAP],
    OmicsType.GENOMICS: [DataSource.ENA, DataSource.SRA, DataSource.EGA, DataSource.GEO],
    OmicsType.SINGLE_CELL: [DataSource.CELLXGENE, DataSource.HCA, DataSource.GEO],
    OmicsType.MULTI_OMICS: [DataSource.GDC, DataSource.CPTAC, DataSource.GEO],
    OmicsType.METAGENOMICS: [DataSource.MGNIFY, DataSource.MGRAST, DataSource.ENA, DataSource.GEO],
}


def get_sources_for_omics(omics_types: List[OmicsType]) -> List[DataSource]:
    """Get relevant data sources for given omics types."""
    sources = []
    for ot in omics_types:
        for source in OMICS_TO_SOURCES.get(ot, []):
            if source not in sources:
                sources.append(source)
    return sources
