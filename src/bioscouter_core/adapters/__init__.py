"""
Multi-omics data source adapters.

Each adapter normalizes data from a specific source to the
:class:`~bioscouter_core.models.unified.UnifiedDataset` shape. The public
orchestrator in ``bioscouter_core.orchestrator`` registers all adapters and dispatches
searches to the appropriate ones based on the requested omics types.
"""

from .base import BaseSourceAdapter
from .cellxgene_adapter import CellxGeneAdapter
from .cptac_adapter import CPTACAdapter
from .ena_adapter import ENAAdapter
from .encode_adapter import ENCODEAdapter
from .gdc_adapter import GDCAdapter
from .geo_adapter import GEOAdapter
from .gtex_adapter import GTExAdapter
from .hca_adapter import HCAAdapter
from .massive_adapter import MassIVEAdapter
from .metabolights_adapter import MetaboLightsAdapter
from .metabolomics_workbench_adapter import MetabolomicsWorkbenchAdapter
from .mgnify_adapter import MGnifyAdapter
from .mgrast_adapter import MGRASTAdapter
from .pride_adapter import PRIDEAdapter

__all__ = [
    "BaseSourceAdapter",
    "CellxGeneAdapter",
    "CPTACAdapter",
    "ENAAdapter",
    "ENCODEAdapter",
    "GDCAdapter",
    "GEOAdapter",
    "GTExAdapter",
    "HCAAdapter",
    "MassIVEAdapter",
    "MetaboLightsAdapter",
    "MetabolomicsWorkbenchAdapter",
    "MGnifyAdapter",
    "MGRASTAdapter",
    "PRIDEAdapter",
]
