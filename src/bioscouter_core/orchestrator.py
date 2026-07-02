"""Standalone federated search orchestration for the public core package.

This module intentionally excludes BioScouter production services such as
authentication, credits, admin settings, commercial routing, vector stores,
and LLM services. It keeps the scientific path needed for reproducibility:
source selection, parallel adapter calls, normalization, deduplication,
metadata-readiness scoring, deterministic concept expansion, faceting, and
traceable result summaries.
"""

from __future__ import annotations

import asyncio
import time
from dataclasses import dataclass
from typing import Any, Iterable

import structlog

from .adapters.cellxgene_adapter import get_cellxgene_adapter
from .adapters.cptac_adapter import get_cptac_adapter
from .adapters.ena_adapter import get_ena_adapter
from .adapters.encode_adapter import get_encode_adapter
from .adapters.gdc_adapter import get_gdc_adapter
from .adapters.geo_adapter import get_geo_adapter
from .adapters.gtex_adapter import get_gtex_adapter
from .adapters.hca_adapter import get_hca_adapter
from .adapters.massive_adapter import get_massive_adapter
from .adapters.metabolights_adapter import get_metabolights_adapter
from .adapters.metabolomics_workbench_adapter import get_metabolomics_workbench_adapter
from .adapters.mgnify_adapter import get_mgnify_adapter
from .adapters.mgrast_adapter import get_mgrast_adapter
from .adapters.pride_adapter import get_pride_adapter
from .core.config import get_settings
from .models.unified import (
    DataSource,
    OmicsType,
    UnifiedDataset,
    UnifiedSearchQuery,
    UnifiedSearchResponse,
    detect_omics_types,
    get_sources_for_omics,
)
from .services.facets import calculate_facets
from .services.quality_score import calculate_and_set_quality
from .services.query_ontology import expand_query

logger = structlog.get_logger(__name__)

VALID_SORT_OPTIONS = frozenset({"relevance", "date_desc", "date_asc", "sample_count"})

ADAPTER_REGISTRY = (
    (DataSource.GEO, get_geo_adapter),
    (DataSource.PRIDE, get_pride_adapter),
    (DataSource.METABOLIGHTS, get_metabolights_adapter),
    (DataSource.ENCODE, get_encode_adapter),
    (DataSource.GDC, get_gdc_adapter),
    (DataSource.CELLXGENE, get_cellxgene_adapter),
    (DataSource.MGNIFY, get_mgnify_adapter),
    (DataSource.MGRAST, get_mgrast_adapter),
    (DataSource.HCA, get_hca_adapter),
    (DataSource.METABOLOMICS_WORKBENCH, get_metabolomics_workbench_adapter),
    (DataSource.CPTAC, get_cptac_adapter),
    (DataSource.ENA, get_ena_adapter),
    (DataSource.MASSIVE, get_massive_adapter),
    (DataSource.GTEX, get_gtex_adapter),
)

SOURCE_PRIORITY = {
    DataSource.PRIDE: 11,
    DataSource.METABOLIGHTS: 10,
    DataSource.HCA: 9,
    DataSource.GEO: 8,
    DataSource.GDC: 8,
    DataSource.CPTAC: 7,
    DataSource.CELLXGENE: 7,
    DataSource.ENA: 6,
    DataSource.METABOLOMICS_WORKBENCH: 6,
    DataSource.MGNIFY: 5,
    DataSource.MASSIVE: 4,
    DataSource.MGRAST: 4,
}


@dataclass(slots=True)
class SearchConfig:
    """Runtime knobs for standalone public-core searches."""

    max_results: int = 50
    min_relevance_score: float = 0.0
    sort_by: str = "relevance"
    source_timeout_seconds: float = 30.0
    concept_expansion: bool = False


class BioScouterCoreSearch:
    """Run federated searches across the public BioScouter adapter set."""

    def __init__(self, config: SearchConfig | None = None, adapters: dict[DataSource, Any] | None = None) -> None:
        settings = get_settings()
        self.config = config or SearchConfig(
            max_results=settings.default_max_results,
            min_relevance_score=settings.min_relevance_score,
            source_timeout_seconds=settings.source_timeout_seconds,
        )
        self._adapters = adapters or self._register_adapters()

    @property
    def adapters(self) -> dict[DataSource, Any]:
        """Registered adapter instances keyed by data source."""

        return dict(self._adapters)

    def available_sources(self) -> list[DataSource]:
        """Return sources supported by the current adapter registry."""

        return list(self._adapters.keys())

    async def close(self) -> None:
        """Close adapter HTTP clients where adapters expose ``close``."""

        for adapter in self._adapters.values():
            close = getattr(adapter, "close", None)
            if close is None:
                continue
            result = close()
            if asyncio.iscoroutine(result):
                await result

    async def search(self, query: UnifiedSearchQuery | str, **overrides: Any) -> UnifiedSearchResponse:
        """Execute a federated search and return normalized results.

        Args:
            query: A plain query string or ``UnifiedSearchQuery``.
            **overrides: Optional fields such as ``sources``, ``omics_types``,
                ``organism``, ``max_results``, ``sort_by``,
                ``min_relevance_score``, ``source_timeout_seconds``, and
                ``concept_expansion``.
        """

        start = time.perf_counter()
        request = self._coerce_query(query, overrides)
        config = self._coerce_config(request, overrides)
        detected_omics = request.omics_types or detect_omics_types(request.query)
        sources_to_search = self._select_sources(detected_omics, request.sources)

        expansion = expand_query(request.query) if config.concept_expansion else None
        search_text = expansion.expanded_query if expansion and expansion.changed else request.query

        per_source = await self._search_sources_parallel(
            search_text=search_text,
            sources=sources_to_search,
            max_results=config.max_results,
            organism=request.organism,
            timeout_seconds=config.source_timeout_seconds,
        )
        merged = self._merge_results(per_source)

        if config.min_relevance_score > 0:
            merged = [d for d in merged if (d.relevance_score or 0.0) >= config.min_relevance_score]

        for dataset in merged:
            calculate_and_set_quality(dataset)

        sorted_results = self._apply_sort(merged, config.sort_by)[: config.max_results]
        source_counts, omics_counts = self._summarize(sorted_results)
        facets = calculate_facets(sorted_results, sources_to_search)

        return UnifiedSearchResponse(
            query=request.query,
            detected_omics_types=detected_omics,
            sources_searched=sources_to_search,
            datasets=sorted_results,
            total_results=len(sorted_results),
            results_by_source=source_counts,
            results_by_omics=omics_counts,
            facets=facets,
            execution_time_ms=round((time.perf_counter() - start) * 1000, 2),
            ai_summary=None,
        )

    async def _search_sources_parallel(
        self,
        *,
        search_text: str,
        sources: list[DataSource],
        max_results: int,
        organism: str | None,
        timeout_seconds: float,
    ) -> dict[DataSource, list[UnifiedDataset]]:
        tasks = [
            self._search_one_source(
                source=source,
                search_text=search_text,
                max_results=max_results,
                organism=organism,
                timeout_seconds=timeout_seconds,
            )
            for source in sources
        ]
        results = await asyncio.gather(*tasks)
        return {source: datasets for source, datasets in results}

    async def _search_one_source(
        self,
        *,
        source: DataSource,
        search_text: str,
        max_results: int,
        organism: str | None,
        timeout_seconds: float,
    ) -> tuple[DataSource, list[UnifiedDataset]]:
        adapter = self._adapters[source]
        try:
            datasets = await asyncio.wait_for(
                adapter.search(search_text, max_results=max_results, organism=organism),
                timeout=timeout_seconds,
            )
            return source, datasets
        except Exception as exc:
            logger.warning("Adapter search failed", source=source.value, error=str(exc))
            return source, []

    def _merge_results(self, results_by_source: dict[DataSource, list[UnifiedDataset]]) -> list[UnifiedDataset]:
        merged: list[UnifiedDataset] = []
        by_id: dict[str, int] = {}
        by_secondary: dict[str, int] = {}

        for source in sorted(results_by_source, key=lambda item: SOURCE_PRIORITY.get(item, 0), reverse=True):
            for dataset in results_by_source[source]:
                dataset_key = self._dataset_key(dataset)
                if dataset_key in by_id:
                    self._merge_dataset(merged[by_id[dataset_key]], dataset)
                    continue

                matched_index = self._find_secondary_match(dataset, by_secondary)
                if matched_index is not None:
                    primary = merged[matched_index]
                    if SOURCE_PRIORITY.get(self._as_source(dataset.source), 0) > SOURCE_PRIORITY.get(self._as_source(primary.source), 0):
                        dataset = self._merge_dataset(dataset, primary)
                        merged[matched_index] = dataset
                    else:
                        self._merge_dataset(primary, dataset)
                    continue

                by_id[dataset_key] = len(merged)
                for accession in self._all_accessions(dataset):
                    by_secondary.setdefault(accession, len(merged))
                merged.append(dataset)

        return merged

    def _select_sources(self, omics_types: list[OmicsType], requested: Iterable[DataSource] | None) -> list[DataSource]:
        if requested:
            selected = []
            for source in requested:
                coerced = self._as_source(source)
                if coerced in self._adapters:
                    selected.append(coerced)
            return selected

        selected = [source for source in get_sources_for_omics(omics_types) if source in self._adapters]
        if not selected and DataSource.GEO in self._adapters:
            return [DataSource.GEO]
        return selected

    def _coerce_query(self, query: UnifiedSearchQuery | str, overrides: dict[str, Any]) -> UnifiedSearchQuery:
        if isinstance(query, UnifiedSearchQuery):
            data = query.model_dump()
            data.update({key: value for key, value in overrides.items() if key in UnifiedSearchQuery.model_fields})
            return UnifiedSearchQuery(**data)

        fields = {
            key: value
            for key, value in overrides.items()
            if key in UnifiedSearchQuery.model_fields
        }
        return UnifiedSearchQuery(query=str(query), **fields)

    def _coerce_config(self, request: UnifiedSearchQuery, overrides: dict[str, Any]) -> SearchConfig:
        return SearchConfig(
            max_results=int(overrides.get("max_results") or request.max_results or self.config.max_results),
            min_relevance_score=float(
                overrides.get("min_relevance_score")
                if overrides.get("min_relevance_score") is not None
                else request.min_relevance_score
                if request.min_relevance_score is not None
                else self.config.min_relevance_score
            ),
            sort_by=str(overrides.get("sort_by") or request.sort_by or self.config.sort_by),
            source_timeout_seconds=float(
                overrides.get("source_timeout_seconds") or self.config.source_timeout_seconds
            ),
            concept_expansion=bool(
                overrides.get("concept_expansion")
                if "concept_expansion" in overrides
                else request.concept_expansion or self.config.concept_expansion
            ),
        )

    def _register_adapters(self) -> dict[DataSource, Any]:
        adapters: dict[DataSource, Any] = {}
        for source, factory in ADAPTER_REGISTRY:
            try:
                adapters[source] = factory()
            except Exception as exc:
                logger.warning("Adapter registration failed", source=source.value, error=str(exc))
        return adapters

    def _apply_sort(self, datasets: list[UnifiedDataset], sort_by: str) -> list[UnifiedDataset]:
        if sort_by not in VALID_SORT_OPTIONS:
            sort_by = "relevance"

        if sort_by == "date_desc":
            return sorted(datasets, key=lambda d: d.submission_date or "", reverse=True)
        if sort_by == "date_asc":
            return sorted(datasets, key=lambda d: d.submission_date or "")
        if sort_by == "sample_count":
            return sorted(datasets, key=lambda d: d.sample_count or 0, reverse=True)
        return sorted(datasets, key=lambda d: d.relevance_score or 0.0, reverse=True)

    def _merge_dataset(self, primary: UnifiedDataset, secondary: UnifiedDataset) -> UnifiedDataset:
        primary_source = self._as_source(primary.source)
        secondary_source = self._as_source(secondary.source)
        if primary_source and primary_source not in primary.merged_sources:
            primary.merged_sources.append(primary_source)
        if secondary_source and secondary_source not in primary.merged_sources:
            primary.merged_sources.append(secondary_source)

        primary.secondary_accession = sorted(set(primary.secondary_accession) | set(secondary.secondary_accession))
        primary.pubmed_ids = sorted(set(primary.pubmed_ids) | set(secondary.pubmed_ids))

        existing_urls = {link.url for link in primary.download_links}
        for link in secondary.download_links:
            if link.url not in existing_urls:
                primary.download_links.append(link)
                existing_urls.add(link.url)

        primary.extensions.update(secondary.extensions or {})
        primary.relevance_score = max(primary.relevance_score or 0.0, secondary.relevance_score or 0.0)

        if not primary.description and secondary.description:
            primary.description = secondary.description
        if primary.sample_count == 0 and secondary.sample_count:
            primary.sample_count = secondary.sample_count
        return primary

    def _dataset_key(self, dataset: UnifiedDataset) -> str:
        return f"{self._source_value(dataset.source)}:{dataset.accession}".lower()

    def _all_accessions(self, dataset: UnifiedDataset) -> set[str]:
        values = {dataset.accession}
        values.update(dataset.secondary_accession or [])
        return {value.lower().strip() for value in values if value}

    def _find_secondary_match(self, dataset: UnifiedDataset, by_secondary: dict[str, int]) -> int | None:
        for accession in self._all_accessions(dataset):
            if accession in by_secondary:
                return by_secondary[accession]
        return None

    def _summarize(self, datasets: list[UnifiedDataset]) -> tuple[dict[str, int], dict[str, int]]:
        source_counts: dict[str, int] = {}
        omics_counts: dict[str, int] = {}
        for dataset in datasets:
            source_key = self._source_value(dataset.source)
            omics_key = dataset.omics_type.value if hasattr(dataset.omics_type, "value") else str(dataset.omics_type)
            source_counts[source_key] = source_counts.get(source_key, 0) + 1
            omics_counts[omics_key] = omics_counts.get(omics_key, 0) + 1
        return source_counts, omics_counts

    def _as_source(self, value: DataSource | str | None) -> DataSource | None:
        if value is None:
            return None
        if isinstance(value, DataSource):
            return value
        try:
            return DataSource(str(value))
        except ValueError:
            return None

    def _source_value(self, value: DataSource | str) -> str:
        source = self._as_source(value)
        return source.value if source else str(value)


def get_core_search(config: SearchConfig | None = None) -> BioScouterCoreSearch:
    """Convenience factory for standalone public-core search."""

    return BioScouterCoreSearch(config=config)
