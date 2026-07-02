"""Command-line entry point for public-core searches."""

from __future__ import annotations

import argparse
import asyncio
import json
from typing import Any

from .models.unified import DataSource, UnifiedSearchQuery
from .orchestrator import BioScouterCoreSearch, SearchConfig


def _parse_sources(values: list[str] | None) -> list[DataSource] | None:
    if not values:
        return None
    return [DataSource(value.lower()) for value in values]


def _to_jsonable(value: Any) -> Any:
    if hasattr(value, "model_dump"):
        return value.model_dump(mode="json")
    return value


async def _run(args: argparse.Namespace) -> int:
    searcher = BioScouterCoreSearch(
        config=SearchConfig(
            max_results=args.max_results,
            sort_by=args.sort_by,
            source_timeout_seconds=args.timeout,
            concept_expansion=args.concept_expansion,
        )
    )
    try:
        response = await searcher.search(
            UnifiedSearchQuery(
                query=args.query,
                sources=_parse_sources(args.source),
                organism=args.organism,
                max_results=args.max_results,
                sort_by=args.sort_by,
                concept_expansion=args.concept_expansion,
            )
        )
    finally:
        await searcher.close()

    if args.json:
        print(json.dumps(_to_jsonable(response), indent=2, ensure_ascii=False))
        return 0

    print(f"{response.total_results} datasets from {len(response.sources_searched)} sources")
    for dataset in response.datasets:
        print(f"{dataset.accession}\t{dataset.source}\t{dataset.relevance_score:.2f}\t{dataset.title}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description="Search public omics repositories through BioScouter core.")
    parser.add_argument("query", help="Natural-language or keyword query")
    parser.add_argument("--source", action="append", help="Limit to a source, e.g. geo, pride, cptac")
    parser.add_argument("--organism", help="Organism filter, e.g. Homo sapiens")
    parser.add_argument("--max-results", type=int, default=20)
    parser.add_argument("--sort-by", default="relevance", choices=["relevance", "date_desc", "date_asc", "sample_count"])
    parser.add_argument("--timeout", type=float, default=30.0, help="Per-source timeout in seconds")
    parser.add_argument("--concept-expansion", action="store_true", help="Use deterministic curated concept expansion")
    parser.add_argument("--json", action="store_true", help="Print full JSON response")
    return asyncio.run(_run(parser.parse_args()))


if __name__ == "__main__":
    raise SystemExit(main())

