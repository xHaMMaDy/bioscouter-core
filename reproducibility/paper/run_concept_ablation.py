"""Run paired searches with and without deterministic concept expansion.

The output is descriptive evidence only. It records candidate overlap and count
changes so the manuscript can discuss concept normalization without claiming
relevance gains before independent labels are complete.
"""

from __future__ import annotations

import argparse
import csv
import json
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from urllib import error, request


ROOT = Path(__file__).resolve().parent
CONTROLLED_QUERIES = ROOT / "controlled-benchmark-queries.csv"


def read_queries(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def request_json(
    url: str,
    payload: dict[str, Any],
    *,
    token: str | None,
    reviewer_token: str | None,
    timeout: float = 180,
) -> tuple[int, dict[str, Any] | str]:
    headers = {
        "Accept": "application/json",
        "Content-Type": "application/json",
        "User-Agent": "BioScouter-concept-ablation/1.0",
    }
    if token:
        headers["Authorization"] = f"Bearer {token}"
    if reviewer_token:
        headers["X-Reviewer-Token"] = reviewer_token

    req = request.Request(
        url,
        data=json.dumps(payload).encode("utf-8"),
        headers=headers,
        method="POST",
    )
    try:
        with request.urlopen(req, timeout=timeout) as response:
            text = response.read().decode("utf-8", errors="replace")
            try:
                return response.status, json.loads(text)
            except json.JSONDecodeError:
                return response.status, text
    except error.HTTPError as exc:
        return exc.code, exc.read().decode("utf-8", errors="replace")
    except (TimeoutError, error.URLError) as exc:
        return 0, str(exc)


def result_key(dataset: dict[str, Any]) -> str:
    return f"{dataset.get('source', '')}|{dataset.get('accession') or dataset.get('id') or dataset.get('title', '')}"


def run_variant(
    endpoint: str,
    row: dict[str, str],
    *,
    concept_expansion: bool,
    max_results: int,
    token: str | None,
    reviewer_token: str | None,
) -> tuple[int, float, dict[str, Any] | str]:
    payload = {
        "query": row["query"],
        "omics_types": [row["omics_type"]],
        "max_results": max_results,
        "sort_by": "relevance",
        "ranking_mode": "keyword",
        "include_vector": False,
        "min_relevance_score": 0.0,
        "auto_index_results": False,
        "concept_expansion": concept_expansion,
    }
    started = time.perf_counter()
    status, data = request_json(endpoint, payload, token=token, reviewer_token=reviewer_token)
    return status, time.perf_counter() - started, data


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-url", default="http://127.0.0.1:8001")
    parser.add_argument("--queries", default=str(CONTROLLED_QUERIES))
    parser.add_argument("--max-results", type=int, default=100)
    parser.add_argument("--token")
    parser.add_argument("--reviewer-token")
    parser.add_argument("--out-dir")
    parser.add_argument("--delay-s", type=float, default=0.5)
    args = parser.parse_args()

    created = datetime.now(timezone.utc)
    out_dir = Path(args.out_dir) if args.out_dir else ROOT / f"concept-ablation-{created.strftime('%Y%m%d-%H%M%S')}"
    out_dir.mkdir(parents=True, exist_ok=True)
    endpoint = args.base_url.rstrip("/") + "/api/omics/search"

    records: list[dict[str, Any]] = []
    for row in read_queries(Path(args.queries)):
        status_plain, elapsed_plain, plain = run_variant(
            endpoint,
            row,
            concept_expansion=False,
            max_results=args.max_results,
            token=args.token,
            reviewer_token=args.reviewer_token,
        )
        if args.delay_s:
            time.sleep(args.delay_s)
        status_expanded, elapsed_expanded, expanded = run_variant(
            endpoint,
            row,
            concept_expansion=True,
            max_results=args.max_results,
            token=args.token,
            reviewer_token=args.reviewer_token,
        )

        plain_datasets = plain.get("datasets", []) if isinstance(plain, dict) else []
        expanded_datasets = expanded.get("datasets", []) if isinstance(expanded, dict) else []
        plain_keys = {result_key(dataset) for dataset in plain_datasets}
        expanded_keys = {result_key(dataset) for dataset in expanded_datasets}
        shared = plain_keys & expanded_keys
        union = plain_keys | expanded_keys

        records.append(
            {
                "query_id": row["query_id"],
                "query": row["query"],
                "omics_type": row["omics_type"],
                "plain_status": status_plain,
                "expanded_status": status_expanded,
                "plain_results": len(plain_datasets),
                "expanded_results": len(expanded_datasets),
                "shared_results": len(shared),
                "jaccard_overlap": round(len(shared) / len(union), 4) if union else 1.0,
                "plain_elapsed_s": round(elapsed_plain, 3),
                "expanded_elapsed_s": round(elapsed_expanded, 3),
            }
        )
        print(
            f"{row['query_id']} plain={len(plain_datasets)} "
            f"expanded={len(expanded_datasets)} shared={len(shared)}"
        )
        if args.delay_s:
            time.sleep(args.delay_s)

    if records:
        csv_path = out_dir / "concept-ablation-summary.csv"
        with csv_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(records[0].keys()))
            writer.writeheader()
            writer.writerows(records)

    (out_dir / "concept-ablation-run.json").write_text(
        json.dumps(
            {
                "created_at": created.isoformat(),
                "interpretation": "Descriptive candidate-set ablation; not an independent relevance metric.",
                "endpoint": endpoint,
                "max_results": args.max_results,
                "queries": str(Path(args.queries)),
                "records": records,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    print(out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
