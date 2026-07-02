"""Run the expanded BioScouter benchmark panel.

This script records a dated live snapshot from the hosted/local API without
overwriting the original 10-query manuscript latency freeze.
"""

from __future__ import annotations

import argparse
import csv
import json
import statistics
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from urllib import error, request


ROOT = Path(__file__).resolve().parent
DEFAULT_QUERIES = ROOT / "expanded-benchmark-queries.csv"


def post_json(url: str, payload: dict, token: str | None) -> tuple[int, dict | str]:
    body = json.dumps(payload).encode("utf-8")
    headers = {"Content-Type": "application/json"}
    if token:
        headers["Authorization"] = f"Bearer {token}"

    req = request.Request(url, data=body, headers=headers, method="POST")
    try:
        with request.urlopen(req, timeout=180) as resp:
            text = resp.read().decode("utf-8")
            try:
                return resp.status, json.loads(text)
            except json.JSONDecodeError:
                return resp.status, text
    except error.HTTPError as exc:
        text = exc.read().decode("utf-8", errors="replace")
        return exc.code, text
    except error.URLError as exc:
        return 0, str(exc)


def as_list(value) -> str:
    if isinstance(value, list):
        return "; ".join(str(v) for v in value if v)
    if value is None:
        return ""
    return str(value)


def dataset_rows(query_row: dict, payload: dict | str, top_n: int) -> list[dict]:
    if not isinstance(payload, dict):
        return []
    datasets = payload.get("datasets") or payload.get("results") or []
    rows = []
    for rank, dataset in enumerate(datasets[:top_n], start=1):
        rows.append({
            "query_id": query_row["id"],
            "query": query_row["query"],
            "expected_omics_type": query_row["omics_type"],
            "rank": rank,
            "dataset_id": dataset.get("id", ""),
            "accession": dataset.get("accession", ""),
            "source": dataset.get("source", ""),
            "reported_omics_type": dataset.get("omics_type", ""),
            "title": dataset.get("title", ""),
            "description": (dataset.get("description") or "")[:1000].replace("\r", " ").replace("\n", " "),
            "organism": as_list(dataset.get("organism")),
            "sample_count": dataset.get("sample_count", ""),
            "disease": as_list(dataset.get("disease")),
            "tissue": as_list(dataset.get("tissue")),
            "relevance_score": dataset.get("relevance_score", ""),
            "quality_score": dataset.get("quality_score", ""),
            "source_url": dataset.get("source_url", ""),
            "annotator_label": "",
            "annotator_notes": "",
        })
    return rows


def summarize(results: list[dict]) -> dict:
    successful = [r for r in results if r["status"] == 200]
    counts = [r["total_results"] for r in successful if isinstance(r["total_results"], int)]
    elapsed = [r["elapsed_s"] for r in successful]
    source_counter: Counter[str] = Counter()
    omics_counter: Counter[str] = Counter()

    for row in successful:
        for source, count in (row.get("results_by_source") or {}).items():
            source_counter[str(source)] += int(count)
        for omics_type, count in (row.get("results_by_omics") or {}).items():
            omics_counter[str(omics_type)] += int(count)

    return {
        "queries": len(results),
        "successful_queries": len(successful),
        "failed_queries": len(results) - len(successful),
        "total_results": sum(counts),
        "median_results_per_query": statistics.median(counts) if counts else None,
        "mean_results_per_query": round(statistics.mean(counts), 2) if counts else None,
        "median_elapsed_s": round(statistics.median(elapsed), 3) if elapsed else None,
        "mean_elapsed_s": round(statistics.mean(elapsed), 3) if elapsed else None,
        "sources_observed": dict(source_counter.most_common()),
        "omics_observed": dict(omics_counter.most_common()),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-url", default="http://127.0.0.1:8001")
    parser.add_argument("--token", default=None)
    parser.add_argument("--queries", default=str(DEFAULT_QUERIES))
    parser.add_argument("--max-results", type=int, default=100)
    parser.add_argument("--top-n", type=int, default=10)
    parser.add_argument("--delay-s", type=float, default=1.0)
    args = parser.parse_args()

    endpoint = args.base_url.rstrip("/") + "/api/omics/search"
    created = datetime.now(timezone.utc)
    stamp = created.strftime("%Y%m%d-%H%M%S")
    out_dir = ROOT / f"expanded-benchmark-{stamp}"
    out_dir.mkdir(parents=True, exist_ok=True)

    query_path = Path(args.queries)
    results: list[dict] = []
    labeling_rows: list[dict] = []

    with query_path.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            if row["role"] != "expanded_coverage":
                continue

            payload = {
                "query": row["query"],
                "omics_types": [row["omics_type"]],
                "max_results": args.max_results,
                "sort_by": "relevance",
            }
            started = time.perf_counter()
            status, data = post_json(endpoint, payload, args.token)
            elapsed = round(time.perf_counter() - started, 3)

            total_results = None
            results_by_source = {}
            results_by_omics = {}
            if isinstance(data, dict):
                total_results = data.get("total_results")
                results_by_source = data.get("results_by_source") or {}
                results_by_omics = data.get("results_by_omics") or {}
                labeling_rows.extend(dataset_rows(row, data, args.top_n))

            results.append({
                "id": row["id"],
                "query": row["query"],
                "omics_type": row["omics_type"],
                "status": status,
                "elapsed_s": elapsed,
                "total_results": total_results,
                "results_by_source": results_by_source,
                "results_by_omics": results_by_omics,
                "response": data,
            })
            print(f"{row['id']} status={status} elapsed={elapsed}s results={total_results}")
            if args.delay_s:
                time.sleep(args.delay_s)

    snapshot = {
        "project": "BioScouter",
        "created_at": created.isoformat(),
        "endpoint": endpoint,
        "query_file": query_path.name,
        "max_results_per_source": args.max_results,
        "top_n_for_labeling": args.top_n,
        "snapshot_type": "expanded live coverage benchmark",
        "notes": [
            "This file is a dated live API snapshot and may differ from future reruns.",
            "It supplements, but does not replace, the original 10-query latency benchmark freeze.",
            "Top-10 labeling CSV files generated with this run should be independently annotated before final submission claims are updated.",
        ],
        "summary": summarize(results),
        "results": results,
    }

    snapshot_path = out_dir / "expanded-benchmark-freeze.json"
    snapshot_path.write_text(json.dumps(snapshot, indent=2), encoding="utf-8")

    summary_path = out_dir / "expanded-benchmark-summary.csv"
    with summary_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=[
            "id", "query", "omics_type", "status", "elapsed_s", "total_results",
            "results_by_source", "results_by_omics",
        ])
        writer.writeheader()
        for row in results:
            writer.writerow({
                **{k: row[k] for k in ["id", "query", "omics_type", "status", "elapsed_s", "total_results"]},
                "results_by_source": json.dumps(row["results_by_source"], sort_keys=True),
                "results_by_omics": json.dumps(row["results_by_omics"], sort_keys=True),
            })

    if labeling_rows:
        fields = list(labeling_rows[0].keys())
        for name in ("annotator_A_top10.csv", "annotator_B_top10.csv"):
            with (out_dir / name).open("w", newline="", encoding="utf-8") as fh:
                writer = csv.DictWriter(fh, fieldnames=fields)
                writer.writeheader()
                writer.writerows(labeling_rows)

    print(out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
