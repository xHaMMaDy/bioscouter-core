"""Measure per-source BioScouter adapter reliability without changing the index."""

from __future__ import annotations

import argparse
import csv
import json
import platform
import statistics
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from controlled_evaluation import request_json, sha256_file


ROOT = Path(__file__).resolve().parent
DEFAULT_PROBES = ROOT / "source-reliability-probes.csv"


def percentile(values: list[float], fraction: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    index = round((len(ordered) - 1) * fraction)
    return ordered[index]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-url", default="http://127.0.0.1:8001")
    parser.add_argument("--token")
    parser.add_argument("--reviewer-token")
    parser.add_argument("--probes", default=str(DEFAULT_PROBES))
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--max-results", type=int, default=20)
    parser.add_argument("--delay-s", type=float, default=0.5)
    parser.add_argument("--out-dir")
    args = parser.parse_args()

    probe_path = Path(args.probes)
    with probe_path.open(newline="", encoding="utf-8") as handle:
        probes = list(csv.DictReader(handle))
    endpoint = args.base_url.rstrip("/") + "/api/omics/search"
    created = datetime.now(timezone.utc)
    out_dir = (
        Path(args.out_dir)
        if args.out_dir
        else ROOT / f"source-reliability-{created.strftime('%Y%m%d-%H%M%S')}"
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    runs = []

    for probe in probes:
        for repetition in range(1, args.repetitions + 1):
            payload = {
                "query": probe["query"],
                "omics_types": [probe["omics_type"]],
                "sources": [probe["source"]],
                "max_results": args.max_results,
                "ranking_mode": "keyword",
                "include_vector": False,
                "min_relevance_score": 0.0,
                "auto_index_results": False,
            }
            started = time.perf_counter()
            status, response = request_json(
                endpoint,
                method="POST",
                payload=payload,
                token=args.token,
                reviewer_token=args.reviewer_token,
            )
            elapsed = round(time.perf_counter() - started, 3)
            datasets = response.get("datasets", []) if isinstance(response, dict) else []
            runs.append(
                {
                    "source": probe["source"],
                    "query": probe["query"],
                    "omics_type": probe["omics_type"],
                    "repetition": repetition,
                    "cache_state": "cold_observation" if repetition == 1 else "warm_observation",
                    "status": status,
                    "success": status == 200,
                    "elapsed_s": elapsed,
                    "result_count": len(datasets),
                    "error": "" if status == 200 else str(response)[:500],
                }
            )
            print(
                f"{probe['source']} run={repetition} status={status} "
                f"elapsed={elapsed}s results={len(datasets)}"
            )
            if args.delay_s:
                time.sleep(args.delay_s)

    summary = []
    for source in [probe["source"] for probe in probes]:
        source_runs = [row for row in runs if row["source"] == source]
        successful = [row for row in source_runs if row["success"]]
        elapsed = [row["elapsed_s"] for row in successful]
        summary.append(
            {
                "source": source,
                "attempts": len(source_runs),
                "successes": len(successful),
                "success_rate": round(len(successful) / len(source_runs), 3),
                "median_elapsed_s": round(statistics.median(elapsed), 3) if elapsed else "",
                "p95_elapsed_s": round(percentile(elapsed, 0.95), 3) if elapsed else "",
                "median_results": (
                    statistics.median(row["result_count"] for row in successful)
                    if successful
                    else ""
                ),
            }
        )

    run_fields = list(runs[0].keys())
    with (out_dir / "source-reliability-runs.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=run_fields)
        writer.writeheader()
        writer.writerows(runs)
    with (out_dir / "source-reliability-summary.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0].keys()))
        writer.writeheader()
        writer.writerows(summary)

    manifest = {
        "created_at": created.isoformat(),
        "endpoint": endpoint,
        "probe_sha256": sha256_file(probe_path),
        "repetitions": args.repetitions,
        "max_results": args.max_results,
        "python": sys.version,
        "platform": platform.platform(),
        "interpretation": (
            "The first request is labeled a cold observation and later requests warm "
            "observations; upstream and server caches are not forcibly purged."
        ),
        "summary": summary,
    }
    (out_dir / "source-reliability-manifest.json").write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8",
    )
    print(out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
