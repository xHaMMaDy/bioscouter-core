"""Rerun the BioScouter manuscript benchmark query set.

This script records live API outputs for comparison against benchmark-freeze.json.
It does not overwrite the frozen manuscript snapshot.
"""

from __future__ import annotations

import argparse
import csv
import json
import time
from datetime import datetime, timezone
from pathlib import Path
from urllib import error, request


ROOT = Path(__file__).resolve().parent
QUERIES = ROOT / "benchmark-queries.csv"


def post_json(url: str, payload: dict, token: str | None) -> tuple[int, dict | str]:
    body = json.dumps(payload).encode("utf-8")
    headers = {"Content-Type": "application/json"}
    if token:
        headers["Authorization"] = f"Bearer {token}"

    req = request.Request(url, data=body, headers=headers, method="POST")
    try:
        with request.urlopen(req, timeout=120) as resp:
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


def result_count(payload: dict | str) -> int | None:
    if not isinstance(payload, dict):
        return None
    for key in ("results", "datasets"):
        value = payload.get(key)
        if isinstance(value, list):
            return len(value)
    if isinstance(payload.get("data"), dict):
        return result_count(payload["data"])
    return None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-url", default="http://127.0.0.1:8001")
    parser.add_argument("--token", default=None)
    parser.add_argument("--max-results", type=int, default=200)
    args = parser.parse_args()

    endpoint = args.base_url.rstrip("/") + "/api/omics/search"
    rows = []

    with QUERIES.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            if row["role"] != "latency_benchmark":
                continue

            payload = {"query": row["query"], "max_results": args.max_results}
            started = time.perf_counter()
            status, data = post_json(endpoint, payload, args.token)
            elapsed = time.perf_counter() - started
            rows.append({
                "id": row["id"],
                "query": row["query"],
                "status": status,
                "elapsed_s": round(elapsed, 3),
                "result_count": result_count(data),
                "response": data,
            })

    out = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "endpoint": endpoint,
        "max_results": args.max_results,
        "results": rows,
    }
    out_path = ROOT / f"live-rerun-{datetime.now(timezone.utc).strftime('%Y%m%d-%H%M%S')}.json"
    out_path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(out_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

