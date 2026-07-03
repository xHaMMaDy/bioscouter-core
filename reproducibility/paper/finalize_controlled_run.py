from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent
DEFAULT_QUERIES = ROOT / "benchmark-queries-80.csv"
DEFAULT_CORPUS_QUERIES = ROOT / "corpus-construction-queries.csv"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def enrich_summary_tables(out_dir: Path, query_count: int) -> None:
    normalized_path = out_dir / "normalized-top10-results.csv"
    if not normalized_path.exists():
        return

    queries_with_results: dict[str, set[str]] = {}
    with normalized_path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            system = row.get("system", "")
            query_id = row.get("query_id", "")
            if system and query_id:
                queries_with_results.setdefault(system, set()).add(query_id)

    for name in ("system-summary.csv", "ablation-summary.csv", "baseline-comparison-summary.csv"):
        path = out_dir / name
        if not path.exists():
            continue
        with path.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        if not rows:
            continue
        fieldnames = list(rows[0].keys())
        if "queries_with_results" not in fieldnames:
            insert_after = fieldnames.index("failed_queries") + 1 if "failed_queries" in fieldnames else len(fieldnames)
            fieldnames[insert_after:insert_after] = ["queries_with_results", "zero_result_queries"]
        if "timing_scope" not in fieldnames:
            fieldnames.append("timing_scope")

        for row in rows:
            system = row.get("system", "")
            row_count = len(queries_with_results.get(system, set()))
            total = int(row.get("queries") or query_count)
            row["queries_with_results"] = str(row_count)
            row["zero_result_queries"] = str(max(0, total - row_count))
            row["timing_scope"] = row.get("timing_scope") or "legacy_source_request_elapsed_only"

        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("out_dir")
    parser.add_argument("--queries", default=str(DEFAULT_QUERIES))
    parser.add_argument("--corpus-queries", default=str(DEFAULT_CORPUS_QUERIES))
    parser.add_argument("--base-url", default="http://127.0.0.1:8001")
    args = parser.parse_args()

    out_dir = Path(args.out_dir).resolve()
    queries = Path(args.queries).resolve()
    corpus_queries = Path(args.corpus_queries).resolve()
    run_path = out_dir / "system_runs.json"
    corpus_path = out_dir / "controlled-corpus.json"
    manifest_path = out_dir / "run_manifest.json"
    if not run_path.exists() or not corpus_path.exists():
        raise FileNotFoundError("The controlled run is incomplete: system_runs.json or controlled-corpus.json is missing.")

    run = read_json(run_path)
    manifest = read_json(manifest_path) if manifest_path.exists() else {}
    query_count = len({row["query_id"] for row in run.get("runs", [])})
    enrich_summary_tables(out_dir, query_count)
    timing_scopes = sorted(
        {
            str(row.get("timing_scope", "")).strip()
            for row in run.get("runs", [])
            if str(row.get("timing_scope", "")).strip()
        }
    )
    config = {
        "created_at": run.get("created_at"),
        "base_url": args.base_url,
        "query_file": str(queries),
        "query_count": query_count,
        "query_sha256": sha256_file(queries),
        "corpus_query_file": str(corpus_queries),
        "corpus_query_sha256": sha256_file(corpus_queries),
        "corpus_file": str(corpus_path),
        "corpus_sha256": sha256_file(corpus_path),
        "top_k": run.get("top_k"),
        "candidate_depth": run.get("candidate_depth"),
        "semantic_depth": run.get("semantic_depth"),
        "embedding_model": run.get("embedding_model"),
        "systems": run.get("systems", []),
        "production_vector_index_used": run.get("production_vector_index_used", False),
        "auto_index_results": run.get("auto_index_results", False),
        "concept_expansion": False,
        "timing_scope": (
            timing_scopes
            if timing_scopes
            else ["legacy_source_request_elapsed_only; excludes local semantic reranking"]
        ),
    }
    config_path = out_dir / "benchmark_run_config.json"
    write_json(config_path, config)

    tracked = [queries, corpus_queries]
    tracked.extend(
        sorted(
            path
            for path in out_dir.iterdir()
            if path.is_file() and path.name != manifest_path.name
        )
    )
    manifest["created_at"] = manifest.get("created_at") or run.get("created_at")
    manifest["files"] = {
        path.name: {"sha256": sha256_file(path), "bytes": path.stat().st_size}
        for path in tracked
    }
    write_json(manifest_path, manifest)
    print(out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
