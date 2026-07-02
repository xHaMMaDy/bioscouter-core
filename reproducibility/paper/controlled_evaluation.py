"""Run the leakage-resistant BioScouter retrieval evaluation.

The controlled hybrid path does not use the production pgvector index. It
builds a frozen disease-agnostic corpus first, then evaluates held-out topic
queries against that corpus. Keyword, embedding, and hybrid rankings therefore
share traceable candidate inputs and cannot auto-index benchmark results.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import random
import re
import sys
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable
from urllib import error, parse, request


ROOT = Path(__file__).resolve().parent
CONTROLLED_QUERIES = ROOT / "controlled-benchmark-queries.csv"
CORPUS_QUERIES = ROOT / "corpus-construction-queries.csv"
SYSTEMS = (
    "bioscouter_keyword",
    "bioscouter_embedding",
    "bioscouter_hybrid",
    "native_source_api",
    "omicsdi",
)
REQUIRED_QUERY_FIELDS = {
    "query_id",
    "query",
    "omics_type",
    "native_source",
    "topic_family",
    "heldout_terms",
    "expected_constraints",
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False), encoding="utf-8")


def normalized_text(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", " ", value.casefold()).strip()


def validate_inputs(
    query_path: Path = CONTROLLED_QUERIES,
    corpus_path: Path = CORPUS_QUERIES,
) -> dict[str, Any]:
    queries = read_csv(query_path)
    corpus_queries = read_csv(corpus_path)
    errors: list[str] = []

    if len(queries) != 30:
        errors.append(f"Expected 30 controlled queries, found {len(queries)}")
    if len({row.get("query_id") for row in queries}) != len(queries):
        errors.append("Controlled query IDs must be unique")
    if len({row.get("topic_family") for row in queries}) != len(queries):
        errors.append("Controlled topic families must be unique")

    expected_omics = {
        "transcriptomics",
        "proteomics",
        "metabolomics",
        "epigenomics",
        "genomics",
        "single_cell",
        "metagenomics",
        "multi_omics",
    }
    observed_omics = {row.get("omics_type", "") for row in queries}
    if observed_omics != expected_omics:
        errors.append(
            f"Omics coverage mismatch: expected {sorted(expected_omics)}, "
            f"found {sorted(observed_omics)}"
        )

    corpus_text = normalized_text(" ".join(row.get("query", "") for row in corpus_queries))
    leaked_terms: list[str] = []
    for row in queries:
        missing = REQUIRED_QUERY_FIELDS - set(row)
        if missing:
            errors.append(f"{row.get('query_id', '?')} missing fields: {sorted(missing)}")
        for term in row.get("heldout_terms", "").split(";"):
            term = normalized_text(term)
            if term and term in corpus_text:
                leaked_terms.append(f"{row.get('query_id')}:{term}")
    if leaked_terms:
        errors.append(
            "Held-out topic terms appear in corpus-construction queries: "
            + ", ".join(leaked_terms)
        )

    result = {
        "valid": not errors,
        "controlled_queries": len(queries),
        "corpus_queries": len(corpus_queries),
        "omics_types": sorted(observed_omics),
        "errors": errors,
        "query_sha256": sha256_file(query_path),
        "corpus_query_sha256": sha256_file(corpus_path),
    }
    return result


def request_json(
    url: str,
    *,
    method: str = "GET",
    payload: dict[str, Any] | None = None,
    token: str | None = None,
    reviewer_token: str | None = None,
    timeout: float = 180,
) -> tuple[int, dict[str, Any] | str]:
    headers = {"Accept": "application/json", "User-Agent": "BioScouter-paper-evaluation/1.0"}
    data = None
    if payload is not None:
        headers["Content-Type"] = "application/json"
        data = json.dumps(payload).encode("utf-8")
    if token:
        headers["Authorization"] = f"Bearer {token}"
    if reviewer_token:
        headers["X-Reviewer-Token"] = reviewer_token
    req = request.Request(url, data=data, headers=headers, method=method)
    try:
        with request.urlopen(req, timeout=timeout) as response:
            text = response.read().decode("utf-8", errors="replace")
            try:
                return response.status, json.loads(text)
            except json.JSONDecodeError:
                return response.status, text
    except error.HTTPError as exc:
        text = exc.read().decode("utf-8", errors="replace")
        return exc.code, text
    except (error.URLError, TimeoutError) as exc:
        return 0, str(exc)


def as_text(value: Any) -> str:
    if isinstance(value, list):
        parts = []
        for item in value:
            if isinstance(item, dict):
                parts.append(str(item.get("name") or item.get("label") or item.get("value") or ""))
            elif item is not None:
                parts.append(str(item))
        return "; ".join(part for part in parts if part)
    if isinstance(value, dict):
        return str(value.get("name") or value.get("label") or value.get("value") or "")
    return "" if value is None else str(value)


def canonical_record(dataset: dict[str, Any], system: str) -> dict[str, Any]:
    source = as_text(dataset.get("source")).casefold()
    accession = as_text(dataset.get("accession") or dataset.get("id"))
    title = as_text(dataset.get("title"))
    key_seed = f"{source}|{accession}" if accession else f"{source}|{normalized_text(title)}"
    return {
        "record_key": key_seed,
        "system": system,
        "accession": accession,
        "source": source,
        "title": title,
        "description": as_text(dataset.get("description"))[:2000].replace("\r", " ").replace("\n", " "),
        "omics_type": as_text(dataset.get("omics_type") or dataset.get("omicsType")),
        "organism": as_text(dataset.get("organism") or dataset.get("organisms")),
        "sample_count": dataset.get("sample_count", ""),
        "disease": as_text(dataset.get("disease")),
        "tissue": as_text(dataset.get("tissue")),
        "source_url": as_text(dataset.get("source_url")),
        "relevance_score": dataset.get("relevance_score", ""),
    }


def bioscouter_search(
    endpoint: str,
    query_row: dict[str, str],
    *,
    max_results: int,
    sources: list[str] | None,
    token: str | None,
    reviewer_token: str | None,
) -> tuple[int, float, list[dict[str, Any]], dict[str, Any] | str]:
    payload: dict[str, Any] = {
        "query": query_row["query"],
        "omics_types": [query_row["omics_type"]],
        "max_results": max_results,
        "sort_by": "relevance",
        "ranking_mode": "keyword",
        "include_vector": False,
        "min_relevance_score": 0.0,
        "auto_index_results": False,
        "concept_expansion": False,
    }
    if sources:
        payload["sources"] = sources
    started = time.perf_counter()
    status, response = request_json(
        endpoint,
        method="POST",
        payload=payload,
        token=token,
        reviewer_token=reviewer_token,
    )
    elapsed = time.perf_counter() - started
    datasets = response.get("datasets", []) if isinstance(response, dict) else []
    return status, elapsed, datasets, response


def omicsdi_search(query: str, max_results: int) -> tuple[int, float, list[dict[str, Any]], Any]:
    url = "https://www.omicsdi.org/ws/dataset/search?" + parse.urlencode(
        {"query": query, "start": 0, "size": min(max_results, 100)}
    )
    started = time.perf_counter()
    status, response = request_json(url, timeout=90)
    elapsed = time.perf_counter() - started
    datasets = response.get("datasets", []) if isinstance(response, dict) else []
    normalized = []
    for dataset in datasets:
        item = dict(dataset)
        item["accession"] = item.get("id", "")
        item["omics_type"] = item.get("omicsType", [])
        item["organism"] = item.get("organisms", [])
        source = str(item.get("source", "")).casefold()
        accession = str(item.get("id", ""))
        item["source_url"] = (
            f"https://www.omicsdi.org/dataset/{parse.quote(source)}/{parse.quote(accession)}"
        )
        normalized.append(item)
    return status, elapsed, normalized, response


@dataclass
class EmbeddingRanker:
    model_name: str = "sentence-transformers/all-MiniLM-L6-v2"

    def __post_init__(self) -> None:
        try:
            import numpy as np
            from sentence_transformers import SentenceTransformer
        except ImportError as exc:
            raise RuntimeError(
                "Embedding evaluation requires numpy and sentence-transformers"
            ) from exc
        self.np = np
        self.model = SentenceTransformer(self.model_name)

    @staticmethod
    def document_text(record: dict[str, Any]) -> str:
        return " ".join(
            str(record.get(field, ""))
            for field in ("title", "description", "omics_type", "organism", "disease", "tissue")
            if record.get(field)
        )

    def encode(self, records: list[dict[str, Any]]) -> Any:
        texts = [self.document_text(record) for record in records]
        if not texts:
            return self.np.empty((0, 384), dtype="float32")
        return self.model.encode(
            texts,
            normalize_embeddings=True,
            show_progress_bar=False,
        )

    def rank(
        self,
        query: str,
        records: list[dict[str, Any]],
        *,
        top_k: int,
        precomputed_embeddings: Any | None = None,
    ) -> list[dict[str, Any]]:
        if not records:
            return []
        query_embedding = self.model.encode(
            [query],
            normalize_embeddings=True,
            show_progress_bar=False,
        )[0]
        embeddings = (
            precomputed_embeddings
            if precomputed_embeddings is not None
            else self.encode(records)
        )
        scores = embeddings @ query_embedding
        order = self.np.argsort(-scores)[:top_k]
        ranked = []
        for index in order:
            record = dict(records[int(index)])
            record["relevance_score"] = round(float(scores[int(index)]), 6)
            ranked.append(record)
        return ranked


def deduplicate(records: Iterable[dict[str, Any]]) -> list[dict[str, Any]]:
    selected: dict[str, dict[str, Any]] = {}
    for record in records:
        key = record["record_key"]
        current = selected.get(key)
        if current is None or len(record.get("description", "")) > len(current.get("description", "")):
            selected[key] = record
    return list(selected.values())


def build_corpus(
    *,
    endpoint: str,
    out_dir: Path,
    max_results: int,
    token: str | None,
    reviewer_token: str | None,
    delay_s: float,
) -> Path:
    corpus_records: list[dict[str, Any]] = []
    runs = []
    for row in read_csv(CORPUS_QUERIES):
        status, elapsed, datasets, raw = bioscouter_search(
            endpoint,
            row,
            max_results=max_results,
            sources=None,
            token=token,
            reviewer_token=reviewer_token,
        )
        records = [canonical_record(dataset, "corpus_builder") for dataset in datasets]
        corpus_records.extend(records)
        runs.append(
            {
                "query_id": row["query_id"],
                "query": row["query"],
                "status": status,
                "elapsed_s": round(elapsed, 3),
                "records": len(records),
                "response": raw,
            }
        )
        print(f"corpus {row['query_id']} status={status} records={len(records)}")
        if delay_s:
            time.sleep(delay_s)
    corpus_records = deduplicate(corpus_records)
    path = out_dir / "controlled-corpus.json"
    write_json(
        path,
        {
            "created_at": datetime.now(timezone.utc).isoformat(),
            "construction_query_sha256": sha256_file(CORPUS_QUERIES),
            "auto_index_results": False,
            "record_count": len(corpus_records),
            "records": corpus_records,
            "runs": runs,
        },
    )
    return path


def pool_results(
    queries: list[dict[str, str]],
    run_rows: list[dict[str, Any]],
    out_dir: Path,
    *,
    top_k: int,
) -> None:
    query_lookup = {row["query_id"]: row for row in queries}
    pooled: dict[tuple[str, str], dict[str, Any]] = {}
    mapping_rows: list[dict[str, Any]] = []

    for run in run_rows:
        query_id = run["query_id"]
        for rank, record in enumerate(run["records"][:top_k], start=1):
            pool_seed = f"{query_id}|{record['record_key']}"
            pool_id = hashlib.sha256(pool_seed.encode("utf-8")).hexdigest()[:16]
            key = (query_id, pool_id)
            pooled.setdefault(key, {**record, "query_id": query_id, "pool_id": pool_id})
            mapping_rows.append(
                {
                    "query_id": query_id,
                    "pool_id": pool_id,
                    "system": run["system"],
                    "rank": rank,
                    "record_key": record["record_key"],
                }
            )

    mapping_path = out_dir / "pool_mapping.csv"
    with mapping_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["query_id", "pool_id", "system", "rank", "record_key"],
        )
        writer.writeheader()
        writer.writerows(mapping_rows)

    fields = [
        "query_id",
        "query",
        "expected_omics_type",
        "expected_constraints",
        "pool_id",
        "accession",
        "source",
        "reported_omics_type",
        "title",
        "description",
        "organism",
        "sample_count",
        "disease",
        "tissue",
        "source_url",
        "annotator_label",
        "annotator_notes",
    ]
    base_rows = []
    for (query_id, _pool_id), record in pooled.items():
        query = query_lookup[query_id]
        base_rows.append(
            {
                "query_id": query_id,
                "query": query["query"],
                "expected_omics_type": query["omics_type"],
                "expected_constraints": query["expected_constraints"],
                "pool_id": record["pool_id"],
                "accession": record.get("accession", ""),
                "source": record.get("source", ""),
                "reported_omics_type": record.get("omics_type", ""),
                "title": record.get("title", ""),
                "description": record.get("description", ""),
                "organism": record.get("organism", ""),
                "sample_count": record.get("sample_count", ""),
                "disease": record.get("disease", ""),
                "tissue": record.get("tissue", ""),
                "source_url": record.get("source_url", ""),
                "annotator_label": "",
                "annotator_notes": "",
            }
        )

    for filename, seed in (("annotator_A_blinded.csv", 1729), ("annotator_B_blinded.csv", 2718)):
        rows = list(base_rows)
        random.Random(seed).shuffle(rows)
        with (out_dir / filename).open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)


def run_controlled(args: argparse.Namespace) -> int:
    validation = validate_inputs(Path(args.queries), Path(args.corpus_queries))
    if not validation["valid"]:
        print(json.dumps(validation, indent=2))
        return 2

    created = datetime.now(timezone.utc)
    stamp = created.strftime("%Y%m%d-%H%M%S")
    out_dir = Path(args.out_dir) if args.out_dir else ROOT / f"controlled-evaluation-{stamp}"
    out_dir.mkdir(parents=True, exist_ok=True)
    endpoint = args.base_url.rstrip("/") + "/api/omics/search"
    corpus_path = Path(args.corpus) if args.corpus else out_dir / "controlled-corpus.json"
    if not corpus_path.exists():
        corpus_path = build_corpus(
            endpoint=endpoint,
            out_dir=out_dir,
            max_results=args.corpus_max_results,
            token=args.token,
            reviewer_token=args.reviewer_token,
            delay_s=args.delay_s,
        )

    corpus_payload = json.loads(corpus_path.read_text(encoding="utf-8"))
    corpus_records = corpus_payload.get("records", [])
    ranker = EmbeddingRanker(args.embedding_model)
    corpus_embeddings = ranker.encode(corpus_records)
    queries = read_csv(Path(args.queries))
    run_rows: list[dict[str, Any]] = []
    raw_responses: list[dict[str, Any]] = []

    for row in queries:
        status, elapsed, datasets, raw = bioscouter_search(
            endpoint,
            row,
            max_results=args.candidate_depth,
            sources=None,
            token=args.token,
            reviewer_token=args.reviewer_token,
        )
        keyword_records = [
            canonical_record(dataset, "bioscouter_keyword") for dataset in datasets
        ]
        keyword_records = deduplicate(keyword_records)
        run_rows.append(
            {
                "query_id": row["query_id"],
                "system": "bioscouter_keyword",
                "status": status,
                "elapsed_s": round(elapsed, 3),
                "records": keyword_records[: args.top_k],
            }
        )

        embedding_records = ranker.rank(
            row["query"],
            keyword_records,
            top_k=args.top_k,
        )
        for record in embedding_records:
            record["system"] = "bioscouter_embedding"
        run_rows.append(
            {
                "query_id": row["query_id"],
                "system": "bioscouter_embedding",
                "status": status,
                "elapsed_s": round(elapsed, 3),
                "records": embedding_records,
            }
        )

        semantic_records = ranker.rank(
            row["query"],
            corpus_records,
            top_k=args.semantic_depth,
            precomputed_embeddings=corpus_embeddings,
        )
        hybrid_candidates = deduplicate(keyword_records + semantic_records)
        hybrid_records = ranker.rank(
            row["query"],
            hybrid_candidates,
            top_k=args.top_k,
        )
        for record in hybrid_records:
            record["system"] = "bioscouter_hybrid"
        run_rows.append(
            {
                "query_id": row["query_id"],
                "system": "bioscouter_hybrid",
                "status": status,
                "elapsed_s": round(elapsed, 3),
                "records": hybrid_records,
            }
        )

        native_status, native_elapsed, native_datasets, native_raw = bioscouter_search(
            endpoint,
            row,
            max_results=args.top_k,
            sources=[row["native_source"]],
            token=args.token,
            reviewer_token=args.reviewer_token,
        )
        native_records = [
            canonical_record(dataset, "native_source_api")
            for dataset in native_datasets[: args.top_k]
        ]
        run_rows.append(
            {
                "query_id": row["query_id"],
                "system": "native_source_api",
                "status": native_status,
                "elapsed_s": round(native_elapsed, 3),
                "records": native_records,
            }
        )

        omicsdi_status, omicsdi_elapsed, omicsdi_datasets, omicsdi_raw = omicsdi_search(
            row["query"],
            args.top_k,
        )
        omicsdi_records = [
            canonical_record(dataset, "omicsdi")
            for dataset in omicsdi_datasets[: args.top_k]
        ]
        run_rows.append(
            {
                "query_id": row["query_id"],
                "system": "omicsdi",
                "status": omicsdi_status,
                "elapsed_s": round(omicsdi_elapsed, 3),
                "records": omicsdi_records,
            }
        )
        raw_responses.append(
            {
                "query_id": row["query_id"],
                "query": row["query"],
                "bioscouter": raw,
                "native_source_api": native_raw,
                "omicsdi": omicsdi_raw,
            }
        )
        print(
            f"{row['query_id']} BioScouter={len(keyword_records)} "
            f"native={len(native_records)} OmicsDI={len(omicsdi_records)}"
        )
        if args.delay_s:
            time.sleep(args.delay_s)

    run_path = out_dir / "system_runs.json"
    write_json(
        run_path,
        {
            "created_at": created.isoformat(),
            "systems": SYSTEMS,
            "top_k": args.top_k,
            "candidate_depth": args.candidate_depth,
            "semantic_depth": args.semantic_depth,
            "embedding_model": args.embedding_model,
            "production_vector_index_used": False,
            "auto_index_results": False,
            "corpus_sha256": sha256_file(corpus_path),
            "runs": run_rows,
        },
    )
    write_json(out_dir / "raw_responses.json", raw_responses)
    pool_results(queries, run_rows, out_dir, top_k=args.top_k)

    manifest_files = [
        Path(args.queries),
        Path(args.corpus_queries),
        corpus_path,
        run_path,
        out_dir / "raw_responses.json",
        out_dir / "pool_mapping.csv",
        out_dir / "annotator_A_blinded.csv",
        out_dir / "annotator_B_blinded.csv",
    ]
    manifest = {
        "created_at": created.isoformat(),
        "validation": validation,
        "files": {
            path.name: {"sha256": sha256_file(path), "bytes": path.stat().st_size}
            for path in manifest_files
        },
    }
    write_json(out_dir / "run_manifest.json", manifest)
    print(out_dir)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--base-url", default="http://127.0.0.1:8001")
    parser.add_argument("--token")
    parser.add_argument("--reviewer-token")
    parser.add_argument("--queries", default=str(CONTROLLED_QUERIES))
    parser.add_argument("--corpus-queries", default=str(CORPUS_QUERIES))
    parser.add_argument("--corpus")
    parser.add_argument("--out-dir")
    parser.add_argument("--top-k", type=int, default=10)
    parser.add_argument("--candidate-depth", type=int, default=100)
    parser.add_argument("--semantic-depth", type=int, default=50)
    parser.add_argument("--corpus-max-results", type=int, default=100)
    parser.add_argument("--delay-s", type=float, default=0.5)
    parser.add_argument(
        "--embedding-model",
        default="sentence-transformers/all-MiniLM-L6-v2",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    validation = validate_inputs(Path(args.queries), Path(args.corpus_queries))
    if args.validate_only:
        print(json.dumps(validation, indent=2))
        return 0 if validation["valid"] else 2
    return run_controlled(args)


if __name__ == "__main__":
    sys.exit(main())
