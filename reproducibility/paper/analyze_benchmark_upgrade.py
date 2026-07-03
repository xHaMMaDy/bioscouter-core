from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any


REQUIRED_INTERNAL = {
    "bioscouter_keyword",
    "bioscouter_embedding",
    "bioscouter_hybrid",
}
REQUIRED_BASELINES = {
    "native_source_api",
    "omicsdi",
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def numeric_row(row: dict[str, str]) -> dict[str, Any]:
    return {
        "system": row["system"],
        "queries": int(row["queries"]),
        "successful_queries": int(row["successful_queries"]),
        "failed_queries": int(row["failed_queries"]),
        "queries_with_results": int(row.get("queries_with_results") or row["successful_queries"]),
        "zero_result_queries": int(row.get("zero_result_queries") or 0),
        "returned_topk_records": int(row["returned_topk_records"]),
        "unique_records": int(row["unique_records"]),
        "observed_sources": int(row["observed_sources"]),
        "median_elapsed_s": float(row["median_elapsed_s"]),
        "mean_elapsed_s": float(row["mean_elapsed_s"]),
        "timing_scope": row.get("timing_scope") or "legacy_request_elapsed",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("controlled_run")
    parser.add_argument("--concept-run")
    args = parser.parse_args()

    controlled_dir = Path(args.controlled_run).resolve()
    rows = [
        numeric_row(row)
        for row in read_csv(controlled_dir / "system-summary.csv")
    ]
    by_system = {row["system"]: row for row in rows}
    hybrid = by_system.get("bioscouter_hybrid", {})
    systems = set(by_system)
    gate = {
        "at_least_60_hybrid_queries_with_results": hybrid.get("queries_with_results", 0) >= 60,
        "required_internal_modes_present": REQUIRED_INTERNAL <= systems,
        "required_baselines_present": REQUIRED_BASELINES <= systems,
        "traceable_outputs_present": all(
            (controlled_dir / filename).exists()
            for filename in (
                "system_runs.json",
                "raw_responses.json",
                "normalized-top10-results.csv",
                "system-summary.csv",
                "ablation-summary.csv",
                "baseline-comparison-summary.csv",
                "pool_mapping.csv",
                "run_manifest.json",
            )
        ),
    }

    concept: dict[str, Any] | None = None
    if args.concept_run:
        concept_dir = Path(args.concept_run).resolve()
        concept_rows = read_csv(concept_dir / "concept-ablation-aggregate.csv")
        concept = concept_rows[0] if concept_rows else None

    payload = {
        "controlled_run": str(controlled_dir),
        "systems": rows,
        "decision_gate": gate,
        "decision_gate_passed": all(gate.values()),
        "concept_ablation": concept,
        "relevance_boundary": (
            "The 80-query controlled outputs are not independently labeled. "
            "P@10, strict P@10, nDCG@10, MRR, and relevant-dataset superiority "
            "must not be inferred from returned counts or latency."
        ),
        "independent_relevance_evidence": (
            "Use the separately completed and adjudicated 30-query, 300-record evaluation."
        ),
    }
    json_path = controlled_dir / "benchmark-analysis.json"
    json_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    lines = [
        "# BioScouter Benchmark Upgrade Analysis",
        "",
        f"Controlled run: `{controlled_dir.name}`",
        "",
        "## Decision Gate",
        "",
    ]
    for key, passed in gate.items():
        lines.append(f"- {'PASS' if passed else 'FAIL'}: {key.replace('_', ' ')}")
    lines.extend(
        [
            "",
            f"Overall gate: **{'PASS' if all(gate.values()) else 'FAIL'}**",
            "",
            "## System Summary",
            "",
            "| System | Requests | With results | Returned top-10 | Unique | Sources |",
            "|---|---:|---:|---:|---:|---:|",
        ]
    )
    for row in rows:
        lines.append(
            f"| {row['system']} | {row['successful_queries']}/{row['queries']} | "
            f"{row['queries_with_results']}/{row['queries']} | "
            f"{row['returned_topk_records']} | {row['unique_records']} | "
            f"{row['observed_sources']} |"
        )
    if concept:
        lines.extend(
            [
                "",
                "## Concept-Normalization Ablation",
                "",
                f"- Paired queries: {concept['queries']}",
                f"- More/fewer/same returned counts: "
                f"{concept['queries_with_more_results']}/"
                f"{concept['queries_with_fewer_results']}/"
                f"{concept['queries_with_same_count']}",
                f"- Median candidate-set Jaccard overlap: {concept['median_jaccard_overlap']}",
            ]
        )
    lines.extend(
        [
            "",
            "## Interpretation Boundary",
            "",
            payload["relevance_boundary"],
            "",
            payload["independent_relevance_evidence"],
        ]
    )
    (controlled_dir / "benchmark-analysis.md").write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )
    print(json_path)
    return 0 if all(gate.values()) else 2


if __name__ == "__main__":
    raise SystemExit(main())
