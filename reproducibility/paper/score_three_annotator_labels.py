"""Score the 80-query three-annotator relevance evaluation.

The evaluator packet contains three complete blinded label sheets plus a
separate adjudication sheet for the 0/1/2 split rows. This script validates the
packet, builds final consensus labels, and writes P@10, strict P@10, nDCG@10,
MRR, agreement, and paired comparison outputs.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
import sys
from collections import Counter
from pathlib import Path
from typing import Any

from score_relevance_labels import (
    SYSTEM_ORDER,
    VALID_LABELS,
    agreement_summary,
    pairwise_comparisons,
    per_query_metrics,
    read_mapping,
    summarize_systems,
    write_csv,
    write_json,
)


def default_packet_dir() -> Path:
    script_path = Path(__file__).resolve()
    candidates = [
        script_path.parents[1]
        / "evaluator_packet"
        / "BioScouter_80Q_Three_Annotator_Evaluation_20260703",
        script_path.parent / "independent_relevance_evaluation_80q_20260703",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


DEFAULT_PACKET_DIR = default_packet_dir()

ANNOTATOR_FILES = {
    "A": "BioScouter_80Q_Annotator_A_Labeling_Sheet_COMPLETED.csv",
    "B": "BioScouter_80Q_Annotator_B_Labeling_Sheet_COMPLETED.csv",
    "C": "BioScouter_80Q_Annotator_C_Labeling_Sheet_COMPLETED.csv",
}

IDENTITY_FIELDS = (
    "query_id",
    "pool_id",
    "accession",
    "source",
    "reported_omics_type",
    "title",
    "source_url",
)


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def read_annotator_sheet(
    path: Path,
    annotator: str,
) -> tuple[dict[tuple[str, str], int], dict[tuple[str, str], dict[str, str]], list[str]]:
    labels: dict[tuple[str, str], int] = {}
    rows_by_key: dict[tuple[str, str], dict[str, str]] = {}
    missing: list[str] = []
    for row_number, row in enumerate(read_csv_rows(path), start=2):
        query_id = row.get("query_id", "").strip()
        pool_id = row.get("pool_id", "").strip()
        if not query_id or not pool_id:
            raise ValueError(f"{path}:{row_number} is missing query_id or pool_id")
        key = (query_id, pool_id)
        if key in labels:
            raise ValueError(f"{path}:{row_number} duplicates {query_id}/{pool_id}")
        raw_label = row.get("annotator_label", "").strip()
        if raw_label not in VALID_LABELS:
            if not raw_label:
                missing.append(f"{query_id}/{pool_id}")
                continue
            raise ValueError(
                f"{path}:{row_number} has invalid annotator_label={raw_label!r}"
            )
        labels[key] = int(raw_label)
        rows_by_key[key] = row
    if missing:
        raise ValueError(f"Annotator {annotator} has {len(missing)} missing labels")
    return labels, rows_by_key, missing


def read_adjudication(path: Path) -> dict[tuple[str, str], dict[str, str]]:
    rows_by_key: dict[tuple[str, str], dict[str, str]] = {}
    for row_number, row in enumerate(read_csv_rows(path), start=2):
        query_id = row.get("query_id", "").strip()
        pool_id = row.get("pool_id", "").strip()
        raw_label = row.get("adjudicated_label", "").strip()
        if not query_id or not pool_id:
            raise ValueError(f"{path}:{row_number} is missing query_id or pool_id")
        if raw_label not in VALID_LABELS:
            raise ValueError(
                f"{path}:{row_number} has invalid adjudicated_label={raw_label!r}"
            )
        key = (query_id, pool_id)
        if key in rows_by_key:
            raise ValueError(f"{path}:{row_number} duplicates {query_id}/{pool_id}")
        rows_by_key[key] = row
    return rows_by_key


def majority_label(labels: list[int]) -> int | None:
    counts = Counter(labels)
    label, count = counts.most_common(1)[0]
    return int(label) if count >= 2 else None


def fleiss_kappa(label_sets: dict[str, dict[tuple[str, str], int]]) -> float:
    keys = sorted(next(iter(label_sets.values())).keys())
    n_raters = len(label_sets)
    categories = [0, 1, 2]
    row_agreements = []
    category_totals = Counter()
    for key in keys:
        counts = Counter(labels[key] for labels in label_sets.values())
        category_totals.update(counts)
        row_agreements.append(
            sum(count * (count - 1) for count in counts.values())
            / (n_raters * (n_raters - 1))
        )
    p_bar = statistics.mean(row_agreements)
    total_ratings = len(keys) * n_raters
    p_e = sum((category_totals[label] / total_ratings) ** 2 for label in categories)
    if math.isclose(1.0, p_e):
        return 1.0 if math.isclose(p_bar, 1.0) else 0.0
    return (p_bar - p_e) / (1 - p_e)


def validate_same_row_sets(
    label_sets: dict[str, dict[tuple[str, str], int]],
    row_sets: dict[str, dict[tuple[str, str], dict[str, str]]],
) -> list[str]:
    keys_by_annotator = {name: set(labels) for name, labels in label_sets.items()}
    reference = keys_by_annotator["A"]
    errors: list[str] = []
    for name, keys in keys_by_annotator.items():
        if keys != reference:
            errors.append(
                f"Annotator {name} row set mismatch: "
                f"missing={len(reference - keys)}, extra={len(keys - reference)}"
            )
    for key in sorted(reference):
        row_a = row_sets["A"][key]
        for name in ("B", "C"):
            row = row_sets[name][key]
            for field in IDENTITY_FIELDS:
                if (row.get(field, "") or "") != (row_a.get(field, "") or ""):
                    errors.append(
                        f"{key[0]}/{key[1]} metadata mismatch in {field} for {name}"
                    )
                    break
    return errors


def build_final_consensus(
    label_sets: dict[str, dict[tuple[str, str], int]],
    row_sets: dict[str, dict[tuple[str, str], dict[str, str]]],
    adjudication: dict[tuple[str, str], dict[str, str]],
) -> tuple[dict[tuple[str, str], int], list[dict[str, Any]], dict[str, Any]]:
    labels: dict[tuple[str, str], int] = {}
    rows: list[dict[str, Any]] = []
    majority_rows = 0
    unanimous_rows = 0
    adjudicated_rows = 0
    unexpected_adjudication_rows = sorted(set(adjudication) - set(label_sets["A"]))
    if unexpected_adjudication_rows:
        raise ValueError(
            "Adjudication contains rows outside the annotator pool: "
            f"{unexpected_adjudication_rows[:5]}"
        )

    for key in sorted(label_sets["A"]):
        values = [label_sets[name][key] for name in ("A", "B", "C")]
        counts = Counter(values)
        agreed = majority_label(values)
        row = dict(row_sets["A"][key])
        row.update(
            {
                "annotator_A_label": values[0],
                "annotator_B_label": values[1],
                "annotator_C_label": values[2],
                "final_label": "",
                "final_label_source": "",
                "adjudication_notes": "",
            }
        )
        if agreed is None:
            adjudicated = adjudication.get(key)
            if not adjudicated:
                raise ValueError(f"Missing adjudication for {key[0]}/{key[1]}")
            row["final_label"] = int(adjudicated["adjudicated_label"])
            row["final_label_source"] = "adjudicated_no_majority"
            row["adjudication_notes"] = adjudicated.get("adjudication_notes", "")
            adjudicated_rows += 1
        else:
            row["final_label"] = agreed
            row["final_label_source"] = (
                "unanimous" if counts[agreed] == 3 else "majority"
            )
            if counts[agreed] == 3:
                unanimous_rows += 1
            majority_rows += 1
            if key in adjudication:
                raise ValueError(
                    f"Adjudication row provided despite majority for {key[0]}/{key[1]}"
                )
        labels[key] = int(row["final_label"])
        rows.append(row)

    summary = {
        "total_rows": len(rows),
        "unanimous_rows": unanimous_rows,
        "majority_rows": majority_rows,
        "adjudicated_no_majority_rows": adjudicated_rows,
        "final_label_counts": dict(Counter(labels.values())),
    }
    return labels, rows, summary


def write_metric_outputs(
    prefix: str,
    mapping: list[dict[str, Any]],
    labels: dict[tuple[str, str], int],
    out_dir: Path,
) -> dict[str, Any]:
    metrics = per_query_metrics(mapping, labels)
    summary = summarize_systems(metrics)
    pairwise = pairwise_comparisons(metrics)
    write_csv(out_dir / f"{prefix}_per_query_metrics.csv", metrics, list(metrics[0]))
    write_csv(out_dir / f"{prefix}_system_summary.csv", summary, list(summary[0]))
    write_csv(out_dir / f"{prefix}_pairwise_tests.csv", pairwise, list(pairwise[0]))
    return {"system_summary": summary, "pairwise_tests": pairwise}


def system_row(summary: list[dict[str, Any]], system: str) -> dict[str, Any]:
    for row in summary:
        if row["system"] == system:
            return row
    raise KeyError(system)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--packet-dir", default=str(DEFAULT_PACKET_DIR))
    parser.add_argument("--adjudication")
    parser.add_argument("--out-dir")
    args = parser.parse_args()

    packet_dir = Path(args.packet_dir)
    returned_dir = packet_dir / "RETURNED_LABELS"
    coordinator_dir = packet_dir / "COORDINATOR_ONLY"
    adjudication_path = (
        Path(args.adjudication)
        if args.adjudication
        else coordinator_dir / "ADJUDICATION_80Q_COMPLETED_Amr.csv"
    )
    mapping_path = coordinator_dir / "pool_mapping_DO_NOT_SEND_TO_ANNOTATORS.csv"
    out_dir = Path(args.out_dir) if args.out_dir else coordinator_dir / "three_annotator_scoring_output"
    out_dir.mkdir(parents=True, exist_ok=True)

    label_sets: dict[str, dict[tuple[str, str], int]] = {}
    row_sets: dict[str, dict[tuple[str, str], dict[str, str]]] = {}
    validation: dict[str, Any] = {"annotators": {}, "complete": False, "errors": []}
    for name, filename in ANNOTATOR_FILES.items():
        labels, rows, missing = read_annotator_sheet(returned_dir / filename, name)
        label_sets[name] = labels
        row_sets[name] = rows
        validation["annotators"][name] = {
            "file": str(returned_dir / filename),
            "labels": len(labels),
            "missing": missing,
            "label_counts": dict(Counter(labels.values())),
        }

    row_set_errors = validate_same_row_sets(label_sets, row_sets)
    validation["errors"].extend(row_set_errors)
    if row_set_errors:
        write_json(out_dir / "three_annotator_label_validation.json", validation)
        raise ValueError("; ".join(row_set_errors[:10]))

    adjudication = read_adjudication(adjudication_path)
    final_labels, consensus_rows, final_summary = build_final_consensus(
        label_sets,
        row_sets,
        adjudication,
    )
    mapping = read_mapping(mapping_path)
    mapped_keys = {(row["query_id"], row["pool_id"]) for row in mapping}
    if mapped_keys - set(final_labels):
        missing = sorted(mapped_keys - set(final_labels))
        raise ValueError(f"Missing final labels for mapped records: {missing[:5]}")
    validation["complete"] = True
    validation["adjudication"] = {
        "file": str(adjudication_path),
        "rows": len(adjudication),
    }
    validation["final_consensus"] = final_summary
    write_json(out_dir / "three_annotator_label_validation.json", validation)

    pairwise = {
        "A_vs_B": agreement_summary(label_sets["A"], label_sets["B"]),
        "A_vs_C": agreement_summary(label_sets["A"], label_sets["C"]),
        "B_vs_C": agreement_summary(label_sets["B"], label_sets["C"]),
    }
    agreement = {
        "total_rows": final_summary["total_rows"],
        "unanimous_rows": final_summary["unanimous_rows"],
        "unanimous_fraction": round(
            final_summary["unanimous_rows"] / final_summary["total_rows"],
            6,
        ),
        "majority_rows": final_summary["majority_rows"],
        "majority_fraction": round(
            final_summary["majority_rows"] / final_summary["total_rows"],
            6,
        ),
        "adjudicated_no_majority_rows": final_summary["adjudicated_no_majority_rows"],
        "adjudicated_no_majority_fraction": round(
            final_summary["adjudicated_no_majority_rows"]
            / final_summary["total_rows"],
            6,
        ),
        "fleiss_kappa_nominal": round(fleiss_kappa(label_sets), 6),
        "pairwise": pairwise,
    }
    write_json(out_dir / "three_annotator_agreement_final.json", agreement)

    fields = list(consensus_rows[0])
    write_csv(out_dir / "three_annotator_final_consensus_all_rows.csv", consensus_rows, fields)

    metrics_payload = write_metric_outputs("final_consensus", mapping, final_labels, out_dir)
    report = {
        "status": "three_annotator_80q_validation_complete",
        "packet_dir": str(packet_dir),
        "mapping": str(mapping_path),
        "agreement": agreement,
        "final_consensus": final_summary,
        "final_consensus_metrics": metrics_payload,
    }
    write_json(out_dir / "three_annotator_independent_validation_report.json", report)

    hybrid = system_row(metrics_payload["system_summary"], "bioscouter_hybrid")
    print(
        json.dumps(
            {
                "status": report["status"],
                "rows": final_summary["total_rows"],
                "queries": hybrid["queries"],
                "hybrid_p_at_10": hybrid["mean_p_at_10"],
                "hybrid_strict_p_at_10": hybrid["mean_strict_p_at_10"],
                "hybrid_ndcg_at_10": hybrid["mean_ndcg_at_10"],
                "fleiss_kappa_nominal": agreement["fleiss_kappa_nominal"],
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
