from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any


VALID_LABELS = {"0", "1", "2"}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2), encoding="utf-8")


def row_key(row: dict[str, str]) -> tuple[str, str, str, str, str]:
    return (
        row.get("query_id", "").strip(),
        row.get("rank", "").strip(),
        row.get("dataset_id", "").strip(),
        row.get("accession", "").strip(),
        row.get("source", "").strip(),
    )


def validate_rows(path: Path, rows: list[dict[str, str]]) -> tuple[dict[tuple[str, str, str, str, str], int], list[dict[str, Any]]]:
    labels: dict[tuple[str, str, str, str, str], int] = {}
    problems: list[dict[str, Any]] = []
    required = {"query_id", "rank", "dataset_id", "accession", "source", "annotator_label"}
    missing_columns = sorted(required - set(rows[0].keys())) if rows else sorted(required)
    if missing_columns:
        problems.append({"type": "missing_columns", "columns": missing_columns})
        return labels, problems

    for row_number, row in enumerate(rows, start=2):
        key = row_key(row)
        if any(not value for value in key[:2]):
            problems.append({"type": "missing_key", "row": row_number, "key": key})
            continue
        if key in labels:
            problems.append({"type": "duplicate_key", "row": row_number, "key": key})
            continue
        raw = str(row.get("annotator_label", "")).strip()
        if raw not in VALID_LABELS:
            problems.append(
                {
                    "type": "invalid_or_missing_label",
                    "row": row_number,
                    "query_id": key[0],
                    "rank": key[1],
                    "value": raw,
                }
            )
            continue
        labels[key] = int(raw)
    return labels, problems


def dcg(labels: list[int], k: int = 10) -> float:
    padded = (labels + [0] * k)[:k]
    return sum(((2**label) - 1) / math.log2(index + 2) for index, label in enumerate(padded))


def per_query_metrics(rows: list[dict[str, str]], labels: dict[tuple[str, str, str, str, str], int], annotator: str) -> list[dict[str, Any]]:
    grouped: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for row in rows:
        key = row_key(row)
        if key not in labels:
            continue
        try:
            rank = int(row.get("rank", "0"))
        except ValueError:
            rank = 0
        grouped[row.get("query_id", "").strip()].append((rank, labels[key]))

    output = []
    for query_id in sorted(grouped):
        ranked = [label for _, label in sorted(grouped[query_id], key=lambda item: item[0])[:10]]
        ideal = sorted(ranked, reverse=True)
        ideal_dcg = dcg(ideal, 10)
        relaxed_relevant = sum(1 for label in ranked if label >= 1)
        strict_relevant = sum(1 for label in ranked if label == 2)
        output.append(
            {
                "annotator": annotator,
                "query_id": query_id,
                "returned_at_10": len(ranked),
                "p_at_10": round(relaxed_relevant / 10, 6),
                "strict_p_at_10": round(strict_relevant / 10, 6),
                "ndcg_at_10": round(dcg(ranked, 10) / ideal_dcg if ideal_dcg else 0.0, 6),
            }
        )
    return output


def mean(values: list[float]) -> float:
    return sum(values) / len(values) if values else 0.0


def summarize(metrics: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in metrics:
        grouped[row["annotator"]].append(row)
    return [
        {
            "annotator": annotator,
            "n_queries": len(rows),
            "mean_p_at_10": round(mean([float(row["p_at_10"]) for row in rows]), 6),
            "mean_strict_p_at_10": round(mean([float(row["strict_p_at_10"]) for row in rows]), 6),
            "mean_ndcg_at_10": round(mean([float(row["ndcg_at_10"]) for row in rows]), 6),
        }
        for annotator, rows in sorted(grouped.items())
    ]


def unweighted_kappa(values_a: list[int], values_b: list[int]) -> float:
    n = len(values_a)
    observed = sum(a == b for a, b in zip(values_a, values_b)) / n
    counts_a = Counter(values_a)
    counts_b = Counter(values_b)
    expected = sum((counts_a[label] / n) * (counts_b[label] / n) for label in (0, 1, 2))
    return 1.0 if expected == 1.0 else (observed - expected) / (1 - expected)


def linear_weighted_kappa(values_a: list[int], values_b: list[int]) -> float:
    n = len(values_a)
    observed_disagreement = sum(abs(a - b) / 2 for a, b in zip(values_a, values_b)) / n
    counts_a = Counter(values_a)
    counts_b = Counter(values_b)
    expected_disagreement = sum(
        (counts_a[a] / n) * (counts_b[b] / n) * (abs(a - b) / 2)
        for a in (0, 1, 2)
        for b in (0, 1, 2)
    )
    if expected_disagreement == 0:
        return 1.0 if observed_disagreement == 0 else 0.0
    return 1 - (observed_disagreement / expected_disagreement)


def make_disagreements(
    rows_a: list[dict[str, str]],
    labels_a: dict[tuple[str, str, str, str, str], int],
    labels_b: dict[tuple[str, str, str, str, str], int],
) -> list[dict[str, Any]]:
    disagreements = []
    for row in rows_a:
        key = row_key(row)
        if key not in labels_a or key not in labels_b:
            continue
        label_a = labels_a[key]
        label_b = labels_b[key]
        if label_a == label_b:
            continue
        disagreements.append(
            {
                "query_id": row.get("query_id", ""),
                "query": row.get("query", ""),
                "rank": row.get("rank", ""),
                "dataset_id": row.get("dataset_id", ""),
                "accession": row.get("accession", ""),
                "source": row.get("source", ""),
                "title": row.get("title", ""),
                "source_url": row.get("source_url", ""),
                "annotator_A_label": label_a,
                "annotator_B_label": label_b,
                "adjudicated_label": "",
                "adjudication_notes": "",
            }
        )
    return disagreements


def read_adjudication(
    path: Path,
    labels_a: dict[tuple[str, str, str, str, str], int],
    labels_b: dict[tuple[str, str, str, str, str], int],
) -> tuple[dict[tuple[str, str, str, str, str], int], list[dict[str, Any]]]:
    adjudicated: dict[tuple[str, str, str, str, str], int] = {}
    problems: list[dict[str, Any]] = []
    for row_number, row in enumerate(read_csv(path), start=2):
        key = row_key(row)
        raw = str(row.get("adjudicated_label", "")).strip()
        if raw not in VALID_LABELS:
            problems.append(
                {
                    "type": "invalid_or_missing_adjudicated_label",
                    "row": row_number,
                    "query_id": row.get("query_id", ""),
                    "rank": row.get("rank", ""),
                    "value": raw,
                }
            )
            continue
        if key not in labels_a or key not in labels_b:
            problems.append(
                {
                    "type": "adjudication_key_not_found",
                    "row": row_number,
                    "query_id": row.get("query_id", ""),
                    "rank": row.get("rank", ""),
                    "accession": row.get("accession", ""),
                    "source": row.get("source", ""),
                }
            )
            continue
        if labels_a[key] == labels_b[key]:
            problems.append(
                {
                    "type": "adjudication_row_not_a_disagreement",
                    "row": row_number,
                    "query_id": row.get("query_id", ""),
                    "rank": row.get("rank", ""),
                }
            )
            continue
        adjudicated[key] = int(raw)

    expected = {key for key in labels_a if labels_a[key] != labels_b[key]}
    missing = sorted(expected - set(adjudicated))
    if missing:
        problems.extend(
            {
                "type": "missing_disagreement_adjudication",
                "query_id": key[0],
                "rank": key[1],
                "accession": key[3],
                "source": key[4],
            }
            for key in missing[:50]
        )
    return adjudicated, problems


def consensus_labels(
    labels_a: dict[tuple[str, str, str, str, str], int],
    labels_b: dict[tuple[str, str, str, str, str], int],
    adjudicated: dict[tuple[str, str, str, str, str], int],
) -> dict[tuple[str, str, str, str, str], int]:
    output: dict[tuple[str, str, str, str, str], int] = {}
    for key in sorted(labels_a):
        if labels_a[key] == labels_b[key]:
            output[key] = labels_a[key]
        else:
            output[key] = adjudicated[key]
    return output


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--annotator-a", required=True)
    parser.add_argument("--annotator-b", required=True)
    parser.add_argument("--adjudication")
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()

    path_a = Path(args.annotator_a)
    path_b = Path(args.annotator_b)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    rows_a = read_csv(path_a)
    rows_b = read_csv(path_b)
    labels_a, problems_a = validate_rows(path_a, rows_a)
    labels_b, problems_b = validate_rows(path_b, rows_b)
    shared_keys = sorted(set(labels_a) & set(labels_b))
    only_a = sorted(set(labels_a) - set(labels_b))
    only_b = sorted(set(labels_b) - set(labels_a))

    validation = {
        "annotator_A_rows": len(rows_a),
        "annotator_B_rows": len(rows_b),
        "annotator_A_valid_labels": len(labels_a),
        "annotator_B_valid_labels": len(labels_b),
        "annotator_A_problems": problems_a,
        "annotator_B_problems": problems_b,
        "keys_only_in_A": only_a[:20],
        "keys_only_in_B": only_b[:20],
        "complete": not problems_a and not problems_b and not only_a and not only_b,
    }
    write_json(out_dir / "label_validation.json", validation)

    if not validation["complete"]:
        print("Labels are incomplete or invalid. See label_validation.json.", file=sys.stderr)
        return 2

    values_a = [labels_a[key] for key in shared_keys]
    values_b = [labels_b[key] for key in shared_keys]
    exact = sum(a == b for a, b in zip(values_a, values_b)) / len(shared_keys)
    within_one = sum(abs(a - b) <= 1 for a, b in zip(values_a, values_b)) / len(shared_keys)
    agreement = {
        "n_paired": len(shared_keys),
        "exact_percent_agreement": round(exact, 6),
        "within_one_category_agreement": round(within_one, 6),
        "cohens_kappa_unweighted": round(unweighted_kappa(values_a, values_b), 6),
        "cohens_kappa_linear_weighted": round(linear_weighted_kappa(values_a, values_b), 6),
    }
    write_json(out_dir / "agreement_summary.json", agreement)

    metrics = per_query_metrics(rows_a, labels_a, "Evaluator_A") + per_query_metrics(rows_b, labels_b, "Evaluator_B")
    write_csv(out_dir / "per_query_metrics.csv", metrics, list(metrics[0].keys()))
    summary = summarize(metrics)
    write_csv(out_dir / "system_summary.csv", summary, list(summary[0].keys()))

    disagreements = make_disagreements(rows_a, labels_a, labels_b)
    if disagreements:
        write_csv(out_dir / "disagreements_for_adjudication.csv", disagreements, list(disagreements[0].keys()))
    else:
        write_csv(
            out_dir / "disagreements_for_adjudication.csv",
            [],
            [
                "query_id",
                "query",
                "rank",
                "dataset_id",
                "accession",
                "source",
                "title",
                "source_url",
                "annotator_A_label",
                "annotator_B_label",
                "adjudicated_label",
                "adjudication_notes",
            ],
        )

    result: dict[str, Any] = {"agreement": agreement, "summary": summary}
    if args.adjudication:
        adjudicated, adjudication_problems = read_adjudication(
            Path(args.adjudication),
            labels_a,
            labels_b,
        )
        adjudication_validation = {
            "adjudication_rows": len(read_csv(Path(args.adjudication))),
            "valid_adjudicated_disagreements": len(adjudicated),
            "problems": adjudication_problems,
            "complete": not adjudication_problems,
        }
        write_json(out_dir / "adjudication_validation.json", adjudication_validation)
        if adjudication_problems:
            print("Adjudication is incomplete or invalid. See adjudication_validation.json.", file=sys.stderr)
            return 3
        final_labels = consensus_labels(labels_a, labels_b, adjudicated)
        final_metrics = per_query_metrics(rows_a, final_labels, "Adjudicated_Consensus")
        write_csv(out_dir / "adjudicated_per_query_metrics.csv", final_metrics, list(final_metrics[0].keys()))
        final_summary = summarize(final_metrics)
        write_csv(out_dir / "adjudicated_system_summary.csv", final_summary, list(final_summary[0].keys()))
        result["adjudicated_summary"] = final_summary

    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
