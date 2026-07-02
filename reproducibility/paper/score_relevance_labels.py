"""Score blinded pooled relevance judgments for the controlled evaluation.

Labels are ordinal:
  2 = directly relevant
  1 = partially relevant
  0 = irrelevant or insufficiently described

The scorer fails closed when either annotator sheet is incomplete. Final
adjudicated metrics are only produced when a complete adjudication sheet is
provided explicitly.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import random
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any


VALID_LABELS = {"0", "1", "2"}
SYSTEM_ORDER = (
    "bioscouter_keyword",
    "bioscouter_embedding",
    "bioscouter_hybrid",
    "native_source_api",
    "omicsdi",
)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2), encoding="utf-8")


def read_labels(
    path: Path,
    *,
    label_field: str = "annotator_label",
    require_complete: bool = True,
) -> tuple[dict[tuple[str, str], int], list[str]]:
    labels: dict[tuple[str, str], int] = {}
    missing: list[str] = []
    for row_number, row in enumerate(read_csv(path), start=2):
        query_id = str(row.get("query_id", "")).strip()
        pool_id = str(row.get("pool_id", "")).strip()
        if not query_id or not pool_id:
            raise ValueError(f"{path}:{row_number} is missing query_id or pool_id")
        key = (query_id, pool_id)
        if key in labels:
            raise ValueError(f"{path}:{row_number} duplicates {query_id}/{pool_id}")
        raw_label = str(row.get(label_field, "")).strip()
        if raw_label not in VALID_LABELS:
            if not raw_label:
                missing.append(f"{query_id}/{pool_id}")
                continue
            raise ValueError(
                f"{path}:{row_number} has invalid {label_field}={raw_label!r}"
            )
        labels[key] = int(raw_label)
    if require_complete and missing:
        return labels, missing
    return labels, missing


def read_mapping(path: Path) -> list[dict[str, Any]]:
    rows = []
    seen: set[tuple[str, str, str]] = set()
    for row_number, row in enumerate(read_csv(path), start=2):
        key = (row["query_id"], row["system"], row["pool_id"])
        if key in seen:
            raise ValueError(f"{path}:{row_number} duplicates mapping {key}")
        seen.add(key)
        if row["system"] not in SYSTEM_ORDER:
            raise ValueError(f"{path}:{row_number} has unknown system {row['system']!r}")
        rank = int(row["rank"])
        if rank < 1:
            raise ValueError(f"{path}:{row_number} has invalid rank {rank}")
        rows.append({**row, "rank": rank})
    return rows


def dcg(labels: list[int], k: int = 10) -> float:
    padded = (labels + [0] * k)[:k]
    return sum(
        ((2**label) - 1) / math.log2(index + 2)
        for index, label in enumerate(padded)
    )


def per_query_metrics(
    mapping: list[dict[str, Any]],
    labels: dict[tuple[str, str], int],
    *,
    k: int = 10,
) -> list[dict[str, Any]]:
    pooled_labels: dict[str, dict[str, int]] = defaultdict(dict)
    rankings: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in mapping:
        key = (row["query_id"], row["pool_id"])
        if key not in labels:
            raise ValueError(f"Missing label for mapped record {key}")
        pooled_labels[row["query_id"]][row["pool_id"]] = labels[key]
        rankings[(row["query_id"], row["system"])].append(row)

    output = []
    query_ids = sorted(pooled_labels)
    for query_id in query_ids:
        ideal_labels = sorted(pooled_labels[query_id].values(), reverse=True)[:k]
        ideal_dcg = dcg(ideal_labels, k)
        total_relaxed_relevant = sum(
            1 for label in pooled_labels[query_id].values() if label >= 1
        )
        for system in SYSTEM_ORDER:
            rows = sorted(
                rankings.get((query_id, system), []),
                key=lambda row: row["rank"],
            )[:k]
            ranked_pool_ids = [row["pool_id"] for row in rows]
            ranked_labels = [pooled_labels[query_id][pool_id] for pool_id in ranked_pool_ids]
            relaxed_relevant = sum(1 for label in ranked_labels if label >= 1)
            strict_relevant = sum(1 for label in ranked_labels if label == 2)
            reciprocal_rank = 0.0
            for rank, label in enumerate(ranked_labels, start=1):
                if label == 2:
                    reciprocal_rank = 1.0 / rank
                    break
            output.append(
                {
                    "query_id": query_id,
                    "system": system,
                    "returned_at_k": len(ranked_labels),
                    "p_at_10": round(relaxed_relevant / k, 6),
                    "strict_p_at_10": round(strict_relevant / k, 6),
                    "ndcg_at_10": round(
                        dcg(ranked_labels, k) / ideal_dcg if ideal_dcg else 0.0,
                        6,
                    ),
                    "mrr_at_10": round(reciprocal_rank, 6),
                    "pooled_relevant_coverage_at_10": round(
                        relaxed_relevant / total_relaxed_relevant
                        if total_relaxed_relevant
                        else 0.0,
                        6,
                    ),
                    "zero_result": int(not ranked_labels),
                }
            )
    return output


def linear_weighted_kappa(labels_a: list[int], labels_b: list[int]) -> float:
    if not labels_a or len(labels_a) != len(labels_b):
        return float("nan")
    n = len(labels_a)
    observed_disagreement = sum(abs(a - b) / 2 for a, b in zip(labels_a, labels_b)) / n
    counts_a = Counter(labels_a)
    counts_b = Counter(labels_b)
    expected_disagreement = sum(
        (counts_a[a] / n) * (counts_b[b] / n) * (abs(a - b) / 2)
        for a in (0, 1, 2)
        for b in (0, 1, 2)
    )
    if expected_disagreement == 0:
        return 1.0 if observed_disagreement == 0 else 0.0
    return 1 - (observed_disagreement / expected_disagreement)


def unweighted_kappa(labels_a: list[int], labels_b: list[int]) -> float:
    n = len(labels_a)
    observed = sum(a == b for a, b in zip(labels_a, labels_b)) / n
    counts_a = Counter(labels_a)
    counts_b = Counter(labels_b)
    expected = sum(
        (counts_a[label] / n) * (counts_b[label] / n)
        for label in (0, 1, 2)
    )
    return 1.0 if expected == 1.0 else (observed - expected) / (1 - expected)


def agreement_summary(
    labels_a: dict[tuple[str, str], int],
    labels_b: dict[tuple[str, str], int],
) -> dict[str, Any]:
    if set(labels_a) != set(labels_b):
        only_a = sorted(set(labels_a) - set(labels_b))
        only_b = sorted(set(labels_b) - set(labels_a))
        raise ValueError(
            f"Annotator sheets do not contain the same pool records; "
            f"only A={only_a[:5]}, only B={only_b[:5]}"
        )
    keys = sorted(labels_a)
    values_a = [labels_a[key] for key in keys]
    values_b = [labels_b[key] for key in keys]
    exact = sum(a == b for a, b in zip(values_a, values_b)) / len(keys)
    adjacent = sum(abs(a - b) <= 1 for a, b in zip(values_a, values_b)) / len(keys)
    return {
        "n_paired": len(keys),
        "exact_percent_agreement": round(exact, 6),
        "within_one_category_agreement": round(adjacent, 6),
        "cohens_kappa_unweighted": round(unweighted_kappa(values_a, values_b), 6),
        "cohens_kappa_linear_weighted": round(
            linear_weighted_kappa(values_a, values_b),
            6,
        ),
    }


def bootstrap_mean_ci(
    values: list[float],
    *,
    iterations: int = 5000,
    seed: int = 20260702,
) -> tuple[float, float]:
    if not values:
        return float("nan"), float("nan")
    rng = random.Random(seed)
    means = []
    for _ in range(iterations):
        sample = [rng.choice(values) for _ in values]
        means.append(statistics.mean(sample))
    means.sort()
    low = means[int(0.025 * (len(means) - 1))]
    high = means[int(0.975 * (len(means) - 1))]
    return low, high


def summarize_systems(metrics: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in metrics:
        grouped[row["system"]].append(row)
    summary = []
    metric_names = (
        "p_at_10",
        "strict_p_at_10",
        "ndcg_at_10",
        "mrr_at_10",
        "pooled_relevant_coverage_at_10",
    )
    for system in SYSTEM_ORDER:
        rows = grouped.get(system, [])
        result: dict[str, Any] = {
            "system": system,
            "queries": len(rows),
            "zero_result_queries": sum(row["zero_result"] for row in rows),
        }
        for metric in metric_names:
            values = [float(row[metric]) for row in rows]
            low, high = bootstrap_mean_ci(values)
            result[f"mean_{metric}"] = round(statistics.mean(values), 6) if values else None
            result[f"{metric}_ci95_low"] = round(low, 6) if values else None
            result[f"{metric}_ci95_high"] = round(high, 6) if values else None
        summary.append(result)
    return summary


def paired_sign_permutation(
    differences: list[float],
    *,
    iterations: int = 20000,
    seed: int = 20260702,
) -> float:
    if not differences:
        return float("nan")
    observed = abs(statistics.mean(differences))
    rng = random.Random(seed)
    extreme = 0
    for _ in range(iterations):
        permuted = [value if rng.random() < 0.5 else -value for value in differences]
        if abs(statistics.mean(permuted)) >= observed - 1e-12:
            extreme += 1
    return (extreme + 1) / (iterations + 1)


def pairwise_comparisons(
    metrics: list[dict[str, Any]],
    *,
    baseline: str = "omicsdi",
) -> list[dict[str, Any]]:
    lookup = {
        (row["query_id"], row["system"]): row
        for row in metrics
    }
    query_ids = sorted({row["query_id"] for row in metrics})
    output = []
    for system in SYSTEM_ORDER:
        if system == baseline:
            continue
        for metric in ("p_at_10", "strict_p_at_10", "ndcg_at_10", "mrr_at_10"):
            differences = [
                float(lookup[(query_id, system)][metric])
                - float(lookup[(query_id, baseline)][metric])
                for query_id in query_ids
                if (query_id, system) in lookup and (query_id, baseline) in lookup
            ]
            output.append(
                {
                    "system": system,
                    "baseline": baseline,
                    "metric": metric,
                    "paired_queries": len(differences),
                    "mean_difference": round(statistics.mean(differences), 6),
                    "sign_permutation_p": round(paired_sign_permutation(differences), 6),
                }
            )
    return output


def make_adjudication_sheet(
    annotator_a_path: Path,
    labels_a: dict[tuple[str, str], int],
    labels_b: dict[tuple[str, str], int],
    out_path: Path,
) -> None:
    source_rows = {
        (row["query_id"], row["pool_id"]): row
        for row in read_csv(annotator_a_path)
    }
    fields = [
        "query_id",
        "query",
        "pool_id",
        "accession",
        "source",
        "title",
        "annotator_A_label",
        "annotator_B_label",
        "adjudicated_label",
        "adjudication_notes",
    ]
    rows = []
    for key in sorted(labels_a):
        source = source_rows[key]
        label_a = labels_a[key]
        label_b = labels_b[key]
        rows.append(
            {
                "query_id": key[0],
                "query": source.get("query", ""),
                "pool_id": key[1],
                "accession": source.get("accession", ""),
                "source": source.get("source", ""),
                "title": source.get("title", ""),
                "annotator_A_label": label_a,
                "annotator_B_label": label_b,
                "adjudicated_label": label_a if label_a == label_b else "",
                "adjudication_notes": "",
            }
        )
    write_csv(out_path, rows, fields)


def write_metric_outputs(
    prefix: str,
    mapping: list[dict[str, Any]],
    labels: dict[tuple[str, str], int],
    out_dir: Path,
) -> dict[str, Any]:
    metrics = per_query_metrics(mapping, labels)
    summary = summarize_systems(metrics)
    pairwise = pairwise_comparisons(metrics)
    write_csv(
        out_dir / f"{prefix}_per_query_metrics.csv",
        metrics,
        list(metrics[0].keys()),
    )
    write_csv(
        out_dir / f"{prefix}_system_summary.csv",
        summary,
        list(summary[0].keys()),
    )
    write_csv(
        out_dir / f"{prefix}_pairwise_tests.csv",
        pairwise,
        list(pairwise[0].keys()),
    )
    return {"system_summary": summary, "pairwise_tests": pairwise}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--annotator-a", required=True)
    parser.add_argument("--annotator-b", required=True)
    parser.add_argument("--mapping", required=True)
    parser.add_argument("--adjudicated")
    parser.add_argument("--out-dir")
    args = parser.parse_args()

    path_a = Path(args.annotator_a)
    path_b = Path(args.annotator_b)
    mapping_path = Path(args.mapping)
    out_dir = Path(args.out_dir) if args.out_dir else path_a.parent / "scoring"
    out_dir.mkdir(parents=True, exist_ok=True)

    labels_a, missing_a = read_labels(path_a)
    labels_b, missing_b = read_labels(path_b)
    validation = {
        "annotator_A_total_labels": len(labels_a),
        "annotator_B_total_labels": len(labels_b),
        "annotator_A_missing": missing_a,
        "annotator_B_missing": missing_b,
        "complete": not missing_a and not missing_b,
    }
    write_json(out_dir / "label_validation.json", validation)
    if missing_a or missing_b:
        print(
            f"Incomplete labels: annotator A missing {len(missing_a)}, "
            f"annotator B missing {len(missing_b)}",
            file=sys.stderr,
        )
        return 2

    mapping = read_mapping(mapping_path)
    agreement = agreement_summary(labels_a, labels_b)
    make_adjudication_sheet(
        path_a,
        labels_a,
        labels_b,
        out_dir / "adjudication_sheet.csv",
    )
    report: dict[str, Any] = {
        "status": "independent_labels_complete_not_adjudicated",
        "agreement": agreement,
        "annotator_A": write_metric_outputs("annotator_A", mapping, labels_a, out_dir),
        "annotator_B": write_metric_outputs("annotator_B", mapping, labels_b, out_dir),
        "final_adjudicated_metrics": None,
    }

    if args.adjudicated:
        adjudicated_path = Path(args.adjudicated)
        adjudicated_labels, missing_adjudicated = read_labels(
            adjudicated_path,
            label_field="adjudicated_label",
        )
        if missing_adjudicated:
            print(
                f"Adjudication incomplete: {len(missing_adjudicated)} labels missing",
                file=sys.stderr,
            )
            write_json(out_dir / "independent_validation_report.json", report)
            return 3
        if set(adjudicated_labels) != set(labels_a):
            raise ValueError("Adjudication sheet does not match the blinded pool")
        report["status"] = "independent_validation_complete"
        report["final_adjudicated_metrics"] = write_metric_outputs(
            "adjudicated",
            mapping,
            adjudicated_labels,
            out_dir,
        )

    write_json(out_dir / "independent_validation_report.json", report)
    print(json.dumps({"status": report["status"], "agreement": agreement}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
