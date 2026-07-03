from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parent
V3 = ROOT.parent
EVALUATION = (
    V3
    / "evaluator_packet"
    / "BioScouter_Independent_Relevance_Evaluation_20260702"
    / "00_For_Ibrahim"
)
DEFAULT_A = EVALUATION / "returned_labels" / "BioScouter_Evaluator_A_Labeling_Sheet.csv"
DEFAULT_B = EVALUATION / "returned_labels" / "BioScouter_Evaluator_B_Labeling_Sheet.csv"
DEFAULT_ADJUDICATION = EVALUATION / "ADJUDICATION_TO_COMPLETE.csv"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def normalize_query(value: str) -> str:
    return " ".join(re.findall(r"[a-z0-9]+", value.casefold()))


def record_key(row: dict[str, str]) -> tuple[str, str, str]:
    return (
        normalize_query(row.get("query", "")),
        row.get("source", "").casefold().strip(),
        (row.get("accession") or row.get("dataset_id") or "").casefold().strip(),
    )


def label_value(row: dict[str, str], field: str) -> int | None:
    value = row.get(field, "").strip()
    return int(value) if value in {"0", "1", "2"} else None


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("controlled_run")
    parser.add_argument("--annotator-a", default=str(DEFAULT_A))
    parser.add_argument("--annotator-b", default=str(DEFAULT_B))
    parser.add_argument("--adjudication", default=str(DEFAULT_ADJUDICATION))
    args = parser.parse_args()

    run_dir = Path(args.controlled_run).resolve()
    new_pool = read_csv(run_dir / "annotator_A_blinded.csv")
    annotator_a = read_csv(Path(args.annotator_a))
    annotator_b = read_csv(Path(args.annotator_b))
    adjudication = read_csv(Path(args.adjudication))

    a_labels = {record_key(row): label_value(row, "annotator_label") for row in annotator_a}
    b_labels = {record_key(row): label_value(row, "annotator_label") for row in annotator_b}
    adjudicated = {
        record_key(row): label_value(row, "adjudicated_label")
        for row in adjudication
    }

    consensus: dict[tuple[str, str, str], int] = {}
    for key in a_labels.keys() & b_labels.keys():
        a_label = a_labels[key]
        b_label = b_labels[key]
        if a_label is None or b_label is None:
            continue
        if a_label == b_label:
            consensus[key] = a_label
        elif adjudicated.get(key) is not None:
            consensus[key] = int(adjudicated[key])

    output_rows: list[dict[str, object]] = []
    reused_rows: list[dict[str, object]] = []
    for row in sorted(new_pool, key=lambda item: (item["query_id"], item["pool_id"])):
        key = record_key(row)
        reused_label = consensus.get(key)
        output = {
            **row,
            "reused_consensus_label": "" if reused_label is None else reused_label,
            "reuse_status": "exact_query_record_match" if reused_label is not None else "requires_new_label",
        }
        output_rows.append(output)
        if reused_label is not None:
            reused_rows.append(output)

    pooled_fields = list(new_pool[0].keys()) + ["reused_consensus_label", "reuse_status"]
    write_csv(run_dir / "pooled_top10_for_labeling.csv", output_rows, pooled_fields)
    write_csv(run_dir / "reused_labels.csv", reused_rows, pooled_fields)

    report = [
        {
            "new_pooled_records": len(output_rows),
            "prior_consensus_records": len(consensus),
            "safely_reused_labels": len(reused_rows),
            "records_requiring_new_labels": len(output_rows) - len(reused_rows),
            "reuse_rule": "exact normalized query plus source plus accession",
            "interpretation": "Relevance labels are query-specific; accession-only reuse is prohibited.",
        }
    ]
    write_csv(
        run_dir / "label_reuse_report.csv",
        report,
        list(report[0].keys()),
    )
    (run_dir / "labeling_status.json").write_text(
        json.dumps(
            {
                **report[0],
                "independent_metrics_available": False,
                "reason": (
                    "The 80-query pool requires new query-specific labels. "
                    "Existing 30-query adjudicated metrics remain separate."
                ),
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(run_dir / "label_reuse_report.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
