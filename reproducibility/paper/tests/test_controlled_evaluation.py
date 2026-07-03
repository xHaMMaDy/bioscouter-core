from __future__ import annotations

import csv
import importlib.util
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def load_module(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / filename)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


controlled = load_module("controlled_evaluation", "controlled_evaluation.py")
scoring = load_module("score_relevance_labels", "score_relevance_labels.py")


def test_controlled_query_panel_is_disjoint_and_complete():
    validation = controlled.validate_inputs()
    assert validation["valid"] is True
    assert validation["controlled_queries"] == 30
    assert len(validation["omics_types"]) == 8


def test_eighty_query_panel_is_disjoint_and_complete():
    validation = controlled.validate_inputs(
        ROOT / "benchmark-queries-80.csv",
        expected_query_count=80,
    )
    assert validation["valid"] is True
    assert validation["controlled_queries"] == 80
    assert len(validation["omics_types"]) == 8


def test_p_at_10_uses_ten_result_denominator_and_pooled_ideal():
    mapping = []
    labels = {}
    for system in scoring.SYSTEM_ORDER:
        for rank in range(1, 11):
            pool_id = f"{system}-{rank}"
            mapping.append(
                {
                    "query_id": "Q1",
                    "pool_id": pool_id,
                    "system": system,
                    "rank": rank,
                }
            )
            labels[("Q1", pool_id)] = 2 if rank <= 3 else 0

    metrics = scoring.per_query_metrics(mapping, labels)
    keyword = next(row for row in metrics if row["system"] == "bioscouter_keyword")
    assert keyword["p_at_10"] == 0.3
    assert keyword["strict_p_at_10"] == 0.3
    assert keyword["ndcg_at_10"] == 0.469
    assert keyword["mrr_at_10"] == 1.0


def test_linear_weighted_kappa_rewards_adjacent_agreement():
    exact = scoring.linear_weighted_kappa([0, 1, 2], [0, 1, 2])
    adjacent = scoring.linear_weighted_kappa([0, 1, 2], [1, 2, 2])
    reversed_labels = scoring.linear_weighted_kappa([0, 1, 2], [2, 1, 0])
    assert exact == 1.0
    assert exact > adjacent > reversed_labels


def test_system_summaries_separate_ablation_and_baseline_rows(tmp_path):
    run_rows = []
    for system in controlled.SYSTEMS:
        run_rows.append(
            {
                "query_id": "Q1",
                "system": system,
                "status": 200,
                "elapsed_s": 1.25,
                "records": [
                    {
                        "record_key": f"{system}:A1",
                        "accession": "A1",
                        "source": "example",
                        "title": "Example dataset",
                        "omics_type": "transcriptomics",
                        "organism": "Homo sapiens",
                        "sample_count": 10,
                        "relevance_score": 0.9,
                        "source_url": "https://example.org/A1",
                    }
                ],
            }
        )

    controlled.write_system_summaries(run_rows, tmp_path, top_k=10)

    with (tmp_path / "system-summary.csv").open(newline="", encoding="utf-8") as handle:
        system_rows = list(csv.DictReader(handle))
    with (tmp_path / "ablation-summary.csv").open(newline="", encoding="utf-8") as handle:
        ablation_rows = list(csv.DictReader(handle))
    with (tmp_path / "baseline-comparison-summary.csv").open(newline="", encoding="utf-8") as handle:
        baseline_rows = list(csv.DictReader(handle))

    assert {row["system"] for row in system_rows} == set(controlled.SYSTEMS)
    assert {row["system"] for row in ablation_rows} == {
        "bioscouter_keyword",
        "bioscouter_embedding",
        "bioscouter_hybrid",
    }
    assert {row["system"] for row in baseline_rows} == {
        "bioscouter_hybrid",
        "native_source_api",
        "omicsdi",
    }
    assert all(row["returned_topk_records"] == "1" for row in system_rows)
