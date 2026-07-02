from __future__ import annotations

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
