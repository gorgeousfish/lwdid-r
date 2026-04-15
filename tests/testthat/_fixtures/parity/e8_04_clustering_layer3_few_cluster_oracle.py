#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[3]
PARITY_DIR = Path(__file__).resolve().parent
PYTHON_SRC = REPO_ROOT / "lwdid-py_v0.2.3" / "src"

if str(PYTHON_SRC) not in sys.path:
    sys.path.insert(0, str(PYTHON_SRC))

from lwdid.clustering_diagnostics import (  # noqa: E402
    check_clustering_consistency,
    diagnose_clustering,
    recommend_clustering_level,
)


def to_builtin(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): to_builtin(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [to_builtin(v) for v in value]
    if isinstance(value, np.ndarray):
        return [to_builtin(v) for v in value.tolist()]
    if isinstance(value, np.generic):
        return value.item()
    if hasattr(value, "value") and value.__class__.__name__.endswith("Level"):
        return value.value
    return value


def few_cluster_data() -> pd.DataFrame:
    np.random.seed(42)

    n_regions = 5
    units_per_region = 200
    periods = 8
    rows = []

    for region in range(n_regions):
        treated = region < 2
        first_treat = 5 if treated else 0
        region_effect = np.random.normal(0, 3)

        for unit in range(units_per_region):
            unit_id = region * units_per_region + unit

            for year in range(1, periods + 1):
                post = year >= 5 if treated else False
                y = 10 + region_effect + 2.0 * post + np.random.normal(0, 1)
                rows.append(
                    {
                        "region": region,
                        "unit": unit_id,
                        "year": year,
                        "first_treat": first_treat,
                        "Y": float(y),
                    }
                )

    return pd.DataFrame(rows)


def write_fixture(df: pd.DataFrame, filename: str) -> str:
    path = PARITY_DIR / filename
    df.to_csv(path, index=False)
    return filename


def serialize_diagnosis(diag: Any) -> dict[str, Any]:
    return {
        "recommended_cluster_var": diag.recommended_cluster_var,
        "recommendation_reason": diag.recommendation_reason,
        "treatment_variation_level": diag.treatment_variation_level,
        "warnings": [str(w) for w in diag.warnings],
    }


def serialize_recommendation(rec: Any) -> dict[str, Any]:
    alternatives = []
    for alt in rec.alternatives:
        alt_dict = asdict(alt) if not isinstance(alt, dict) else dict(alt)
        alternatives.append(
            {
                "var": alt_dict.get("var_name", alt_dict.get("var")),
                "n_clusters": int(alt_dict["n_clusters"]),
                "reliability_score": float(alt_dict["reliability_score"]),
                "reason": alt_dict["reason"],
            }
        )

    return {
        "recommended_var": rec.recommended_var,
        "n_clusters": int(rec.n_clusters),
        "n_treated_clusters": int(rec.n_treated_clusters),
        "n_control_clusters": int(rec.n_control_clusters),
        "confidence": float(rec.confidence),
        "reasons": [str(x) for x in rec.reasons],
        "warnings": [str(x) for x in rec.warnings],
        "use_wild_bootstrap": bool(rec.use_wild_bootstrap),
        "wild_bootstrap_reason": rec.wild_bootstrap_reason,
        "alternatives": alternatives,
    }


def serialize_consistency(result: Any) -> dict[str, Any]:
    return {
        "is_consistent": bool(result.is_consistent),
        "treatment_variation_level": result.treatment_variation_level,
        "cluster_level": result.cluster_level,
        "n_clusters": int(result.n_clusters),
        "n_treatment_changes_within_cluster": int(
            result.n_treatment_changes_within_cluster
        ),
        "pct_clusters_with_variation": float(result.pct_clusters_with_variation),
        "recommendation": result.recommendation,
        "details": result.details,
    }


def serialize_diag_stats(diag: Any, cluster_var: str) -> dict[str, Any]:
    stats = diag.cluster_structure[cluster_var]
    return {
        "n_clusters": int(stats.n_clusters),
        "level_relative_to_unit": stats.level_relative_to_unit.value,
        "treatment_varies_within": bool(stats.treatment_varies_within_cluster),
        "reliability_score": float(stats.reliability_score),
    }


def main() -> None:
    fixture = few_cluster_data()
    fixture_name = write_fixture(
        fixture,
        "e8_04_clustering_layer3_few_cluster_fixture.csv",
    )

    diag = diagnose_clustering(
        fixture,
        ivar="unit",
        potential_cluster_vars=["region"],
        gvar="first_treat",
        verbose=False,
    )
    rec = recommend_clustering_level(
        fixture,
        ivar="unit",
        tvar="year",
        potential_cluster_vars=["region"],
        gvar="first_treat",
        min_clusters=20,
        verbose=False,
    )
    consistency = check_clustering_consistency(
        fixture,
        ivar="unit",
        cluster_var="region",
        gvar="first_treat",
        verbose=False,
    )

    oracle = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "story": "story-E8-04",
        "role": "qa-parity",
        "status": "layer-3-few-cluster-contract-established",
        "source_tests": [
            "lwdid-py_v0.2.3/tests/_archived_dev_tests/test_clustering_simulated.py::TestFewClusters::test_recommends_wild_bootstrap",
            "lwdid-py_v0.2.3/tests/_archived_dev_tests/test_clustering_empirical.py::TestRecommendClusteringLevelEmpirical::test_wild_bootstrap_recommendation",
            "lwdid-py_v0.2.3/tests/_archived_dev_tests/test_clustering_empirical.py::TestWorkflowIntegrationEmpirical::test_workflow_with_small_clusters",
        ],
        "cases": {
            "few_cluster": {
                "fixture_csv": fixture_name,
                "ivar": "unit",
                "tvar": "year",
                "gvar": "first_treat",
                "potential_cluster_vars": ["region"],
                "min_clusters": 20,
                "consistency_cluster_var": "region",
                "python_diagnose_clustering": serialize_diagnosis(diag),
                "diagnosis_cluster_structure": {
                    "region": serialize_diag_stats(diag, "region"),
                },
                "python_recommend_clustering": serialize_recommendation(rec),
                "python_consistency": serialize_consistency(consistency),
            }
        },
    }

    output_path = PARITY_DIR / "20260324-qa-parity-e8-04-layer3-few-cluster.json"
    with output_path.open("w", encoding="utf-8") as fh:
        json.dump(to_builtin(oracle), fh, indent=2, sort_keys=True)
        fh.write("\n")


if __name__ == "__main__":
    main()
