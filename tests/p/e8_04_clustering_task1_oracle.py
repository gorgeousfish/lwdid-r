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
    _analyze_cluster_var,
    check_clustering_consistency,
    diagnose_clustering,
    recommend_clustering_level,
)


def hierarchical_data() -> pd.DataFrame:
    rng = np.random.default_rng(42)
    rows = []
    n_states = 50
    counties_per_state = 10
    periods = 10

    for state in range(n_states):
        for county in range(counties_per_state):
            county_id = state * counties_per_state + county
            first_treat = 6 if state < 25 else 0
            for year in range(1, periods + 1):
                rows.append(
                    {
                        "state": state,
                        "county": county_id,
                        "year": year,
                        "first_treat": first_treat,
                        "Y": float(rng.normal(10 + 0.5 * year, 1)),
                    }
                )

    return pd.DataFrame(rows)


def small_cluster_data() -> pd.DataFrame:
    rng = np.random.default_rng(42)
    rows = []
    n_states = 5
    counties_per_state = 100

    for state in range(n_states):
        for county in range(counties_per_state):
            county_id = state * counties_per_state + county
            first_treat = 6 if state < 2 else 0
            for year in range(1, 11):
                rows.append(
                    {
                        "state": state,
                        "county": county_id,
                        "year": year,
                        "first_treat": first_treat,
                        "Y": float(rng.normal(10, 1)),
                    }
                )

    return pd.DataFrame(rows)


def treatment_varies_within_cluster_data() -> pd.DataFrame:
    rng = np.random.default_rng(42)
    return pd.DataFrame(
        {
            "state": [1, 1, 1, 1, 2, 2, 2, 2],
            "county": [1, 2, 3, 4, 5, 6, 7, 8],
            "year": [1, 1, 1, 1, 1, 1, 1, 1],
            "first_treat": [3, 3, 4, 4, 3, 4, 3, 4],
            "Y": rng.normal(10, 1, 8),
        }
    )


def never_treated_data() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "id": np.repeat(np.arange(1, 7), 5),
            "time": np.tile(np.arange(1, 6), 6),
            "gvar": np.repeat([3, 3, np.nan, 0, np.inf, 4], 5),
            "cluster": np.repeat(["A", "A", "B", "B", "C", "C"], 5),
        }
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


def compute_balance_score(n_clusters: int, n_treated: int, n_control: int) -> float:
    if n_clusters <= 0:
        return 0.0
    return float(min(min(n_treated, n_control) / (n_clusters / 2), 1.0))


def compute_cv_score(cv: float) -> float:
    if np.isnan(cv):
        return 0.0
    return float(max(0.0, 1.0 - cv / 2.0))


def serialize_cluster_stats(
    data: pd.DataFrame,
    cluster_var: str,
    ivar: str,
    gvar: str | None,
    d: str | None,
) -> dict[str, Any]:
    stats = _analyze_cluster_var(data, ivar=ivar, cluster_var=cluster_var, gvar=gvar, d=d)
    cluster_sizes = data.groupby(cluster_var).size().to_dict()
    balance_score = compute_balance_score(
        stats.n_clusters,
        stats.n_treated_clusters,
        stats.n_control_clusters,
    )
    cv_score = compute_cv_score(stats.cluster_size_cv)

    return {
        "var_name": stats.var_name,
        "n_clusters": int(stats.n_clusters),
        "cluster_sizes": {str(k): int(v) for k, v in cluster_sizes.items()},
        "min_size": int(stats.min_cluster_size),
        "max_size": int(stats.max_cluster_size),
        "mean_size": float(stats.mean_cluster_size),
        "median_size": float(stats.median_cluster_size),
        "cv": float(stats.cluster_size_cv),
        "n_treated_clusters": int(stats.n_treated_clusters),
        "n_control_clusters": int(stats.n_control_clusters),
        "treatment_varies_within": bool(stats.treatment_varies_within_cluster),
        "is_nested_in_unit": bool(stats.is_nested_in_unit),
        "level_relative_to_unit": stats.level_relative_to_unit.value,
        "units_per_cluster": float(stats.units_per_cluster),
        "n_clusters_with_variation": int(stats.n_clusters_with_treatment_variation),
        "balance_score": balance_score,
        "cv_score": cv_score,
        "reliability_score": float(stats.reliability_score),
        "is_valid": bool(stats.is_valid_cluster),
    }


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
        if isinstance(alt, dict):
            alt_dict = dict(alt)
        else:
            alt_dict = asdict(alt)
            if hasattr(alt, "level_relative_to_unit"):
                alt_dict["level_relative_to_unit"] = alt.level_relative_to_unit.value
        alternatives.append(to_builtin(alt_dict))

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
        "pct_clusters_with_variation": float(result.pct_clusters_with_variation),
        "n_treatment_changes_within_cluster": int(
            result.n_treatment_changes_within_cluster
        ),
        "n_clusters": int(result.n_clusters),
        "cluster_level": result.cluster_level,
        "treatment_variation_level": result.treatment_variation_level,
        "recommendation": result.recommendation,
        "details": result.details,
    }


def write_fixture(df: pd.DataFrame, filename: str) -> str:
    path = PARITY_DIR / filename
    df.to_csv(path, index=False)
    return filename


def main() -> None:
    hierarchical = hierarchical_data()
    small_cluster = small_cluster_data()
    within_cluster = treatment_varies_within_cluster_data()
    never_treated = never_treated_data()

    hierarchical_fixture = write_fixture(
        hierarchical,
        "e8_04_clustering_hierarchical_fixture.csv",
    )
    small_cluster_fixture = write_fixture(
        small_cluster,
        "e8_04_clustering_small_cluster_fixture.csv",
    )
    within_cluster_fixture = write_fixture(
        within_cluster,
        "e8_04_clustering_within_cluster_variation_fixture.csv",
    )
    never_treated_fixture = write_fixture(
        never_treated,
        "e8_04_clustering_never_treated_fixture.csv",
    )

    oracle = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "story": "story-E8-04",
        "role": "qa-parity",
        "status": "task-1-helper-parity-established",
        "source_tests": [
            "lwdid-py_v0.2.3/tests/_archived_dev_tests/root_level/test_clustering_diagnostics.py::TestDiagnoseClustering::test_detects_higher_level_cluster",
            "lwdid-py_v0.2.3/tests/_archived_dev_tests/root_level/test_clustering_diagnostics.py::TestDiagnoseClustering::test_recommends_treatment_level",
            "lwdid-py_v0.2.3/tests/_archived_dev_tests/root_level/test_clustering_diagnostics.py::TestDiagnoseClustering::test_detects_treatment_variation_within_cluster",
            "lwdid-py_v0.2.3/tests/_archived_dev_tests/root_level/test_clustering_diagnostics.py::TestRecommendClusteringLevel::test_recommends_wild_bootstrap_few_clusters",
        ],
        "cases": {
            "hierarchical": {
                "fixture_csv": hierarchical_fixture,
                "ivar": "county",
                "gvar": "first_treat",
                "cluster_levels": {
                    "state": "higher",
                    "county": "same",
                },
                "analyze_cluster_var": {
                    "state": serialize_cluster_stats(
                        hierarchical, "state", "county", "first_treat", None
                    ),
                    "county": serialize_cluster_stats(
                        hierarchical, "county", "county", "first_treat", None
                    ),
                },
                "python_diagnose_clustering": serialize_diagnosis(
                    diagnose_clustering(
                        hierarchical,
                        ivar="county",
                        potential_cluster_vars=["state", "county"],
                        gvar="first_treat",
                        verbose=False,
                    )
                ),
                "python_consistency_state": serialize_consistency(
                    check_clustering_consistency(
                        hierarchical,
                        ivar="county",
                        cluster_var="state",
                        gvar="first_treat",
                        verbose=False,
                    )
                ),
            },
            "within_cluster_variation": {
                "fixture_csv": within_cluster_fixture,
                "ivar": "county",
                "gvar": "first_treat",
                "analyze_cluster_var": {
                    "state": serialize_cluster_stats(
                        within_cluster, "state", "county", "first_treat", None
                    ),
                },
                "python_diagnose_clustering": serialize_diagnosis(
                    diagnose_clustering(
                        within_cluster,
                        ivar="county",
                        potential_cluster_vars=["state"],
                        gvar="first_treat",
                        verbose=False,
                    )
                ),
            },
            "never_treated": {
                "fixture_csv": never_treated_fixture,
                "ivar": "id",
                "gvar": "gvar",
                "analyze_cluster_var": {
                    "cluster": serialize_cluster_stats(
                        never_treated, "cluster", "id", "gvar", None
                    ),
                },
            },
            "small_cluster": {
                "fixture_csv": small_cluster_fixture,
                "ivar": "county",
                "gvar": "first_treat",
                "python_diagnose_clustering": serialize_diagnosis(
                    diagnose_clustering(
                        small_cluster,
                        ivar="county",
                        potential_cluster_vars=["state"],
                        gvar="first_treat",
                        verbose=False,
                    )
                ),
                "python_recommend_clustering_level": serialize_recommendation(
                    recommend_clustering_level(
                        small_cluster,
                        ivar="county",
                        tvar="year",
                        potential_cluster_vars=["state"],
                        gvar="first_treat",
                        min_clusters=20,
                        verbose=False,
                    )
                ),
            },
        },
    }

    oracle_path = PARITY_DIR / "20260324-qa-parity-e8-04-clustering-task1-python-oracle.json"
    with oracle_path.open("w", encoding="utf-8") as fh:
        json.dump(to_builtin(oracle), fh, indent=2, sort_keys=True)
        fh.write("\n")


if __name__ == "__main__":
    main()
