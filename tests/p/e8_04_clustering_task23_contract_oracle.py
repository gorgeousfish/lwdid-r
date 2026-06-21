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
    return value


def read_fixture(filename: str) -> pd.DataFrame:
    path = PARITY_DIR / filename
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite fixture: {path}")
    return pd.read_csv(path)


def cluster_unique_counts(
    data: pd.DataFrame,
    cluster_var: str,
    treatment_var: str,
) -> dict[str, int]:
    counts = data.groupby(cluster_var)[treatment_var].nunique().to_dict()
    return {str(k): int(v) for k, v in counts.items()}


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
            alternatives.append({str(k): to_builtin(v) for k, v in alt.items()})
        else:
            alternatives.append(to_builtin(asdict(alt)))

    return {
        "recommended_var": rec.recommended_var,
        "n_clusters": int(rec.n_clusters),
        "n_treated_clusters": int(rec.n_treated_clusters),
        "n_control_clusters": int(rec.n_control_clusters),
        "confidence": float(rec.confidence),
        "reasons": [str(x) for x in rec.reasons],
        "alternatives": alternatives,
        "warnings": [str(x) for x in rec.warnings],
        "use_wild_bootstrap": bool(rec.use_wild_bootstrap),
        "wild_bootstrap_reason": rec.wild_bootstrap_reason,
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


def main() -> None:
    small_cluster = read_fixture("e8_04_clustering_small_cluster_fixture.csv")
    hierarchical = read_fixture("e8_04_clustering_hierarchical_fixture.csv")
    never_treated = read_fixture("e8_04_clustering_never_treated_fixture.csv")

    oracle = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "story": "story-E8-04",
        "role": "qa-parity",
        "status": "task-3-public-contract-frozen",
        "source_tests": [
            "lwdid-py_v0.2.3/tests/_archived_dev_tests/root_level/test_clustering_diagnostics.py::TestDiagnoseClustering::test_recommends_wild_bootstrap_few_clusters",
            "lwdid-py_v0.2.3/tests/_archived_dev_tests/root_level/test_clustering_diagnostics.py::TestCheckClusteringConsistency::test_consistent_clustering",
            "lwdid-py_v0.2.3/tests/_archived_dev_tests/root_level/test_clustering_diagnostics.py::TestCheckClusteringConsistency::test_inconsistent_clustering",
        ],
        "cases": {
            "small_cluster": {
                "fixture_csv": "e8_04_clustering_small_cluster_fixture.csv",
                "ivar": "county",
                "gvar": "first_treat",
                "potential_cluster_vars": ["state"],
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
            "hierarchical_consistency": {
                "fixture_csv": "e8_04_clustering_hierarchical_fixture.csv",
                "ivar": "county",
                "gvar": "first_treat",
                "cluster_var": "state",
                "cluster_unique_counts": cluster_unique_counts(
                    hierarchical,
                    cluster_var="state",
                    treatment_var="first_treat",
                ),
                "python_consistency": serialize_consistency(
                    check_clustering_consistency(
                        hierarchical,
                        ivar="county",
                        cluster_var="state",
                        gvar="first_treat",
                        verbose=False,
                    )
                ),
            },
            "never_treated_consistency": {
                "fixture_csv": "e8_04_clustering_never_treated_fixture.csv",
                "ivar": "id",
                "gvar": "gvar",
                "cluster_var": "cluster",
                "cluster_unique_counts": cluster_unique_counts(
                    never_treated,
                    cluster_var="cluster",
                    treatment_var="gvar",
                ),
                "python_consistency": serialize_consistency(
                    check_clustering_consistency(
                        never_treated,
                        ivar="id",
                        cluster_var="cluster",
                        gvar="gvar",
                        verbose=False,
                    )
                ),
            },
        },
    }

    output_path = (
        PARITY_DIR / "20260324-qa-parity-e8-04-clustering-task23-contract.json"
    )
    with output_path.open("w", encoding="utf-8") as fh:
        json.dump(to_builtin(oracle), fh, indent=2, sort_keys=True)
        fh.write("\n")


if __name__ == "__main__":
    main()
