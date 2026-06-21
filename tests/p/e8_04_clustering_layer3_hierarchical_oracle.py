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


def hierarchical_panel_data() -> pd.DataFrame:
    np.random.seed(42)

    n_regions = 4
    states_per_region = 5
    industries_per_state = 10
    individuals_per_industry = 5
    n_periods = 10

    rows = []
    idcode = 0

    for region in range(n_regions):
        for state_idx in range(states_per_region):
            state = region * states_per_region + state_idx

            for industry_idx in range(industries_per_state):
                industry = state * industries_per_state + industry_idx
                first_treat = 6 if state < 10 else 0

                for _ in range(individuals_per_industry):
                    idcode += 1

                    for year in range(1, n_periods + 1):
                        state_effect = np.random.normal(0, 2)
                        treat_effect = 2.0 if (first_treat > 0 and year >= first_treat) else 0
                        y = 10 + 0.5 * year + state_effect + treat_effect + np.random.normal(0, 1)

                        rows.append(
                            {
                                "idcode": idcode,
                                "year": year,
                                "region": region,
                                "state": state,
                                "industry": industry,
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
    fixture = hierarchical_panel_data()
    fixture_name = write_fixture(
        fixture,
        "e8_04_clustering_layer3_hierarchical_fixture.csv",
    )

    potential_cluster_vars = ["state", "region"]
    min_clusters = 20

    diag = diagnose_clustering(
        fixture,
        ivar="idcode",
        potential_cluster_vars=potential_cluster_vars,
        gvar="first_treat",
        verbose=False,
    )
    rec = recommend_clustering_level(
        fixture,
        ivar="idcode",
        tvar="year",
        potential_cluster_vars=potential_cluster_vars,
        gvar="first_treat",
        min_clusters=min_clusters,
        verbose=False,
    )
    consistency = check_clustering_consistency(
        fixture,
        ivar="idcode",
        cluster_var=rec.recommended_var,
        gvar="first_treat",
        verbose=False,
    )

    oracle = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "story": "story-E8-04",
        "role": "qa-parity",
        "status": "layer-3-hierarchical-contract-established",
        "source_tests": [
            "lwdid-py_v0.2.3/tests/_archived_dev_tests/test_clustering_empirical.py::TestDiagnoseClusteringEmpirical::test_diagnose_with_hierarchical_data",
            "lwdid-py_v0.2.3/tests/_archived_dev_tests/test_clustering_empirical.py::TestRecommendClusteringLevelEmpirical::test_recommends_treatment_level",
            "lwdid-py_v0.2.3/tests/_archived_dev_tests/test_clustering_empirical.py::TestCheckClusteringConsistencyEmpirical::test_consistency_with_treatment_level_cluster",
            "lwdid-py_v0.2.3/tests/_archived_dev_tests/test_clustering_empirical.py::TestIntegrationEmpirical::test_full_workflow",
            "lwdid-py_v0.2.3/tests/_archived_dev_tests/test_clustering_empirical.py::TestIntegrationEmpirical::test_recommendation_matches_diagnostics",
        ],
        "cases": {
            "hierarchical": {
                "fixture_csv": fixture_name,
                "ivar": "idcode",
                "tvar": "year",
                "gvar": "first_treat",
                "potential_cluster_vars": potential_cluster_vars,
                "min_clusters": min_clusters,
                "consistency_cluster_var": rec.recommended_var,
                "python_diagnose_clustering": serialize_diagnosis(diag),
                "diagnosis_cluster_structure": {
                    "state": serialize_diag_stats(diag, "state"),
                    "region": serialize_diag_stats(diag, "region"),
                },
                "python_recommend_clustering": serialize_recommendation(rec),
                "python_consistency": serialize_consistency(consistency),
            }
        },
    }

    output_path = PARITY_DIR / "20260324-qa-parity-e8-04-layer3-hierarchical.json"
    with output_path.open("w", encoding="utf-8") as fh:
        json.dump(to_builtin(oracle), fh, indent=2, sort_keys=True)
        fh.write("\n")


if __name__ == "__main__":
    main()
