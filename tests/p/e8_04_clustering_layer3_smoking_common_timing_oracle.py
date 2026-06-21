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
SMOKING_CSV = REPO_ROOT / "lwdid-py_v0.2.3" / "data" / "smoking.csv"

if str(PYTHON_SRC) not in sys.path:
    sys.path.insert(0, str(PYTHON_SRC))

from lwdid.clustering_diagnostics import (  # noqa: E402
    _get_alternative_reason,
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


def smoking_common_timing_data() -> pd.DataFrame:
    return pd.read_csv(SMOKING_CSV)


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
        if "reason" in alt_dict:
            reason = alt_dict["reason"]
        else:
            reason = _get_alternative_reason(alt)
        alternatives.append(
            {
                "var": alt_dict.get("var_name", alt_dict.get("var")),
                "n_clusters": int(alt_dict["n_clusters"]),
                "reliability_score": float(alt_dict["reliability_score"]),
                "reason": reason,
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
    fixture = smoking_common_timing_data()

    diag = diagnose_clustering(
        fixture,
        ivar="state",
        potential_cluster_vars=["state"],
        gvar=None,
        d="d",
        verbose=False,
    )
    rec = recommend_clustering_level(
        fixture,
        ivar="state",
        tvar="year",
        potential_cluster_vars=["state"],
        gvar=None,
        d="d",
        min_clusters=20,
        verbose=False,
    )
    consistency = check_clustering_consistency(
        fixture,
        ivar="state",
        cluster_var="state",
        gvar=None,
        d="d",
        verbose=False,
    )

    oracle = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "story": "story-E8-04",
        "role": "story-worker",
        "status": "layer-3-smoking-common-timing-contract-established",
        "data_sources": [
            "lwdid-py_v0.2.3/data/smoking.csv",
            "lwdid-r/data/smoking.rda",
        ],
        "source_tests": [
            "lwdid-py_v0.2.3/tests/_archived_dev_tests/test_to_stata.py::common_timing_data",
            "lwdid-py_v0.2.3/tests/_archived_dev_tests/test_to_stata_numerical.py::smoking_data",
        ],
        "cases": {
            "smoking_common_timing": {
                "r_dataset": "smoking",
                "python_data_csv": "data/smoking.csv",
                "ivar": "state",
                "tvar": "year",
                "d": "d",
                "potential_cluster_vars": ["state"],
                "min_clusters": 20,
                "consistency_cluster_var": "state",
                "python_diagnose_clustering": serialize_diagnosis(diag),
                "diagnosis_cluster_structure": {
                    "state": serialize_diag_stats(diag, "state"),
                },
                "python_recommend_clustering": serialize_recommendation(rec),
                "python_consistency": serialize_consistency(consistency),
            }
        },
    }

    output_path = (
        PARITY_DIR / "20260324-story-worker-e8-04-layer3-smoking-common-timing.json"
    )
    with output_path.open("w", encoding="utf-8") as fh:
        json.dump(to_builtin(oracle), fh, indent=2, sort_keys=True)
        fh.write("\n")


if __name__ == "__main__":
    main()
