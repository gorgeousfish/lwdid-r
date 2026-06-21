#!/usr/bin/env python3
"""Shared-data real-data comparator for story E8-03 on castle.csv."""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import tempfile
import textwrap
import warnings
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path("/Users/cxy/Desktop/lwdid_r")
PY_SRC = ROOT / "lwdid-py_v0.2.3" / "src"
CASTLE_CSV = ROOT / "lwdid-py_v0.2.3" / "data" / "castle.csv"

ATT_TOL = 1.0e-6
SE_TOL = 1.0e-4
SR_TOL = 1.0e-6

PAPER_REFERENCE = {
    "table_3": {
        "demean_att": 0.092,
        "detrend_att": 0.067,
    }
}


R_SNIPPET = textwrap.dedent(
    r"""
    suppressPackageStartupMessages(devtools::load_all("/Users/cxy/Desktop/lwdid_r/lwdid-r", quiet = TRUE))

    csv_path <- commandArgs(trailingOnly = TRUE)[[1]]
    data <- read.csv(csv_path, stringsAsFactors = FALSE)
    if (!("gvar" %in% names(data)) && ("effyear" %in% names(data))) {
      data$gvar <- data$effyear
    }

    warnings <- character(0)

    run_with_warnings <- function(expr) {
      withCallingHandlers(
        expr,
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )
    }

    extract_result <- function(res) {
      list(
        overall_assessment = res$overall_assessment,
        recommendations = unname(as.character(res$recommendations)),
        pre_period = if (!is.null(res$pre_period_result)) {
          list(
            sensitivity_ratio = unname(res$pre_period_result$sensitivity_ratio),
            is_robust = isTRUE(res$pre_period_result$is_robust),
            robustness_level = as.character(res$pre_period_result$robustness_level)
          )
        } else {
          NULL
        },
        anticipation = if (!is.null(res$no_anticipation_result)) {
          list(
            anticipation_detected = isTRUE(res$no_anticipation_result$anticipation_detected),
            recommended_exclusion = unname(res$no_anticipation_result$recommended_exclusion)
          )
        } else {
          NULL
        },
        transformation = if (!is.null(res$transformation_comparison)) {
          as.list(res$transformation_comparison)
        } else {
          NULL
        },
        estimator = if (!is.null(res$estimator_comparison)) {
          as.list(res$estimator_comparison)
        } else {
          NULL
        }
      )
    }

    payload <- list(
      warnings = character(0),
      error = NULL,
      result = NULL
    )

    result <- tryCatch(
      run_with_warnings(
        lwdid_sensitivity(
          data = data,
          y = "lhomicide",
          ivar = "sid",
          tvar = "year",
          gvar = "gvar",
          type = "all",
          verbose = FALSE
        )
      ),
      error = function(e) e
    )

    payload$warnings <- unname(warnings)

    if (inherits(result, "error")) {
      payload$error <- conditionMessage(result)
    } else {
      payload$result <- extract_result(result)
    }

    cat(jsonlite::toJSON(payload, auto_unbox = TRUE, null = "null", digits = 17))
    """
)


def load_castle_fixture() -> pd.DataFrame:
    return pd.read_csv(CASTLE_CSV).copy()


def convert_scalar(value):
    if value is None:
        return None
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        if math.isnan(value):
            return None
        return float(value)
    if hasattr(value, "value"):
        return value.value
    return value


def python_payload(df: pd.DataFrame) -> dict:
    sys.path.insert(0, str(PY_SRC))
    from lwdid.sensitivity import sensitivity_analysis

    working = df.copy()
    working["gvar"] = working["effyear"].fillna(0).astype(int)

    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        result = sensitivity_analysis(
            working,
            y="lhomicide",
            ivar="sid",
            tvar="year",
            gvar="gvar",
            verbose=False,
        )

    return {
        "warnings": [str(item.message) for item in captured],
        "result": {
            "overall_assessment": result.overall_assessment,
            "recommendations": [str(item) for item in result.recommendations],
            "pre_period": {
                "sensitivity_ratio": convert_scalar(result.pre_period_result.sensitivity_ratio),
                "is_robust": convert_scalar(result.pre_period_result.is_robust),
                "robustness_level": convert_scalar(result.pre_period_result.robustness_level),
            }
            if result.pre_period_result is not None
            else None,
            "anticipation": {
                "anticipation_detected": convert_scalar(
                    result.anticipation_result.anticipation_detected
                ),
                "recommended_exclusion": convert_scalar(
                    result.anticipation_result.recommended_exclusion
                ),
            }
            if result.anticipation_result is not None
            else None,
            "transformation": {
                key: convert_scalar(value)
                for key, value in result.transformation_comparison.items()
            }
            if result.transformation_comparison is not None
            else None,
            "estimator": {
                key: convert_scalar(value)
                for key, value in result.estimator_comparison.items()
            }
            if result.estimator_comparison is not None
            else None,
        },
    }


def run_r_payload(csv_path: Path) -> dict:
    with tempfile.NamedTemporaryFile("w", suffix=".R", delete=False) as handle:
        handle.write(R_SNIPPET)
        script_path = Path(handle.name)

    try:
        completed = subprocess.run(
            ["Rscript", str(script_path), str(csv_path)],
            check=True,
            capture_output=True,
            text=True,
            cwd=str(ROOT),
        )
    finally:
        script_path.unlink(missing_ok=True)

    stdout = completed.stdout.strip()
    start = stdout.find("{")
    end = stdout.rfind("}")
    if start == -1 or end == -1 or end < start:
        raise RuntimeError(f"Unable to locate R JSON payload in stdout: {stdout}")
    payload = json.loads(stdout[start : end + 1])
    payload.setdefault("error", None)
    payload.setdefault("result", None)
    payload.setdefault("warnings", [])
    payload["stderr"] = completed.stderr.strip()
    return payload


def compare_payloads(python_data: dict, r_data: dict) -> dict:
    exact_findings = []
    numeric_findings = []
    blockers = []

    def add_exact(name: str, expected, actual) -> None:
        status = "pass" if expected == actual else "fail"
        exact_findings.append(
            {
                "name": name,
                "status": status,
                "expected": expected,
                "actual": actual,
            }
        )

    def add_numeric(name: str, expected, actual, tolerance: float) -> None:
        status = "pass" if abs(actual - expected) <= tolerance else "fail"
        numeric_findings.append(
            {
                "name": name,
                "status": status,
                "expected": expected,
                "actual": actual,
                "abs_diff": abs(actual - expected),
                "tolerance": tolerance,
            }
        )

    if r_data["error"] is not None:
        blockers.append(f"R comparator call failed: {r_data['error']}")
        return {
            "status": "blocked",
            "exact_status": "blocked",
            "numeric_status": "blocked",
            "exact_findings": exact_findings,
            "numeric_findings": numeric_findings,
            "blockers": blockers,
        }

    py_result = python_data["result"]
    r_result = r_data["result"]

    add_exact(
        "pre_period.is_robust",
        py_result["pre_period"]["is_robust"],
        r_result["pre_period"]["is_robust"],
    )
    add_exact(
        "pre_period.robustness_level",
        py_result["pre_period"]["robustness_level"],
        r_result["pre_period"]["robustness_level"],
    )
    add_exact(
        "anticipation.anticipation_detected",
        py_result["anticipation"]["anticipation_detected"],
        r_result["anticipation"]["anticipation_detected"],
    )
    add_exact(
        "anticipation.recommended_exclusion",
        py_result["anticipation"]["recommended_exclusion"],
        r_result["anticipation"]["recommended_exclusion"],
    )
    add_exact(
        "estimator.is_null",
        py_result["estimator"] is None,
        r_result["estimator"] is None,
    )

    add_numeric(
        "pre_period.sensitivity_ratio",
        py_result["pre_period"]["sensitivity_ratio"],
        r_result["pre_period"]["sensitivity_ratio"],
        SR_TOL,
    )

    for key in ("demean_att", "detrend_att", "difference"):
        add_numeric(
            f"transformation.{key}",
            py_result["transformation"][key],
            r_result["transformation"][key],
            ATT_TOL,
        )
    for key in ("demean_se", "detrend_se"):
        add_numeric(
            f"transformation.{key}",
            py_result["transformation"][key],
            r_result["transformation"][key],
            SE_TOL,
        )

    exact_failed = any(item["status"] != "pass" for item in exact_findings)
    numeric_failed = any(item["status"] != "pass" for item in numeric_findings)
    status = "failed" if exact_failed or numeric_failed else "passed"

    return {
        "status": status,
        "exact_status": "failed" if exact_failed else "passed",
        "numeric_status": "failed" if numeric_failed else "passed",
        "exact_findings": exact_findings,
        "numeric_findings": numeric_findings,
        "blockers": blockers,
    }


def paper_alignment(r_data: dict) -> dict:
    if r_data["result"] is None or r_data["result"]["transformation"] is None:
        return {"status": "blocked", "details": "missing R transformation result"}

    trans = r_data["result"]["transformation"]
    ref = PAPER_REFERENCE["table_3"]
    return {
        "demean_att_abs_diff_vs_table3": abs(trans["demean_att"] - ref["demean_att"]),
        "detrend_att_abs_diff_vs_table3": abs(trans["detrend_att"] - ref["detrend_att"]),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        required=True,
        help="Path to the JSON report to write.",
    )
    args = parser.parse_args()

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fixture = load_castle_fixture()
    with tempfile.NamedTemporaryFile("w", suffix=".csv", delete=False) as handle:
        fixture_path = Path(handle.name)
        fixture.to_csv(handle.name, index=False)

    try:
        python_data = python_payload(fixture)
        r_data = run_r_payload(fixture_path)
    finally:
        fixture_path.unlink(missing_ok=True)

    cohort_series = fixture.groupby("sid")["effyear"].first()
    report = {
        "story": "story-E8-03",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "dataset": {
            "source": str(CASTLE_CSV),
            "rows": int(len(fixture)),
            "units": int(fixture["sid"].nunique()),
            "treated_units": int(cohort_series.notna().sum()),
            "never_treated_units": int(cohort_series.isna().sum()),
            "years": [int(fixture["year"].min()), int(fixture["year"].max())],
            "cohorts": [int(value) for value in sorted(cohort_series.dropna().unique())],
        },
        "paper_reference": PAPER_REFERENCE,
        "python_reference": python_data,
        "r_result": r_data,
        "comparison": compare_payloads(python_data, r_data),
        "paper_alignment": paper_alignment(r_data),
    }

    output_path.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    print(output_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
