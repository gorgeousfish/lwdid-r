#!/usr/bin/env python3
"""Layer 4 Monte Carlo comparator for story E8-03 common-timing sensitivity."""

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
PARITY_DIR = ROOT / "_automation" / "test-artifacts" / "parity"

ATT_TOL = 1.0e-6
SE_TOL = 1.0e-4
BIAS_TOL = 5.0e-2
COVERAGE_TOL = 2.0e-2

N_PERIODS = 6
FIRST_TREAT = 4
BETA_T = np.array([1.0, 1.5, 0.8, 1.5, 2.0, 2.5])
LAMBDA_R = {4: 0.5, 5: 0.6, 6: 1.0}
BASE_EFFECT = (N_PERIODS - FIRST_TREAT + 1) * sum(
    1.0 / (r - FIRST_TREAT + 1) for r in range(FIRST_TREAT, N_PERIODS + 1)
)


R_SNIPPET = textwrap.dedent(
    r"""
    suppressPackageStartupMessages(devtools::load_all("/Users/cxy/Desktop/lwdid_r/lwdid-r", quiet = TRUE))

    fixture_path <- commandArgs(trailingOnly = TRUE)[[1]]
    fixture <- read.csv(fixture_path, stringsAsFactors = FALSE)

    replication_ids <- sort(unique(fixture$replication))

    run_with_warnings <- function(expr, collector) {
      withCallingHandlers(
        expr,
        warning = function(w) {
          collector$warnings <- c(collector$warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )
    }

    runs <- vector("list", length(replication_ids))

    for (idx in seq_along(replication_ids)) {
      rep_id <- replication_ids[[idx]]
      rep_df <- fixture[fixture$replication == rep_id, ]
      collector <- list(warnings = character(0))

      result <- tryCatch(
        run_with_warnings(
          lwdid_sensitivity(
            data = rep_df[c("id", "year", "y", "d", "post", "x1", "x2")],
            y = "y",
            ivar = "id",
            tvar = "year",
            d = "d",
            post = "post",
            controls = c("x1", "x2"),
            type = "all",
            verbose = FALSE
          ),
          collector
        ),
        error = function(e) e
      )

      treated_units <- length(unique(rep_df$id[rep_df$d == 1]))
      true_att <- unique(rep_df$true_att)
      seed <- unique(rep_df$seed)

      if (inherits(result, "error")) {
        runs[[idx]] <- list(
          replication = as.integer(rep_id),
          seed = as.integer(seed[[1L]]),
          treated_units = as.integer(treated_units),
          true_att = unname(true_att[[1L]]),
          error = conditionMessage(result),
          warning_count = length(collector$warnings),
          warnings = unname(collector$warnings)
        )
        next
      }

      demean_att <- unname(result$transformation_comparison$demean_att)
      demean_se <- unname(result$transformation_comparison$demean_se)
      ci_lower <- demean_att - 1.96 * demean_se
      ci_upper <- demean_att + 1.96 * demean_se

      runs[[idx]] <- list(
        replication = as.integer(rep_id),
        seed = as.integer(seed[[1L]]),
        treated_units = as.integer(treated_units),
        true_att = unname(true_att[[1L]]),
        demean_att = demean_att,
        demean_se = demean_se,
        bias = demean_att - unname(true_att[[1L]]),
        covered = (ci_lower <= unname(true_att[[1L]])) &&
          (unname(true_att[[1L]]) <= ci_upper),
        warning_count = length(collector$warnings),
        warnings = unname(collector$warnings),
        error = NULL
      )
    }

    successful <- Filter(function(item) is.null(item$error), runs)

    extract_numeric <- function(name) {
      vapply(successful, function(item) item[[name]], numeric(1))
    }

    summary <- if (length(successful) == 0L) {
      list(
        replications = 0L,
        mean_true_att = NULL,
        mean_att = NULL,
        mean_se = NULL,
        mean_bias = NULL,
        mean_abs_bias = NULL,
        coverage = NULL,
        finite_att = FALSE,
        finite_se = FALSE
      )
    } else {
      demean_att <- extract_numeric("demean_att")
      demean_se <- extract_numeric("demean_se")
      bias <- extract_numeric("bias")
      covered <- vapply(successful, function(item) isTRUE(item$covered), logical(1))
      true_att <- extract_numeric("true_att")

      list(
        replications = length(successful),
        mean_true_att = mean(true_att),
        mean_att = mean(demean_att),
        mean_se = mean(demean_se),
        mean_bias = mean(bias),
        mean_abs_bias = mean(abs(bias)),
        coverage = mean(covered),
        finite_att = all(is.finite(demean_att)),
        finite_se = all(is.finite(demean_se))
      )
    }

    payload <- list(
      summary = summary,
      by_replication = runs
    )

    cat(jsonlite::toJSON(payload, auto_unbox = TRUE, null = "null", digits = 17))
    """
)


def scenario_hx(x1: np.ndarray, x2: np.ndarray) -> np.ndarray:
    return (x1 - 4.0) / 2.0 + x2 / 3.0


def scenario_fx(x1: np.ndarray, x2: np.ndarray) -> np.ndarray:
    return (x1 - 4.0) / 3.0 + x2 / 2.0


def scenario_ps_index(x1: np.ndarray, x2: np.ndarray) -> np.ndarray:
    return -1.2 + (x1 - 4.0) / 2.0 - x2


def generate_replication(seed: int, n_units: int) -> pd.DataFrame:
    np.random.seed(seed)

    x1 = np.random.gamma(2.0, 2.0, n_units)
    x2 = np.random.binomial(1, 0.6, n_units)
    ps = 1.0 / (1.0 + np.exp(-scenario_ps_index(x1, x2)))
    d = np.random.binomial(1, ps, n_units)
    c = np.random.normal(2.0, 1.0, n_units)
    h_x = scenario_hx(x1, x2)
    f_x = scenario_fx(x1, x2)

    treated_mask = d == 1
    if treated_mask.sum() == 0 or treated_mask.sum() == n_units:
        raise ValueError("replication has degenerate treatment split")

    mean_h_treated = float(h_x[treated_mask].mean())
    lambda_mean = float(np.mean([LAMBDA_R[t] for t in range(FIRST_TREAT, N_PERIODS + 1)]))
    true_att = BASE_EFFECT + lambda_mean * mean_h_treated

    records: list[dict] = []
    for i in range(n_units):
        for t in range(1, N_PERIODS + 1):
            delta_t = float(t)
            u_0 = float(np.random.normal(0.0, 2.0))
            y_0 = delta_t + float(c[i]) + float(BETA_T[t - 1] * f_x[i]) + u_0
            post = 1 if t >= FIRST_TREAT else 0

            if d[i] == 1 and post == 1:
                tau = BASE_EFFECT + LAMBDA_R[t] * float(h_x[i])
                u_1 = float(np.random.normal(0.0, 2.0))
                y = y_0 + tau + u_1 - u_0
            else:
                y = y_0

            records.append(
                {
                    "id": i + 1,
                    "year": 2000 + t,
                    "y": y,
                    "d": int(d[i]),
                    "post": post,
                    "x1": float(x1[i]),
                    "x2": float(x2[i]),
                    "true_att": true_att,
                    "seed": seed,
                }
            )

    return pd.DataFrame.from_records(records)


def build_fixture(replications: int, n_units: int, seed_start: int, min_treated: int) -> pd.DataFrame:
    fixture_parts: list[pd.DataFrame] = []
    candidate_seed = seed_start

    while len(fixture_parts) < replications:
        rep_df = generate_replication(candidate_seed, n_units)
        treated_units = rep_df.loc[rep_df["d"] == 1, "id"].nunique()
        control_units = rep_df["id"].nunique() - treated_units

        if treated_units >= min_treated and control_units >= min_treated:
            rep_df = rep_df.copy()
            rep_df.insert(0, "replication", len(fixture_parts) + 1)
            fixture_parts.append(rep_df)

        candidate_seed += 1

    return pd.concat(fixture_parts, ignore_index=True)


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
    return value


def summarize_runs(runs: list[dict]) -> dict:
    successful = [run for run in runs if run.get("error") is None]

    if not successful:
        return {
            "summary": {
                "replications": 0,
                "mean_true_att": None,
                "mean_att": None,
                "mean_se": None,
                "mean_bias": None,
                "mean_abs_bias": None,
                "coverage": None,
                "finite_att": False,
                "finite_se": False,
            },
            "by_replication": runs,
        }

    mean_true_att = float(np.mean([run["true_att"] for run in successful]))
    mean_att = float(np.mean([run["demean_att"] for run in successful]))
    mean_se = float(np.mean([run["demean_se"] for run in successful]))
    mean_bias = float(np.mean([run["bias"] for run in successful]))
    mean_abs_bias = float(np.mean([abs(run["bias"]) for run in successful]))
    coverage = float(np.mean([1.0 if run["covered"] else 0.0 for run in successful]))
    finite_att = all(math.isfinite(run["demean_att"]) for run in successful)
    finite_se = all(math.isfinite(run["demean_se"]) for run in successful)

    return {
        "summary": {
            "replications": len(successful),
            "mean_true_att": mean_true_att,
            "mean_att": mean_att,
            "mean_se": mean_se,
            "mean_bias": mean_bias,
            "mean_abs_bias": mean_abs_bias,
            "coverage": coverage,
            "finite_att": finite_att,
            "finite_se": finite_se,
        },
        "by_replication": runs,
    }


def run_python_reference(fixture: pd.DataFrame) -> dict:
    sys.path.insert(0, str(PY_SRC))
    from lwdid.sensitivity import sensitivity_analysis

    runs: list[dict] = []

    for rep_id, rep_df in fixture.groupby("replication"):
        working = rep_df[["id", "year", "y", "d", "post", "x1", "x2"]].copy()
        true_att = float(rep_df["true_att"].iloc[0])
        seed = int(rep_df["seed"].iloc[0])
        treated_units = int(rep_df.loc[rep_df["d"] == 1, "id"].nunique())

        try:
            with warnings.catch_warnings(record=True) as captured:
                warnings.simplefilter("always")
                result = sensitivity_analysis(
                    working,
                    y="y",
                    ivar="id",
                    tvar="year",
                    d="d",
                    post="post",
                    controls=["x1", "x2"],
                    verbose=False,
                )

            transform = result.transformation_comparison
            demean_att = convert_scalar(transform["demean_att"])
            demean_se = convert_scalar(transform["demean_se"])
            ci_lower = demean_att - 1.96 * demean_se
            ci_upper = demean_att + 1.96 * demean_se

            runs.append(
                {
                    "replication": int(rep_id),
                    "seed": seed,
                    "treated_units": treated_units,
                    "true_att": true_att,
                    "demean_att": demean_att,
                    "demean_se": demean_se,
                    "bias": demean_att - true_att,
                    "covered": ci_lower <= true_att <= ci_upper,
                    "warning_count": len(captured),
                    "warnings": [str(item.message) for item in captured],
                    "error": None,
                }
            )
        except Exception as exc:  # pragma: no cover - comparator artifact
            runs.append(
                {
                    "replication": int(rep_id),
                    "seed": seed,
                    "treated_units": treated_units,
                    "true_att": true_att,
                    "error": str(exc),
                    "warning_count": 0,
                    "warnings": [],
                }
            )

    return summarize_runs(runs)


def run_r_reference(fixture_path: Path) -> dict:
    with tempfile.NamedTemporaryFile("w", suffix=".R", delete=False) as handle:
        handle.write(R_SNIPPET)
        script_path = Path(handle.name)

    try:
        completed = subprocess.run(
            ["Rscript", str(script_path), str(fixture_path)],
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
    return payload


def compare_summaries(python_data: dict, r_data: dict, expected_replications: int) -> dict:
    exact_findings: list[dict] = []
    numeric_findings: list[dict] = []
    blockers: list[str] = []

    py_summary = python_data["summary"]
    r_summary = r_data["summary"]

    py_errors = [item for item in python_data["by_replication"] if item.get("error") is not None]
    r_errors = [item for item in r_data["by_replication"] if item.get("error") is not None]

    if py_errors:
        blockers.append(f"Python Monte Carlo had {len(py_errors)} failed replications")
    if r_errors:
        blockers.append(f"R Monte Carlo had {len(r_errors)} failed replications")

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

    def add_numeric(name: str, expected: float, actual: float, tolerance: float) -> None:
        abs_diff = abs(actual - expected)
        status = "pass" if abs_diff <= tolerance else "fail"
        numeric_findings.append(
            {
                "name": name,
                "status": status,
                "expected": expected,
                "actual": actual,
                "abs_diff": abs_diff,
                "tolerance": tolerance,
            }
        )

    add_exact("dataset.replications", expected_replications, py_summary["replications"])
    add_exact("r.summary.replications", expected_replications, r_summary["replications"])
    add_exact("python.summary.finite_att", True, py_summary["finite_att"])
    add_exact("python.summary.finite_se", True, py_summary["finite_se"])
    add_exact("r.summary.finite_att", True, r_summary["finite_att"])
    add_exact("r.summary.finite_se", True, r_summary["finite_se"])

    add_numeric("summary.mean_att", py_summary["mean_att"], r_summary["mean_att"], ATT_TOL)
    add_numeric("summary.mean_se", py_summary["mean_se"], r_summary["mean_se"], SE_TOL)
    add_numeric("summary.mean_bias", py_summary["mean_bias"], r_summary["mean_bias"], ATT_TOL)
    add_numeric(
        "summary.mean_abs_bias",
        py_summary["mean_abs_bias"],
        r_summary["mean_abs_bias"],
        ATT_TOL,
    )
    add_numeric("summary.coverage", py_summary["coverage"], r_summary["coverage"], 1.0e-12)
    add_numeric(
        "python.acceptance.abs_mean_bias",
        0.0,
        abs(py_summary["mean_bias"]),
        BIAS_TOL,
    )
    add_numeric(
        "r.acceptance.abs_mean_bias",
        0.0,
        abs(r_summary["mean_bias"]),
        BIAS_TOL,
    )
    add_numeric(
        "python.acceptance.coverage_gap",
        0.0,
        abs(py_summary["coverage"] - 0.95),
        COVERAGE_TOL,
    )
    add_numeric(
        "r.acceptance.coverage_gap",
        0.0,
        abs(r_summary["coverage"] - 0.95),
        COVERAGE_TOL,
    )

    exact_status = "passed" if not blockers and all(
        item["status"] == "pass" for item in exact_findings
    ) else "failed"
    numeric_status = "passed" if not blockers and all(
        item["status"] == "pass" for item in numeric_findings
    ) else "failed"
    status = "passed" if exact_status == "passed" and numeric_status == "passed" else "failed"

    return {
        "status": status,
        "exact_status": exact_status,
        "numeric_status": numeric_status,
        "exact_findings": exact_findings,
        "numeric_findings": numeric_findings,
        "blockers": blockers,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--replications", type=int, default=50)
    parser.add_argument("--n-units", type=int, default=160)
    parser.add_argument("--seed-start", type=int, default=2401)
    parser.add_argument("--min-treated", type=int, default=25)
    parser.add_argument(
        "--fixture-output",
        type=Path,
        default=PARITY_DIR / "e8_03_monte_carlo_common_timing_fixture.csv",
    )
    parser.add_argument(
        "--json-output",
        type=Path,
        default=PARITY_DIR / "20260323-qa-parity-e8-03-monte-carlo-common-timing.json",
    )
    args = parser.parse_args()

    fixture = build_fixture(
        replications=args.replications,
        n_units=args.n_units,
        seed_start=args.seed_start,
        min_treated=args.min_treated,
    )
    args.fixture_output.parent.mkdir(parents=True, exist_ok=True)
    fixture.to_csv(args.fixture_output, index=False)

    python_data = run_python_reference(fixture)
    r_data = run_r_reference(args.fixture_output)
    comparison = compare_summaries(python_data, r_data, args.replications)

    payload = {
        "story": "story-E8-03",
        "dataset": {
            "name": "paper-dgp-common-timing-scenario-1c",
            "layer": "layer_4_monte_carlo",
            "scenario": "1C",
            "n_units": args.n_units,
            "n_periods": N_PERIODS,
            "first_treat_period": FIRST_TREAT,
            "replications": args.replications,
            "seed_start": args.seed_start,
            "fixture_csv": str(args.fixture_output),
            "metric_target": "transformation_comparison.demean_att",
        },
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "paper_reference": {
            "source": "Docs/lw2025.md Appendix C / lwdid-py test_paper_table_7_validation.py Scenario 1C",
            "notes": [
                "Uses the paper-backed common-timing DGP with controls x1/x2.",
                "Tracks Monte Carlo bias and 95% CI coverage for the comprehensive sensitivity demean branch.",
            ],
        },
        "tolerances": {
            "att": ATT_TOL,
            "se": SE_TOL,
            "monte_carlo_bias": BIAS_TOL,
            "monte_carlo_coverage": COVERAGE_TOL,
        },
        "python_reference": python_data,
        "r_result": r_data,
        "comparison": comparison,
    }

    args.json_output.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n")

    print(json.dumps({"status": comparison["status"], "json": str(args.json_output)}))
    return 0 if comparison["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
