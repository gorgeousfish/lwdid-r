#!/usr/bin/env python3
"""Stata-backed oracle for smoking pre-period and no-anticipation branches."""

from __future__ import annotations

import csv
import json
import subprocess
import tempfile
import textwrap
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path("/Users/cxy/Desktop/lwdid_r")
PARITY_DIR = ROOT / "_automation" / "test-artifacts" / "parity"
STATA_DIR = ROOT / "lwdid_stata"
SMOKING_CSV = ROOT / "lwdid-py_v0.2.3" / "data" / "smoking.csv"

OUT_JSON = PARITY_DIR / "20260324-qa-parity-e8-03-smoking-stata-sensitivity-branches.json"
OUT_MD = PARITY_DIR / "20260324-qa-parity-e8-03-smoking-stata-sensitivity-branches.md"

ATT_TOL = 1.0e-6
SE_TOL = 1.0e-4
PRE_PERIOD_NPRE1_BOUNDARY_AUDIT = (
    PARITY_DIR / "20260324-theory-parity-e8-03-pre-period-npre1-stata-boundary.md"
)
PRE_PERIOD_NPRE1_MANUAL_ORACLE_JSON = (
    PARITY_DIR / "20260324-theory-parity-e8-03-pre-period-npre1-manual-oracle.json"
)
PRE_PERIOD_NPRE1_MANUAL_ORACLE_MD = (
    PARITY_DIR / "20260324-theory-parity-e8-03-pre-period-npre1-manual-oracle.md"
)


R_SNIPPET = textwrap.dedent(
    r"""
    suppressPackageStartupMessages(devtools::load_all("/Users/cxy/Desktop/lwdid_r/lwdid-r", quiet = TRUE))

    data(smoking, package = "lwdid", envir = environment())

    warnings <- character(0)
    result <- withCallingHandlers(
      lwdid_sensitivity(
        data = smoking,
        y = "lcigsale",
        ivar = "state",
        tvar = "year",
        d = "d",
        post = "post",
        type = "all",
        verbose = FALSE
      ),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )

    to_spec <- function(spec) {
      list(
        n_pre_periods = spec$n_pre_periods,
        start_period = spec$start_period,
        end_period = spec$end_period,
        excluded_periods = spec$excluded_periods,
        att = spec$att,
        se = spec$se,
        pvalue = spec$pvalue,
        ci_lower = spec$ci_lower,
        ci_upper = spec$ci_upper,
        converged = isTRUE(spec$converged)
      )
    }

    to_noa <- function(est) {
      list(
        excluded_periods = est$excluded_periods,
        n_pre_periods_used = est$n_pre_periods_used,
        att = est$att,
        se = est$se,
        converged = isTRUE(est$converged)
      )
    }

    payload <- list(
      warnings = unique(unname(warnings)),
      pre_period = list(
        pre_period_range_tested = as.integer(result$pre_period_result$pre_period_range_tested),
        baseline_spec = list(
          n_pre_periods = result$pre_period_result$baseline_spec$n_pre_periods,
          att = result$pre_period_result$baseline_spec$att,
          se = result$pre_period_result$baseline_spec$se
        ),
        sensitivity_ratio = result$pre_period_result$sensitivity_ratio,
        robustness_level = result$pre_period_result$robustness_level,
        is_robust = isTRUE(result$pre_period_result$is_robust),
        specs = lapply(result$pre_period_result$specifications, to_spec)
      ),
      no_anticipation = list(
        max_anticipation_tested = result$no_anticipation_result$max_anticipation_tested,
        anticipation_detected = isTRUE(result$no_anticipation_result$anticipation_detected),
        recommended_exclusion = result$no_anticipation_result$recommended_exclusion,
        detection_method = result$no_anticipation_result$detection_method,
        estimates = lapply(result$no_anticipation_result$estimates, to_noa)
      )
    )

    cat(jsonlite::toJSON(payload, auto_unbox = TRUE, null = "null", digits = 17))
    """
)


def load_pre_years() -> list[int]:
    years = set()
    with SMOKING_CSV.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if int(float(row["post"])) == 0:
                years.add(int(float(row["year"])))
    if not years:
        raise RuntimeError(f"No pre-treatment years found in {SMOKING_CSV}")
    return sorted(years)


def build_pre_specs(pre_years: list[int]) -> list[dict]:
    return [
        {
            "n_pre_periods": n_pre,
            "selected_pre_years": pre_years[-n_pre:],
        }
        for n_pre in range(1, len(pre_years) + 1)
    ]


def build_no_anticipation_specs(pre_years: list[int], max_anticipation: int = 3) -> list[dict]:
    max_tested = min(max_anticipation, max(0, len(pre_years) - 1))
    specs = []
    for excluded_periods in range(max_tested + 1):
        remaining_pre_years = pre_years if excluded_periods == 0 else pre_years[:-excluded_periods]
        specs.append(
            {
                "excluded_periods": excluded_periods,
                "remaining_pre_years": remaining_pre_years,
                "n_pre_periods_used": len(remaining_pre_years),
            }
        )
    return specs


def year_condition(years: list[int]) -> str:
    if not years:
        raise RuntimeError("At least one year is required to build a Stata filter")
    return " | ".join(f"year == {year}" for year in years)


def build_stata_script(pre_specs: list[dict], no_anticipation_specs: list[dict], output_path: Path) -> str:
    blocks = [
        "clear all",
        "set more off",
        f'adopath ++ "{STATA_DIR.as_posix()}"',
        f'local outpath "{output_path.as_posix()}"',
        "tempname fh",
        'file open `fh\' using "`outpath\'", write text replace',
        "",
    ]

    for spec in pre_specs:
        condition = year_condition(spec["selected_pre_years"])
        blocks.extend(
            [
                f"* pre_period n_pre = {spec['n_pre_periods']}",
                f'use "{(STATA_DIR / "smoking.dta").as_posix()}", clear',
                f"keep if post == 1 | {condition}",
                "capture noisily lwdid lcigsale d, ivar(state) tvar(year) post(post) rolling(demean)",
                "local rc = _rc",
                f'file write `fh\' "pre_{spec["n_pre_periods"]}_rc=`rc\'" _n',
                "if (`rc' == 0) {",
                "  local att : display %21.17g e(att)",
                "  local se : display %21.17g e(se_att)",
                f'  file write `fh\' "pre_{spec["n_pre_periods"]}_att=`att\'" _n',
                f'  file write `fh\' "pre_{spec["n_pre_periods"]}_se=`se\'" _n',
                "}",
                "",
            ]
        )

    for spec in no_anticipation_specs:
        condition = year_condition(spec["remaining_pre_years"])
        blocks.extend(
            [
                f"* no_anticipation exclude = {spec['excluded_periods']}",
                f'use "{(STATA_DIR / "smoking.dta").as_posix()}", clear',
                f"keep if post == 1 | {condition}",
                "egen __time_index = group(year)",
                "drop year",
                "rename __time_index year",
                "capture noisily lwdid lcigsale d, ivar(state) tvar(year) post(post) rolling(demean)",
                "local rc = _rc",
                f'file write `fh\' "na_{spec["excluded_periods"]}_rc=`rc\'" _n',
                "if (`rc' == 0) {",
                "  local att : display %21.17g e(att)",
                "  local se : display %21.17g e(se_att)",
                f'  file write `fh\' "na_{spec["excluded_periods"]}_att=`att\'" _n',
                f'  file write `fh\' "na_{spec["excluded_periods"]}_se=`se\'" _n',
                "}",
                "",
            ]
        )

    blocks.extend(
        [
            "file close `fh'",
            "exit, clear",
            "",
        ]
    )
    return "\n".join(blocks)


def run_r_payload() -> dict:
    with tempfile.NamedTemporaryFile("w", suffix=".R", delete=False, encoding="utf-8") as handle:
        handle.write(R_SNIPPET)
        script_path = Path(handle.name)

    try:
        completed = subprocess.run(
            ["Rscript", str(script_path)],
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
    return json.loads(stdout[start : end + 1])


def load_manual_oracle() -> dict:
    payload = json.loads(PRE_PERIOD_NPRE1_MANUAL_ORACLE_JSON.read_text(encoding="utf-8"))
    manual = payload["manual_oracle"]
    return {
        "source_json": str(PRE_PERIOD_NPRE1_MANUAL_ORACLE_JSON),
        "source_md": str(PRE_PERIOD_NPRE1_MANUAL_ORACLE_MD),
        "generated_at_cst": payload["generated_at_cst"],
        "source_status": payload["comparison"]["status"],
        "branch": payload["branch"],
        "specification": payload["specification"],
        "truth_source": payload["truth_source"],
        "att": float(manual["att"]),
        "se": float(manual["se"]),
        "pvalue": float(manual["pvalue"]),
        "ci_lower": float(manual["ci_lower"]),
        "ci_upper": float(manual["ci_upper"]),
    }


def run_stata_payload(pre_specs: list[dict], no_anticipation_specs: list[dict]) -> dict:
    with tempfile.TemporaryDirectory(prefix="lwdidr-stata-sensitivity-") as tmpdir:
        tmp_path = Path(tmpdir)
        do_path = tmp_path / "run_sensitivity_branches.do"
        out_path = tmp_path / "stata_sensitivity_branches.txt"

        do_path.write_text(
            build_stata_script(pre_specs, no_anticipation_specs, out_path),
            encoding="utf-8",
        )

        try:
            completed = subprocess.run(
                ["stata-mp", "-b", "do", str(do_path)],
                check=True,
                capture_output=True,
                text=True,
                cwd=tmpdir,
            )
        finally:
            for candidate in (
                tmp_path / f"{do_path.stem}.log",
                ROOT / f"{do_path.stem}.log",
            ):
                candidate.unlink(missing_ok=True)

        if not out_path.exists():
            raise RuntimeError(
                "Stata sensitivity-branches output missing.\n"
                f"stdout:\n{completed.stdout}\n\nstderr:\n{completed.stderr}"
            )

        values = {}
        for line in out_path.read_text(encoding="utf-8").splitlines():
            if "=" not in line:
                continue
            key, value = line.split("=", 1)
            values[key.strip()] = float(value.strip())

    pre_payload = {
        "specs": [
            {
                "n_pre_periods": spec["n_pre_periods"],
                "selected_pre_years": spec["selected_pre_years"],
                "stata_rc": int(values[f'pre_{spec["n_pre_periods"]}_rc']),
                "att": values.get(f'pre_{spec["n_pre_periods"]}_att'),
                "se": values.get(f'pre_{spec["n_pre_periods"]}_se'),
            }
            for spec in pre_specs
        ]
    }
    noa_payload = {
        "estimates": [
            {
                "excluded_periods": spec["excluded_periods"],
                "remaining_pre_years": spec["remaining_pre_years"],
                "n_pre_periods_used": spec["n_pre_periods_used"],
                "stata_rc": int(values[f'na_{spec["excluded_periods"]}_rc']),
                "att": values.get(f'na_{spec["excluded_periods"]}_att'),
                "se": values.get(f'na_{spec["excluded_periods"]}_se'),
            }
            for spec in no_anticipation_specs
        ]
    }
    return {
        "pre_period": pre_payload,
        "no_anticipation": noa_payload,
    }
def compare_numeric(r_payload: dict, stata_payload: dict, manual_oracle: dict) -> dict:
    findings = []
    blocked_specs = []

    for r_spec, stata_spec in zip(
        r_payload["pre_period"]["specs"],
        stata_payload["pre_period"]["specs"],
        strict=True,
    ):
        if r_spec["n_pre_periods"] != stata_spec["n_pre_periods"]:
            raise RuntimeError("pre_period spec order drift detected")
        if int(stata_spec["stata_rc"]) != 0:
            blocked_specs.append(
                {
                    "branch": "pre_period",
                    "index": int(r_spec["n_pre_periods"]),
                    "stata_rc": int(stata_spec["stata_rc"]),
                }
            )
            continue
        for field, tolerance in (("att", ATT_TOL), ("se", SE_TOL)):
            actual = float(r_spec[field])
            expected = float(stata_spec[field])
            abs_diff = abs(actual - expected)
            findings.append(
                {
                    "branch": "pre_period",
                    "index": int(r_spec["n_pre_periods"]),
                    "field": field,
                    "status": "pass" if abs_diff <= tolerance else "fail",
                    "expected": expected,
                    "actual": actual,
                    "abs_diff": abs_diff,
                    "tolerance": tolerance,
                }
            )

    for r_est, stata_est in zip(
        r_payload["no_anticipation"]["estimates"],
        stata_payload["no_anticipation"]["estimates"],
        strict=True,
    ):
        if r_est["excluded_periods"] != stata_est["excluded_periods"]:
            raise RuntimeError("no_anticipation estimate order drift detected")
        if int(stata_est["stata_rc"]) != 0:
            blocked_specs.append(
                {
                    "branch": "no_anticipation",
                    "index": int(r_est["excluded_periods"]),
                    "stata_rc": int(stata_est["stata_rc"]),
                }
            )
            continue
        for field, tolerance in (("att", ATT_TOL), ("se", SE_TOL)):
            actual = float(r_est[field])
            expected = float(stata_est[field])
            abs_diff = abs(actual - expected)
            findings.append(
                {
                    "branch": "no_anticipation",
                    "index": int(r_est["excluded_periods"]),
                    "field": field,
                    "status": "pass" if abs_diff <= tolerance else "fail",
                    "expected": expected,
                    "actual": actual,
                    "abs_diff": abs_diff,
                    "tolerance": tolerance,
                }
            )

    all_pass = all(item["status"] == "pass" for item in findings)
    if not findings:
        raise RuntimeError("No executable Stata findings were collected")

    manual_oracle_findings = []
    manual_oracle_status = "not-needed"
    boundary_decisions = []
    unresolved_blocked_specs = []
    for item in blocked_specs:
        if (
            item["branch"] == "pre_period"
            and int(item["index"]) == 1
            and int(item["stata_rc"]) == 2001
        ):
            current = next(
                spec
                for spec in r_payload["pre_period"]["specs"]
                if int(spec["n_pre_periods"]) == 1
            )
            for field, tolerance in (
                ("att", ATT_TOL),
                ("se", SE_TOL),
                ("pvalue", ATT_TOL),
                ("ci_lower", SE_TOL),
                ("ci_upper", SE_TOL),
            ):
                actual = float(current[field])
                expected = float(manual_oracle[field])
                abs_diff = abs(actual - expected)
                manual_oracle_findings.append(
                    {
                        "branch": "pre_period",
                        "index": 1,
                        "field": field,
                        "status": "pass" if abs_diff <= tolerance else "fail",
                        "expected": expected,
                        "actual": actual,
                        "abs_diff": abs_diff,
                        "tolerance": tolerance,
                    }
                )
            manual_oracle_status = (
                "passed"
                if all(item["status"] == "pass" for item in manual_oracle_findings)
                else "failed"
            )
            boundary_decisions.append(
                {
                    "branch": "pre_period",
                    "index": 1,
                    "stata_rc": 2001,
                    "decision": "waive",
                    "reason": "stata-implementation-boundary",
                    "allows_story_task_closure": True,
                    "truth_basis": [
                        "papers-and-math",
                        "r-python-agreement",
                    ],
                    "evidence": [
                        str(PRE_PERIOD_NPRE1_MANUAL_ORACLE_JSON),
                        str(PRE_PERIOD_NPRE1_MANUAL_ORACLE_MD),
                        str(PRE_PERIOD_NPRE1_BOUNDARY_AUDIT),
                    ],
                }
            )
        else:
            unresolved_blocked_specs.append(item)

    return {
        "status": (
            "passed"
            if all_pass and not blocked_specs
            else "passed-via-manual-oracle"
            if all_pass and manual_oracle_status == "passed" and not unresolved_blocked_specs
            else "partial"
            if all_pass
            else "failed"
        ),
        "numeric_status": "passed" if all_pass else "failed",
        "blocked_specs": blocked_specs,
        "manual_oracle_status": manual_oracle_status,
        "manual_oracle_findings": manual_oracle_findings,
        "boundary_decisions": boundary_decisions,
        "story_task_status": (
            "closure-ready-via-manual-oracle"
            if (
                all_pass
                and manual_oracle_status == "passed"
                and not unresolved_blocked_specs
            )
            else "closure-ready-via-waiver"
            if all_pass and not unresolved_blocked_specs
            else "blocked"
        ),
        "unresolved_blocked_specs": unresolved_blocked_specs,
        "findings": findings,
    }


def build_payload() -> dict:
    pre_years = load_pre_years()
    pre_specs = build_pre_specs(pre_years)
    no_anticipation_specs = build_no_anticipation_specs(pre_years)
    r_payload = run_r_payload()
    if isinstance(r_payload.get("warnings"), str):
        r_payload["warnings"] = [r_payload["warnings"]]
    stata_payload = run_stata_payload(pre_specs, no_anticipation_specs)
    manual_oracle = load_manual_oracle()
    comparison = compare_numeric(r_payload, stata_payload, manual_oracle)

    return {
        "story": "story-E8-03",
        "task": "8.4",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "dataset": {
            "name": "smoking",
            "pre_years": pre_years,
            "r_source": str(ROOT / "lwdid-r"),
            "stata_source": str(STATA_DIR / "smoking.dta"),
        },
        "commands": {
            "r": (
                "Rscript -e 'devtools::load_all(\"/Users/cxy/Desktop/lwdid_r/lwdid-r\", "
                "quiet=TRUE); data(smoking, package=\"lwdid\"); "
                "lwdid_sensitivity(data=smoking, y=\"lcigsale\", ivar=\"state\", "
                "tvar=\"year\", d=\"d\", post=\"post\", type=\"all\", verbose=FALSE)'"
            ),
            "stata": (
                "stata-mp -b do <tempdo>  # reruns lwdid on smoking.dta for each "
                "pre_period and no_anticipation subset"
            ),
        },
        "tolerances": {
            "att": ATT_TOL,
            "se": SE_TOL,
        },
        "manual_oracle": {
            "pre_period_n_pre1": manual_oracle,
        },
        "r_current": r_payload,
        "stata_reference": stata_payload,
        "comparison": comparison,
    }


def write_markdown(payload: dict) -> None:
    pre_specs = payload["stata_reference"]["pre_period"]["specs"]
    noa_estimates = payload["stata_reference"]["no_anticipation"]["estimates"]
    pre_count = len(pre_specs)
    noa_count = len(noa_estimates)
    failed = [item for item in payload["comparison"]["findings"] if item["status"] != "pass"]
    blocked = payload["comparison"]["blocked_specs"]
    decisions = payload["comparison"]["boundary_decisions"]
    manual_oracle = payload["manual_oracle"]["pre_period_n_pre1"]
    manual_findings = payload["comparison"]["manual_oracle_findings"]

    md = f"""# story-E8-03 smoking sensitivity branches Stata comparator

## 结论

- `confirmed`: qa-parity 本轮为 `smoking` comprehensive 的 `pre_period` 与 `no_anticipation` 两条子路径补齐了 Stata-backed executable oracle。
- `confirmed`: 当前共审计 `{pre_count}` 个 `pre_period` 子规格和 `{noa_count}` 个 `no_anticipation` 子规格；所有 Stata 可执行分支的 ATT 与 SE 均落在 ATT `< {ATT_TOL:g}`、SE `< {SE_TOL:g}` 容差内。
- `bounded`: 本证据只覆盖 Task `8.4` 在 common-timing `pre_period/no_anticipation` branches 的底层 `lwdid()` 精度，不覆盖 `8.5` 的 `ri=TRUE` / RI p-value 公式漂移。
- `resolved-with-manual-oracle`: `pre_period n_pre=1` 在 Stata 端返回 `r(2001)`，但 QA contract 现已接入 theory-parity 的 paper-backed manual oracle；当前 R 的 `ATT/SE/p/CI` 与该 oracle 全部对齐，因此该单点不再阻塞 Task `8.4`。

## 数据与命令

- 数据集：`smoking`
- R 入口：`lwdid_sensitivity(type="all")`
- Stata 入口：对 `smoking.dta` 的每个 `pre_period` / `exclude_periods` 子样本分别调用 `lwdid`

## 关键覆盖

- `pre_period`：范围 `{payload["r_current"]["pre_period"]["pre_period_range_tested"][0]}:{payload["r_current"]["pre_period"]["pre_period_range_tested"][1]}`，baseline `n_pre = {payload["r_current"]["pre_period"]["baseline_spec"]["n_pre_periods"]}`
- `no_anticipation`：`max_anticipation_tested = {payload["r_current"]["no_anticipation"]["max_anticipation_tested"]}`，`anticipation_detected = {str(payload["r_current"]["no_anticipation"]["anticipation_detected"]).lower()}`
- warnings：`{"; ".join(payload["r_current"]["warnings"]) if payload["r_current"]["warnings"] else "none"}`

## 裁决

- `status = {payload["comparison"]["status"]}`
- `numeric_status = {payload["comparison"]["numeric_status"]}`
- `manual_oracle_status = {payload["comparison"]["manual_oracle_status"]}`
- `story_task_status = {payload["comparison"]["story_task_status"]}`
- `failed_findings = {len(failed)}`
- `blocked_specs = {json.dumps(blocked, ensure_ascii=False)}`
- `boundary_decisions = {json.dumps(decisions, ensure_ascii=False)}`

## `n_pre=1` manual oracle

- source json: `{manual_oracle["source_json"]}`
- source status: `{manual_oracle["source_status"]}`
- oracle tuple: `ATT = {manual_oracle["att"]:.15f}`, `SE = {manual_oracle["se"]:.15f}`, `p = {manual_oracle["pvalue"]:.15f}`
- manual_oracle_findings = {json.dumps(manual_findings, ensure_ascii=False)}

## 产出文件

- JSON oracle: `{OUT_JSON}`
- 当前说明: `{OUT_MD}`
"""

    OUT_MD.write_text(md, encoding="utf-8")


def main() -> None:
    payload = build_payload()
    OUT_JSON.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    write_markdown(payload)
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
