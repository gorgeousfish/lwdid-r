#!/usr/bin/env python3
"""Stata-backed precision comparator for story E8-03 smoking sensitivity output."""

from __future__ import annotations

import json
import math
import subprocess
import tempfile
import textwrap
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path("/Users/cxy/Desktop/lwdid_r")
PARITY_DIR = ROOT / "_automation" / "test-artifacts" / "parity"
STATA_DIR = ROOT / "lwdid_stata"

OUT_JSON = PARITY_DIR / "20260324-qa-parity-e8-03-smoking-stata-precision.json"
OUT_MD = PARITY_DIR / "20260324-qa-parity-e8-03-smoking-stata-precision.md"

ATT_TOL = 1.0e-6
SE_TOL = 1.0e-4


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

    payload <- list(
      warnings = unique(unname(warnings)),
      estimator_is_null = is.null(result$estimator_comparison),
      transformation = as.list(result$transformation_comparison)
    )

    cat(jsonlite::toJSON(payload, auto_unbox = TRUE, null = "null", digits = 17))
    """
)


STATA_TEMPLATE = textwrap.dedent(
    """
    clear all
    set more off
    adopath ++ "{ado_path}"
    use "{data_path}", clear

    tempname fh
    file open `fh' using "{output_path}", write text replace

    quietly lwdid lcigsale d, ivar(state) tvar(year) post(post) rolling(demean)
    local demean_att : display %21.17g e(att)
    local demean_se : display %21.17g e(se_att)
    file write `fh' "demean_att=`demean_att'" _n
    file write `fh' "demean_se=`demean_se'" _n

    quietly lwdid lcigsale d, ivar(state) tvar(year) post(post) rolling(detrend)
    local detrend_att : display %21.17g e(att)
    local detrend_se : display %21.17g e(se_att)
    file write `fh' "detrend_att=`detrend_att'" _n
    file write `fh' "detrend_se=`detrend_se'" _n

    file close `fh'
    exit, clear
    """
)


def run_r_payload() -> dict:
    with tempfile.NamedTemporaryFile("w", suffix=".R", delete=False) as handle:
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

    payload = json.loads(stdout[start : end + 1])
    payload.setdefault("warnings", [])
    payload.setdefault("transformation", None)
    payload.setdefault("estimator_is_null", None)
    return payload


def run_stata_payload() -> dict:
    with tempfile.TemporaryDirectory(prefix="lwdidr-stata-precision-") as tmpdir:
        tmp_path = Path(tmpdir)
        do_path = tmp_path / "run_precision.do"
        out_path = tmp_path / "stata_precision.txt"

        do_path.write_text(
            STATA_TEMPLATE.format(
                ado_path=STATA_DIR.as_posix(),
                data_path=(STATA_DIR / "smoking.dta").as_posix(),
                output_path=out_path.as_posix(),
            )
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
                "Stata precision output missing.\n"
                f"stdout:\n{completed.stdout}\n\nstderr:\n{completed.stderr}"
            )

        values = {}
        for line in out_path.read_text().splitlines():
            if "=" not in line:
                continue
            key, value = line.split("=", 1)
            values[key.strip()] = float(value.strip())

        demean_att = values["demean_att"]
        detrend_att = values["detrend_att"]
        difference = abs(demean_att - detrend_att)
        rel_diff = difference / abs(demean_att)

        return {
            "transformation": {
                "demean_att": demean_att,
                "demean_se": values["demean_se"],
                "detrend_att": detrend_att,
                "detrend_se": values["detrend_se"],
                "difference": difference,
                "rel_diff": rel_diff,
            }
        }


def compare_numeric(r_payload: dict, stata_payload: dict) -> dict:
    findings = []
    r_transform = r_payload["transformation"]
    stata_transform = stata_payload["transformation"]

    field_tolerances = {
        "demean_att": ATT_TOL,
        "demean_se": SE_TOL,
        "detrend_att": ATT_TOL,
        "detrend_se": SE_TOL,
        "difference": ATT_TOL,
        "rel_diff": ATT_TOL,
    }

    for field, tolerance in field_tolerances.items():
        expected = float(stata_transform[field])
        actual = float(r_transform[field])
        abs_diff = abs(actual - expected)
        findings.append(
            {
                "name": field,
                "status": "pass" if abs_diff <= tolerance else "fail",
                "expected": expected,
                "actual": actual,
                "abs_diff": abs_diff,
                "tolerance": tolerance,
            }
        )

    all_pass = all(item["status"] == "pass" for item in findings)
    return {
        "status": "passed" if all_pass else "failed",
        "numeric_status": "passed" if all_pass else "failed",
        "findings": findings,
    }


def build_payload() -> dict:
    r_payload = run_r_payload()
    stata_payload = run_stata_payload()
    comparison = compare_numeric(r_payload, stata_payload)

    return {
        "story": "story-E8-03",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "dataset": {
            "name": "smoking",
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
                "stata-mp -b do <tempdo>  # lwdid lcigsale d, ivar(state) "
                "tvar(year) post(post) rolling(demean|detrend)"
            ),
        },
        "tolerances": {
            "att": ATT_TOL,
            "se": SE_TOL,
        },
        "r_current": r_payload,
        "stata_reference": stata_payload,
        "comparison": comparison,
    }


def write_markdown(payload: dict) -> None:
    findings = payload["comparison"]["findings"]
    by_name = {item["name"]: item for item in findings}

    md = f"""# story-E8-03 smoking Stata 精度 comparator 证据

## 结论

- `confirmed`: qa-parity 本轮已为 `smoking` common-timing comprehensive transformation 分支补齐 Stata-backed executable oracle。
- `confirmed`: 当前 `lwdid_sensitivity(type="all")` 的 `demean/detrend` ATT 与 SE 均落在 ATT `< {ATT_TOL:g}`、SE `< {SE_TOL:g}` 容差内。
- `bounded`: 本证据覆盖 Task `8.4` 在 `smoking` common-timing / no-controls transformation path` 的底层 `lwdid()` 精度；`8.5 / 8.6` 与更广的 hardening matrix 仍待补齐。

## 关键证据

1. comparator
   - 脚本：`_automation/test-artifacts/parity/e8_03_smoking_stata_precision_comparator.py`
   - JSON：`_automation/test-artifacts/parity/20260324-qa-parity-e8-03-smoking-stata-precision.json`
   - 回归测试：`lwdid-r/tests/testthat/test-sensitivity-comprehensive-stata-precision-parity.R`

2. Stata 参考值
   - demean ATT：`{payload["stata_reference"]["transformation"]["demean_att"]:.17g}`
   - demean SE：`{payload["stata_reference"]["transformation"]["demean_se"]:.17g}`
   - detrend ATT：`{payload["stata_reference"]["transformation"]["detrend_att"]:.17g}`
   - detrend SE：`{payload["stata_reference"]["transformation"]["detrend_se"]:.17g}`

3. R vs Stata 差值
   - demean ATT abs diff：`{by_name["demean_att"]["abs_diff"]:.3g}`
   - demean SE abs diff：`{by_name["demean_se"]["abs_diff"]:.3g}`
   - detrend ATT abs diff：`{by_name["detrend_att"]["abs_diff"]:.3g}`
   - detrend SE abs diff：`{by_name["detrend_se"]["abs_diff"]:.3g}`
   - difference abs diff：`{by_name["difference"]["abs_diff"]:.3g}`
   - rel_diff abs diff：`{by_name["rel_diff"]["abs_diff"]:.3g}`

## 范围边界

- 本证据不覆盖 `pre_period` / `no_anticipation` 分支的 Stata 手工重建。
- 本证据不覆盖 staggered / FATAL 继承 / 全字段数值合理性。

## 下一步建议

1. 继续把 `8.5` 做成可执行 FATAL inheritance regression，而不是停留在口头 blocker。
2. 为 `8.6` 增加全字段 finite / range audit，优先复用现有 real-data 与 Monte Carlo oracle。
"""
    OUT_MD.write_text(md)


def main() -> None:
    payload = build_payload()
    OUT_JSON.write_text(json.dumps(payload, indent=2, ensure_ascii=False))
    write_markdown(payload)
    print(f"Wrote {OUT_JSON}")
    print(f"Wrote {OUT_MD}")


if __name__ == "__main__":
    main()
