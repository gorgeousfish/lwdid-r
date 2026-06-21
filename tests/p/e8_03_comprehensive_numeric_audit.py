#!/usr/bin/env python3
"""Hardening numeric audit for story E8-03 comprehensive sensitivity outputs."""

from __future__ import annotations

import json
import subprocess
import tempfile
import textwrap
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path("/Users/cxy/Desktop/lwdid_r")
PARITY_DIR = ROOT / "_automation" / "test-artifacts" / "parity"

OUT_JSON = PARITY_DIR / "20260324-qa-parity-e8-03-comprehensive-numeric-audit.json"
OUT_MD = PARITY_DIR / "20260324-qa-parity-e8-03-comprehensive-numeric-audit.md"

ATT_BOUND_MULTIPLIER = 10.0


R_SNIPPET = textwrap.dedent(
    rf"""
    suppressPackageStartupMessages(devtools::load_all("{(ROOT / 'lwdid-r').as_posix()}", quiet = TRUE))

    collect_numeric_scalars <- function(x, prefix = NULL) {{
      out <- list()

      if (is.list(x) && !is.data.frame(x)) {{
        nms <- names(x)
        if (is.null(nms)) {{
          nms <- rep("", length(x))
        }}

        for (i in seq_along(x)) {{
          child_name <- nms[[i]]
          if (!nzchar(child_name)) {{
            child_name <- sprintf("[[%d]]", i)
          }}
          child_prefix <- if (is.null(prefix)) {{
            child_name
          }} else {{
            paste0(prefix, "$", child_name)
          }}
          out <- c(out, collect_numeric_scalars(x[[i]], child_prefix))
        }}

        return(out)
      }}

      if (is.numeric(x) && length(x) == 1L) {{
        out[[prefix]] <- unname(as.numeric(x))
      }}

      out
    }}

    is_att_like_path <- function(path) {{
      grepl("\\$att$", path) ||
        grepl("^transformation_comparison\\$(demean_att|detrend_att)$", path) ||
        grepl("^estimator_comparison\\$(ra|ipw|ipwra|baseline_att)$", path)
    }}

    compute_numeric_audit <- function(result, y_sd, att_bound_multiplier = 10) {{
      numeric_map <- collect_numeric_scalars(result)
      numeric_values <- unlist(numeric_map, use.names = TRUE)

      nonfinite_paths <- names(numeric_values)[!is.finite(numeric_values)]
      att_paths <- names(numeric_values)[vapply(names(numeric_values), is_att_like_path, logical(1))]
      att_values <- numeric_values[att_paths]
      att_ratio_by_path <- abs(att_values) / y_sd
      out_of_range_paths <- names(att_ratio_by_path)[att_ratio_by_path >= att_bound_multiplier]

      list(
        numeric_field_count = length(numeric_values),
        nonfinite_paths = sort(unname(nonfinite_paths)),
        att_ratio_by_path = as.list(att_ratio_by_path),
        max_att_ratio = if (length(att_ratio_by_path) == 0L) {{
          0
        }} else {{
          max(att_ratio_by_path)
        }},
        out_of_range_paths = sort(unname(out_of_range_paths))
      )
    }}

    load_smoking_frozen_controls <- function() {{
      smoking_raw <- utils::read.csv(
        "{(ROOT / 'lwdid-py_v0.2.3' / 'data' / 'smoking.csv').as_posix()}",
        stringsAsFactors = FALSE
      )
      fixture <- utils::read.csv(
        "{(PARITY_DIR / 'e8_03_smoking_frozen_controls_fixture.csv').as_posix()}",
        stringsAsFactors = FALSE
      )
      controls <- c("lnincome", "beer", "age15to24", "lretprice")

      smoking_frozen <- merge(
        smoking_raw[, setdiff(names(smoking_raw), controls), drop = FALSE],
        fixture,
        by = "state",
        all.x = TRUE,
        sort = FALSE
      )

      smoking_frozen <- smoking_frozen[
        order(match(smoking_frozen$state, unique(smoking_raw$state)), smoking_frozen$year),
        ,
        drop = FALSE
      ]

      list(
        data = smoking_frozen,
        controls = controls
      )
    }}

    run_case <- function(case_id) {{
      warnings_seen <- character(0)

      case <- switch(
        case_id,
        "smoking-no-controls" = {{
          data("smoking", package = "lwdid", envir = environment())
          list(
            dataset_name = "smoking",
            description = "common-timing built-in data without controls",
            data = smoking,
            args = list(
              y = "lcigsale",
              ivar = "state",
              tvar = "year",
              d = "d",
              post = "post",
              controls = NULL
            )
          )
        }},
        "castle-no-controls" = {{
          data("castle", package = "lwdid", envir = environment())
          list(
            dataset_name = "castle",
            description = "staggered built-in data without controls",
            data = castle,
            args = list(
              y = "lhomicide",
              ivar = "sid",
              tvar = "year",
              gvar = "gvar",
              controls = NULL
            )
          )
        }},
        "smoking-frozen-controls" = {{
          frozen <- load_smoking_frozen_controls()
          list(
            dataset_name = "smoking-frozen-controls",
            description = "common-timing smoking with frozen unit-level controls",
            data = frozen$data,
            args = list(
              y = "lcigsale",
              ivar = "state",
              tvar = "year",
              d = "d",
              post = "post",
              controls = frozen$controls
            )
          )
        }},
        stop(sprintf("unknown hardening case: %s", case_id))
      )

      result <- withCallingHandlers(
        do.call(
          lwdid_sensitivity,
          c(
            list(
              data = case$data,
              type = "all",
              verbose = FALSE
            ),
            case$args
          )
        ),
        warning = function(w) {{
          warnings_seen <<- c(warnings_seen, conditionMessage(w))
          invokeRestart("muffleWarning")
        }}
      )

      y_sd <- stats::sd(case$data[[case$args$y]], na.rm = TRUE)
      audit <- compute_numeric_audit(
        result,
        y_sd = y_sd,
        att_bound_multiplier = {ATT_BOUND_MULTIPLIER}
      )

      status <- if (
        length(audit$nonfinite_paths) == 0L &&
        length(audit$out_of_range_paths) == 0L
      ) {{
        "passed"
      }} else {{
        "failed"
      }}

      list(
        id = case_id,
        dataset_name = case$dataset_name,
        description = case$description,
        att_bound_multiplier = {ATT_BOUND_MULTIPLIER},
        y_sd = y_sd,
        numeric_field_count = audit$numeric_field_count,
        nonfinite_paths = audit$nonfinite_paths,
        att_ratio_by_path = audit$att_ratio_by_path,
        max_att_ratio = audit$max_att_ratio,
        out_of_range_paths = audit$out_of_range_paths,
        warnings = sort(unique(unname(warnings_seen))),
        status = status
      )
    }}

    cases <- lapply(
      c("smoking-no-controls", "castle-no-controls", "smoking-frozen-controls"),
      run_case
    )
    all_pass <- all(vapply(cases, function(x) identical(x$status, "passed"), logical(1)))

    payload <- list(
      att_bound_multiplier = {ATT_BOUND_MULTIPLIER},
      comparison = list(
        status = if (all_pass) "passed" else "failed",
        numeric_status = if (all_pass) "passed" else "failed"
      ),
      cases = cases
    )

    cat(jsonlite::toJSON(payload, auto_unbox = TRUE, null = "null", digits = 17))
    """
)


def build_payload() -> dict:
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
        raise RuntimeError(f"Unable to locate audit JSON payload in stdout: {stdout}")

    payload = json.loads(stdout[start : end + 1])
    payload["story"] = "story-E8-03"
    payload["generated_at_utc"] = datetime.now(timezone.utc).isoformat()
    payload["commands"] = {
        "audit": f"python {OUT_JSON.with_name('e8_03_comprehensive_numeric_audit.py').as_posix()}",
        "test": (
            "Rscript -e 'devtools::test("
            f"pkg=\"{(ROOT / 'lwdid-r').as_posix()}\", "
            "filter=\"sensitivity-comprehensive-hardening-parity\", reporter=\"summary\")'"
        ),
    }
    return payload


def write_markdown(payload: dict) -> None:
    lines = [
        "# story-E8-03 comprehensive 数值 hardening audit",
        "",
        "## 结论",
        "",
    ]

    if payload["comparison"]["status"] == "passed":
        lines.append("- `confirmed`: 三条 comprehensive real-data 路径当前未出现非有限数值。")
        lines.append(
            f"- `confirmed`: 所有 ATT-like 标量当前均落在 `|ATT| < {ATT_BOUND_MULTIPLIER:g} * sd(y)` 的 hardening 边界内。"
        )
    else:
        lines.append("- `blocker`: comprehensive hardening audit 仍存在非有限数值或越界 ATT。")

    lines.extend(
        [
            "- `bounded`: 本证据覆盖 `smoking`、`castle` 与 frozen-controls 三条 `type=\"all\"` 返回路径，主要用于收敛 Task `8.6`。",
            "",
            "## 工件",
            "",
            "1. hardening audit",
            f"   - 脚本：`{(PARITY_DIR / 'e8_03_comprehensive_numeric_audit.py').relative_to(ROOT).as_posix()}`",
            f"   - JSON：`{OUT_JSON.relative_to(ROOT).as_posix()}`",
            f"   - 测试：`lwdid-r/tests/testthat/test-sensitivity-comprehensive-hardening-parity.R`",
            "",
            "## 案例摘要",
            "",
        ]
    )

    for case in payload["cases"]:
        lines.extend(
            [
                f"### {case['id']}",
                "",
                f"- 数据路径：`{case['dataset_name']}`",
                f"- 数值标量数：`{case['numeric_field_count']}`",
                f"- 最大 ATT/sd(y) 比值：`{case['max_att_ratio']:.6f}`",
                f"- 非有限路径数：`{len(case['nonfinite_paths'])}`",
                f"- 越界路径数：`{len(case['out_of_range_paths'])}`",
                "",
            ]
        )

    lines.extend(
        [
            "## 下一步建议",
            "",
            "1. 继续推进 Task `8.5`，把 FATAL inheritance 也固化为可执行 regression。",
            "2. 若后续补齐 `8.4` 其余 Stata 覆盖，应把同一 hardening audit 扩展到新增 oracle 所覆盖的分支。",
            "",
        ]
    )

    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    payload = build_payload()
    OUT_JSON.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    write_markdown(payload)
    print(f"Wrote {OUT_JSON}")
    print(f"Wrote {OUT_MD}")


if __name__ == "__main__":
    main()
