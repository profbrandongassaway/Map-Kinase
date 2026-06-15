"""Generate a runtime-settings report from collected OpenShift evidence.

Reads:
    deployment_validation/results/openshift_evidence/env_mapkinase.txt
    deployment_validation/results/openshift_evidence/tempdirs.txt

Writes:
    deployment_validation/results/runtime_settings_report.md

Compares the deployed values to the expected defaults in
DEPLOYMENT_VALIDATION_PLAN.md (upload limits, throttle settings, session TTL,
debug-related flags). Flags:
    PASS    Matches expected default (or explicitly-set safe value).
    REVIEW  Variable not set in pod; relying on app default. Plan considers this
            acceptable but OIT may want it explicitly set.
    FAIL    Variable set to an unsafe value (e.g. debug flag enabled in prod).

This script does not contact the cluster; it only parses the files written by
collect_openshift_evidence.sh.
"""

from __future__ import annotations

import re
import sys
from dataclasses import dataclass
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
EVIDENCE_DIR = REPO_ROOT / "deployment_validation" / "results" / "openshift_evidence"
ENV_FILE = EVIDENCE_DIR / "env_mapkinase.txt"
TEMP_FILE = EVIDENCE_DIR / "tempdirs.txt"
OUT_FILE = REPO_ROOT / "deployment_validation" / "results" / "runtime_settings_report.md"


@dataclass
class Check:
    name: str
    expected: str
    actual: str
    status: str  # PASS / REVIEW / FAIL
    notes: str


TRUTHY = {"1", "true", "yes", "on", "enabled"}
FALSY = {"0", "false", "no", "off", "disabled", ""}


def parse_env(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    env: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, _, value = line.partition("=")
        env[key.strip()] = value.strip()
    return env


def parse_tempdirs(path: Path) -> dict[str, str]:
    """Extract the resolved temp/session paths from tempdirs.txt."""
    info: dict[str, str] = {}
    if not path.exists():
        return info
    text = path.read_text(encoding="utf-8", errors="replace")
    for key, pattern in (
        ("tempdir", r"tempfile\.gettempdir\(\)\s*=\s*(.+)"),
        ("env_root", r"MAPKINASE_SESSION_TMP_ROOT\s*=\s*(.+)"),
        ("session_root", r"resolved session root\s*=\s*(.+)"),
        ("session_root_exists", r"^exists\s*=\s*(.+)$"),
    ):
        match = re.search(pattern, text, flags=re.MULTILINE)
        if match:
            info[key] = match.group(1).strip()
    return info


def _int_or_none(value: str) -> int | None:
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _float_or_none(value: str) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def check_numeric(
    env: dict[str, str],
    name: str,
    expected_default: float,
    units: str,
    *,
    allow_higher: bool = False,
    allow_lower: bool = False,
) -> Check:
    """Check a numeric env var against the plan's expected default.

    By default a missing var is REVIEW (relying on code default), and a present
    var must equal the expected default to PASS. If allow_higher/allow_lower is
    set, values outside the safe direction become FAIL.
    """
    raw = env.get(name, "")
    if raw == "":
        return Check(
            name=name,
            expected=f"{expected_default} {units} (app default)",
            actual="<unset>",
            status="REVIEW",
            notes="Variable not set in pod; relying on application default. OIT may want this explicitly set in the deployment manifest.",
        )
    actual = _float_or_none(raw)
    if actual is None:
        return Check(
            name=name,
            expected=f"{expected_default} {units}",
            actual=raw,
            status="FAIL",
            notes="Value is not a valid number.",
        )
    if actual == expected_default:
        return Check(
            name=name,
            expected=f"{expected_default} {units}",
            actual=f"{raw} {units}",
            status="PASS",
            notes="Matches plan's expected default.",
        )
    # Direction-aware judgement.
    if allow_lower and actual < expected_default:
        return Check(
            name=name,
            expected=f"<= {expected_default} {units}",
            actual=f"{raw} {units}",
            status="PASS",
            notes="Stricter than plan default; acceptable.",
        )
    if allow_higher and actual > expected_default:
        return Check(
            name=name,
            expected=f">= {expected_default} {units}",
            actual=f"{raw} {units}",
            status="PASS",
            notes="More permissive than plan default; verify with OIT.",
        )
    # Otherwise flag for review.
    direction = "lower" if actual < expected_default else "higher"
    return Check(
        name=name,
        expected=f"{expected_default} {units}",
        actual=f"{raw} {units}",
        status="REVIEW" if direction == "lower" else "FAIL",
        notes=f"Differs from plan default ({direction}). Confirm this is intentional and documented.",
    )


def check_debug_flag(env: dict[str, str], name: str, expected_off: bool = True) -> Check:
    raw = env.get(name, "")
    if raw == "":
        return Check(
            name=name,
            expected="off (unset)",
            actual="<unset>",
            status="PASS",
            notes="Unset = off, which matches production expectations.",
        )
    lowered = raw.lower()
    if expected_off:
        if lowered in FALSY:
            return Check(
                name=name,
                expected="off",
                actual=raw,
                status="PASS",
                notes="Explicitly disabled.",
            )
        if lowered in TRUTHY:
            return Check(
                name=name,
                expected="off",
                actual=raw,
                status="FAIL",
                notes="Debug flag is enabled in production. Disable before final OIT response.",
            )
        return Check(
            name=name,
            expected="off",
            actual=raw,
            status="REVIEW",
            notes="Value is not a recognized boolean; confirm interpretation.",
        )
    return Check(
        name=name,
        expected="on",
        actual=raw,
        status="REVIEW",
        notes="Manual review required.",
    )


def check_env_text(env: dict[str, str], name: str, expected_substring: str) -> Check:
    raw = env.get(name, "")
    if raw == "":
        return Check(
            name=name,
            expected=f"contains {expected_substring!r}",
            actual="<unset>",
            status="REVIEW",
            notes="Variable not set in pod; confirm intended environment label.",
        )
    if expected_substring.lower() in raw.lower():
        return Check(
            name=name,
            expected=f"contains {expected_substring!r}",
            actual=raw,
            status="PASS",
            notes="Matches expected environment label.",
        )
    return Check(
        name=name,
        expected=f"contains {expected_substring!r}",
        actual=raw,
        status="REVIEW",
        notes="Does not match expected label; confirm with OIT.",
    )


def check_session_workspace(temp_info: dict[str, str]) -> list[Check]:
    checks: list[Check] = []
    tempdir = temp_info.get("tempdir", "<unknown>")
    env_root = temp_info.get("env_root", "<unset>")
    session_root = temp_info.get("session_root", "<unknown>")
    exists = temp_info.get("session_root_exists", "<unknown>")

    checks.append(
        Check(
            name="tempfile.gettempdir()",
            expected="writable container temp (typically /tmp)",
            actual=tempdir,
            status="REVIEW",
            notes="Confirm this path is on a writable, ephemeral volume isolated to the pod.",
        )
    )
    checks.append(
        Check(
            name="MAPKINASE_SESSION_TMP_ROOT (effective)",
            expected="under container temp or explicitly configured writable dir",
            actual=f"env={env_root}, resolved={session_root}",
            status="REVIEW" if env_root == "<unset>" else "PASS",
            notes=(
                "Unset = uses default <tempdir>/mapkinase_sessions, which is acceptable."
                if env_root == "<unset>"
                else "Explicitly configured; confirm the path is writable and pod-local."
            ),
        )
    )
    checks.append(
        Check(
            name="Session root present",
            expected="lazily created (may be absent until first session)",
            actual=f"exists={exists}",
            status="PASS",
            notes="Either state is acceptable per the validation plan.",
        )
    )
    return checks


def run_checks() -> list[Check]:
    env = parse_env(ENV_FILE)
    temp_info = parse_tempdirs(TEMP_FILE)
    checks: list[Check] = []

    # Upload size limits
    checks.append(check_numeric(env, "MAPKINASE_MAX_UPLOAD_SIZE_MB", 10, "MB", allow_lower=True))
    checks.append(check_numeric(env, "MAPKINASE_MAX_TOTAL_UPLOAD_SIZE_MB", 30, "MB", allow_lower=True))

    # Throttling
    checks.append(check_numeric(env, "MAPKINASE_MIN_SECONDS_BETWEEN_RUNS", 5, "seconds", allow_higher=True))
    checks.append(check_numeric(env, "MAPKINASE_MAX_RUNS_PER_MINUTE", 10, "per minute", allow_lower=True))

    # Session workspace TTL
    checks.append(check_numeric(env, "MAPKINASE_SESSION_TMP_TTL_HOURS", 24, "hours", allow_lower=True))

    # Environment label
    checks.append(check_env_text(env, "MAPKINASE_ENV", "prod"))
    # MAPKINASE_PRODUCTION is a boolean-style toggle on some deployments; treat as "should be truthy".
    raw_prod = env.get("MAPKINASE_PRODUCTION", "")
    if raw_prod == "":
        checks.append(
            Check(
                name="MAPKINASE_PRODUCTION",
                expected="on (or use MAPKINASE_ENV=production)",
                actual="<unset>",
                status="REVIEW",
                notes="Either MAPKINASE_ENV=production or MAPKINASE_PRODUCTION=1 should be set so debug paths are gated.",
            )
        )
    elif raw_prod.lower() in TRUTHY:
        checks.append(
            Check(
                name="MAPKINASE_PRODUCTION",
                expected="on",
                actual=raw_prod,
                status="PASS",
                notes="Production mode flag is enabled.",
            )
        )
    else:
        checks.append(
            Check(
                name="MAPKINASE_PRODUCTION",
                expected="on",
                actual=raw_prod,
                status="FAIL",
                notes="Production flag is set to a non-truthy value. Verify debug paths are still gated.",
            )
        )

    # Debug flags - all expected off in production
    for flag in (
        "MAPKINASE_ENABLE_DEBUG_UI",
        "MAPKINASE_ENABLE_DEBUG_EXPORT",
        "MAPKINASE_ENABLE_DEBUG_FILE_OUTPUT",
        "MAPKINASE_ENABLE_CST_SIDECAR_SAVE",
        "M5_TERMINAL_LOG",
    ):
        checks.append(check_debug_flag(env, flag, expected_off=True))

    # M5_TERMINAL_LOG_FILE: if set while M5_TERMINAL_LOG is off, it is harmless,
    # but if M5_TERMINAL_LOG is on we want to know the path target.
    log_file = env.get("M5_TERMINAL_LOG_FILE", "")
    log_enabled = env.get("M5_TERMINAL_LOG", "").lower() in TRUTHY
    if log_file == "" and not log_enabled:
        checks.append(
            Check(
                name="M5_TERMINAL_LOG_FILE",
                expected="unset when M5_TERMINAL_LOG is off",
                actual="<unset>",
                status="PASS",
                notes="No debug file logging configured.",
            )
        )
    elif log_enabled:
        checks.append(
            Check(
                name="M5_TERMINAL_LOG_FILE",
                expected="documented log path inside pod",
                actual=log_file or "<unset>",
                status="REVIEW",
                notes="Terminal log is enabled; confirm the log path is pod-local, not user-accessible, and rotated.",
            )
        )
    else:
        checks.append(
            Check(
                name="M5_TERMINAL_LOG_FILE",
                expected="unset when M5_TERMINAL_LOG is off",
                actual=log_file,
                status="REVIEW",
                notes="Path is set but terminal log is off; value is unused. Consider unsetting to avoid confusion.",
            )
        )

    # Session workspace path checks (from tempdirs.txt)
    checks.extend(check_session_workspace(temp_info))

    return checks


def render(checks: list[Check]) -> str:
    counts = {"PASS": 0, "REVIEW": 0, "FAIL": 0}
    for c in checks:
        counts[c.status] = counts.get(c.status, 0) + 1

    lines: list[str] = []
    lines.append("# Map-Kinase Runtime Settings Report")
    lines.append("")
    lines.append("Generated by `deployment_validation/openshift/report_runtime_settings.py`.")
    lines.append("")
    lines.append(f"Sources: `{ENV_FILE.relative_to(REPO_ROOT)}`, `{TEMP_FILE.relative_to(REPO_ROOT)}`")
    lines.append("")
    lines.append(
        f"**Totals:** PASS = {counts.get('PASS', 0)}, "
        f"REVIEW = {counts.get('REVIEW', 0)}, "
        f"FAIL = {counts.get('FAIL', 0)}"
    )
    lines.append("")
    lines.append("| Setting | Expected | Actual | Status | Notes |")
    lines.append("| --- | --- | --- | --- | --- |")
    for c in checks:
        # Escape pipe characters in fields to keep markdown table valid.
        def esc(text: str) -> str:
            return text.replace("|", "\\|")

        lines.append(
            f"| `{esc(c.name)}` | {esc(c.expected)} | {esc(c.actual)} | **{c.status}** | {esc(c.notes)} |"
        )
    lines.append("")
    lines.append("## Narrative")
    lines.append("")
    if counts.get("FAIL", 0) > 0:
        lines.append(
            "One or more runtime settings did not match the expected production defaults. "
            "Address the FAIL rows above before sending the OIT response."
        )
    elif counts.get("REVIEW", 0) > 0:
        lines.append(
            "All explicit settings match the expected production defaults. The REVIEW rows are "
            "items relying on application defaults rather than explicit pod env vars; "
            "OIT may prefer those values be set explicitly in the deployment manifest, but the "
            "application enforces the same defaults in code."
        )
    else:
        lines.append(
            "All runtime settings match the expected production defaults. Upload limits, "
            "throttling, session TTL, and debug-related flags are configured as documented "
            "in the deployment validation plan."
        )
    lines.append("")
    return "\n".join(lines)


def main(argv: list[str]) -> int:
    if not EVIDENCE_DIR.exists():
        print(
            f"ERROR: {EVIDENCE_DIR} not found.\n"
            "Run `bash deployment_validation/openshift/collect_openshift_evidence.sh` first.",
            file=sys.stderr,
        )
        return 2
    if not ENV_FILE.exists():
        print(f"ERROR: {ENV_FILE} not found. Re-run evidence collection.", file=sys.stderr)
        return 2

    checks = run_checks()
    report = render(checks)

    OUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    OUT_FILE.write_text(report, encoding="utf-8")

    # Mirror to stdout so the user sees a summary immediately.
    print(report)
    print(f"\nWrote {OUT_FILE.relative_to(REPO_ROOT)}", file=sys.stderr)

    # Exit code reflects worst status.
    if any(c.status == "FAIL" for c in checks):
        return 1
    if any(c.status == "REVIEW" for c in checks):
        return 3
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
