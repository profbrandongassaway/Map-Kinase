import logging
import os
from pathlib import Path
from typing import Iterable, List, Optional, Tuple


def _int_env(name: str, default: int, minimum: int = 1) -> int:
    raw = os.getenv(name, str(default)).strip()
    try:
        value = int(raw)
    except ValueError:
        value = default
    return value if value >= minimum else minimum


MAX_UPLOAD_SIZE_MB = _int_env("MAPKINASE_MAX_UPLOAD_SIZE_MB", 10, minimum=1)
MAX_TOTAL_UPLOAD_SIZE_MB = _int_env("MAPKINASE_MAX_TOTAL_UPLOAD_SIZE_MB", 30, minimum=1)
MIN_SECONDS_BETWEEN_RUNS = _int_env("MAPKINASE_MIN_SECONDS_BETWEEN_RUNS", 5, minimum=1)
MAX_RUNS_PER_MINUTE = _int_env("MAPKINASE_MAX_RUNS_PER_MINUTE", 10, minimum=1)

MAX_UPLOAD_SIZE_BYTES = MAX_UPLOAD_SIZE_MB * 1024 * 1024
MAX_TOTAL_UPLOAD_SIZE_BYTES = MAX_TOTAL_UPLOAD_SIZE_MB * 1024 * 1024

_RUN_GUARD_LOGGER = logging.getLogger("mapkinase.run_guard")
if not _RUN_GUARD_LOGGER.handlers:
    _handler = logging.StreamHandler()
    _handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(name)s - %(message)s"))
    _RUN_GUARD_LOGGER.addHandler(_handler)
_RUN_GUARD_LOGGER.setLevel(logging.INFO)


def evaluate_run_throttle(
    *,
    now_monotonic: float,
    last_accepted_run_monotonic: Optional[float],
    prior_attempt_times_monotonic: Iterable[float],
    min_seconds_between_runs: int = MIN_SECONDS_BETWEEN_RUNS,
    max_runs_per_minute: int = MAX_RUNS_PER_MINUTE,
) -> Tuple[bool, str, List[float]]:
    window_start = now_monotonic - 60.0
    recent_attempts = [t for t in prior_attempt_times_monotonic if t >= window_start]
    recent_attempts.append(now_monotonic)

    if last_accepted_run_monotonic is not None and (now_monotonic - last_accepted_run_monotonic) < float(min_seconds_between_runs):
        return False, "Please wait a few seconds before starting another analysis.", recent_attempts
    if len(recent_attempts) > max_runs_per_minute:
        return False, "Too many analysis runs were started recently. Please wait before trying again.", recent_attempts
    return True, "", recent_attempts


def total_upload_size_bytes(file_paths: Iterable[str]) -> Optional[int]:
    total = 0
    try:
        for path in file_paths:
            cleaned = str(path or "").strip()
            if not cleaned:
                continue
            total += int(Path(cleaned).stat().st_size)
    except OSError:
        return None
    return total


def validate_total_upload_size(file_paths: Iterable[str], max_total_upload_size_bytes: int = MAX_TOTAL_UPLOAD_SIZE_BYTES) -> Tuple[bool, Optional[int]]:
    total = total_upload_size_bytes(file_paths)
    if total is None:
        return False, None
    return total <= int(max_total_upload_size_bytes), total


def validate_total_upload_size_from_sizes(
    size_bytes: Iterable[int | float | str],
    max_total_upload_size_bytes: int = MAX_TOTAL_UPLOAD_SIZE_BYTES,
) -> Tuple[bool, Optional[int]]:
    total = 0
    try:
        for raw in size_bytes:
            if raw in (None, "", False):
                continue
            value = int(float(raw))
            if value < 0:
                return False, None
            total += value
    except (TypeError, ValueError):
        return False, None
    return total <= int(max_total_upload_size_bytes), total


def log_run_guard_event(*, session_id: str, reason: str, total_size_bytes: Optional[int] = None, run_attempt_count: Optional[int] = None) -> None:
    _RUN_GUARD_LOGGER.warning(
        "Rejected analysis run | session=%s | reason=%s | total_size_bytes=%s | attempts_last_min=%s",
        session_id or "unknown",
        reason,
        str(total_size_bytes) if total_size_bytes is not None else "na",
        str(run_attempt_count) if run_attempt_count is not None else "na",
    )
