import hashlib
import logging
import os
import shutil
import tempfile
import time
import uuid
from pathlib import Path
from typing import Any, Optional


_WORKSPACE_LOGGER = logging.getLogger("mapkinase.session_workspace")
if not _WORKSPACE_LOGGER.handlers:
    _handler = logging.StreamHandler()
    _handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(name)s - %(message)s"))
    _WORKSPACE_LOGGER.addHandler(_handler)
_WORKSPACE_LOGGER.setLevel(logging.INFO)


_FORBIDDEN_ROOT_TOKENS = {"www", "static", "public", "jsonfiles", "cache", "stored_pathways", "output"}


def _int_env(name: str, default: int, minimum: int = 1) -> int:
    raw = os.getenv(name, str(default)).strip()
    try:
        value = int(raw)
    except ValueError:
        value = default
    return value if value >= minimum else minimum


def _default_workspace_root() -> Path:
    return Path(tempfile.gettempdir()) / "mapkinase_sessions"


def _configured_workspace_root() -> Path:
    raw = str(os.getenv("MAPKINASE_SESSION_TMP_ROOT", "") or "").strip()
    return Path(raw).expanduser() if raw else _default_workspace_root()


def _root_is_forbidden(path: Path) -> bool:
    lowered = {part.lower() for part in path.parts}
    return bool(lowered & _FORBIDDEN_ROOT_TOKENS)


def _ensure_workspace_root() -> Path:
    configured = _configured_workspace_root()
    resolved = configured.resolve(strict=False)
    if _root_is_forbidden(resolved):
        fallback = _default_workspace_root().resolve(strict=False)
        _WORKSPACE_LOGGER.warning(
            "Configured MAPKINASE_SESSION_TMP_ROOT points to a forbidden area. Falling back to the default temp workspace root.",
        )
        resolved = fallback
    resolved.mkdir(parents=True, exist_ok=True)
    return resolved


SESSION_TMP_ROOT = _ensure_workspace_root()
SESSION_TMP_TTL_HOURS = _int_env("MAPKINASE_SESSION_TMP_TTL_HOURS", 24, minimum=1)


def _session_token(session_id: Optional[str] = None) -> str:
    sid = str(session_id or "").strip()
    sid_hash = hashlib.sha256(sid.encode("utf-8")).hexdigest()[:12] if sid else "anon"
    return f"{sid_hash}_{uuid.uuid4().hex}"


def create_session_workspace(session_id: Optional[str] = None) -> Path:
    root = _ensure_workspace_root()
    workspace = root / _session_token(session_id)
    workspace.mkdir(parents=True, exist_ok=False)
    _WORKSPACE_LOGGER.info("Created session workspace | workspace_id=%s", workspace.name)
    return workspace


def _session_attr_name() -> str:
    return "_mapkinase_session_workspace"


def _get_existing_session_workspace(session: Any) -> Optional[Path]:
    existing = getattr(session, _session_attr_name(), None)
    if not existing:
        return None
    try:
        workspace = Path(str(existing)).resolve(strict=False)
    except Exception:
        return None
    return workspace


def get_session_workspace(session: Any) -> Path:
    existing = _get_existing_session_workspace(session)
    if existing and existing.exists():
        return existing

    raw_session_id = str(getattr(session, "id", "") or getattr(session, "_session_id", "") or "").strip()
    workspace = create_session_workspace(raw_session_id or None)
    try:
        setattr(session, _session_attr_name(), str(workspace))
    except Exception:
        pass
    return workspace


def safe_session_path(session: Any, filename_or_subpath: str) -> Path:
    workspace = get_session_workspace(session).resolve(strict=False)
    relative = Path(str(filename_or_subpath or "").strip())
    if relative.is_absolute():
        raise ValueError("Session-relative paths cannot be absolute.")
    candidate = (workspace / relative).resolve(strict=False)
    if workspace not in candidate.parents and candidate != workspace:
        raise ValueError("Session path escapes the workspace root.")
    candidate.parent.mkdir(parents=True, exist_ok=True)
    return candidate


def cleanup_session_workspace(session: Any) -> bool:
    workspace = _get_existing_session_workspace(session)
    if workspace is None:
        return False

    root = _ensure_workspace_root().resolve(strict=False)
    target = workspace.resolve(strict=False)
    if root not in target.parents:
        _WORKSPACE_LOGGER.warning(
            "Refusing to remove workspace outside root | target_name=%s",
            target.name,
        )
        return False

    if not target.exists():
        return True

    try:
        shutil.rmtree(target, ignore_errors=False)
        _WORKSPACE_LOGGER.info("Deleted session workspace | workspace_id=%s", target.name)
        return True
    except Exception as exc:
        _WORKSPACE_LOGGER.warning(
            "Failed to delete session workspace | workspace_id=%s | error=%s",
            target.name,
            str(exc),
        )
        return False


def cleanup_old_session_workspaces(max_age_hours: int = SESSION_TMP_TTL_HOURS) -> int:
    root = _ensure_workspace_root().resolve(strict=False)
    cutoff = time.time() - (max(1, int(max_age_hours)) * 3600)
    removed = 0
    for child in root.iterdir():
        if not child.is_dir():
            continue
        try:
            modified = child.stat().st_mtime
        except OSError:
            continue
        if modified >= cutoff:
            continue
        try:
            shutil.rmtree(child, ignore_errors=False)
            removed += 1
            _WORKSPACE_LOGGER.info(
                "Removed stale session workspace | workspace_id=%s | max_age_hours=%s",
                child.name,
                int(max_age_hours),
            )
        except Exception as exc:
            _WORKSPACE_LOGGER.warning(
                "Failed to remove stale session workspace | workspace_id=%s | error=%s",
                child.name,
                str(exc),
            )
    return removed


def safely_delete_temp_file(path: str | Path, *, reason: str = "") -> bool:
    cleaned = str(path or "").strip()
    if not cleaned:
        return False
    target = Path(cleaned).resolve(strict=False)
    if not target.exists() or not target.is_file():
        return False

    temp_root = Path(tempfile.gettempdir()).resolve(strict=False)
    session_root = _ensure_workspace_root().resolve(strict=False)
    in_temp_root = temp_root in target.parents
    in_session_root = session_root in target.parents
    if not (in_temp_root or in_session_root):
        _WORKSPACE_LOGGER.warning(
            "Refusing to delete file outside allowed temp roots | filename=%s | reason=%s",
            target.name,
            reason or "unspecified",
        )
        return False

    try:
        target.unlink(missing_ok=True)
        _WORKSPACE_LOGGER.info(
            "Deleted temporary file | file=%s | reason=%s",
            str(target.name),
            reason or "unspecified",
        )
        return True
    except Exception as exc:
        _WORKSPACE_LOGGER.warning(
            "Failed to delete temporary file | file=%s | reason=%s | error=%s",
            str(target.name),
            reason or "unspecified",
            str(exc),
        )
        return False
