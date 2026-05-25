import os


_TRUE_VALUES = {"1", "true", "yes", "on"}


def env_truthy(name: str, default: bool = False) -> bool:
    raw = os.environ.get(name)
    if raw is None:
        return default
    return str(raw).strip().lower() in _TRUE_VALUES


def runtime_environment() -> str:
    raw = str(os.environ.get("MAPKINASE_ENV", "production") or "").strip().lower()
    if raw in {"prod", "production"}:
        return "production"
    if raw in {"dev", "development"}:
        return "development"
    if raw in {"test", "testing"}:
        return "testing"
    return raw or "production"


def is_production_mode() -> bool:
    if "MAPKINASE_PRODUCTION" in os.environ:
        return env_truthy("MAPKINASE_PRODUCTION", False)
    return runtime_environment() == "production"
