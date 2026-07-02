"""Environment-backed configuration for public adapters.

No production secrets are bundled in this repository. Optional API keys are
read from environment variables at runtime.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from functools import lru_cache


def _env_int(name: str, default: int) -> int:
    value = os.getenv(name)
    if not value:
        return default
    try:
        return int(value)
    except ValueError:
        return default


def _env_float(name: str, default: float) -> float:
    value = os.getenv(name)
    if not value:
        return default
    try:
        return float(value)
    except ValueError:
        return default


@dataclass(frozen=True)
class CoreSettings:
    """Settings required by public source adapters."""

    ncbi_email: str = os.getenv("BIOSCOUTER_NCBI_EMAIL", os.getenv("NCBI_EMAIL", "anonymous@example.com"))
    ncbi_api_key: str | None = os.getenv("BIOSCOUTER_NCBI_API_KEY", os.getenv("NCBI_API_KEY"))
    ncbi_tool_name: str = os.getenv("BIOSCOUTER_NCBI_TOOL_NAME", "bioscouter-core")
    mgrast_api_key: str | None = os.getenv("BIOSCOUTER_MGRAST_API_KEY", os.getenv("MGRAST_API_KEY"))
    default_max_results: int = _env_int("BIOSCOUTER_DEFAULT_MAX_RESULTS", 50)
    source_timeout_seconds: float = _env_float("BIOSCOUTER_SOURCE_TIMEOUT_SECONDS", 30.0)
    min_relevance_score: float = _env_float("BIOSCOUTER_MIN_RELEVANCE_SCORE", 0.0)


@lru_cache(maxsize=1)
def get_settings() -> CoreSettings:
    """Return cached public-core settings."""

    return CoreSettings()

