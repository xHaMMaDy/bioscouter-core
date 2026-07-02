"""
Base adapter interface for every data source adapter.

Each source (GEO, PRIDE, ENCODE, ...) implements this interface and the base
class provides shared infrastructure:

* HTTP retry helper (:meth:`BaseSourceAdapter._fetch_with_retry`) so every
  adapter doesn't reimplement the same exponential-backoff loop.
* Field normalizers (organisms, lists, dates).
* Shared keyword-relevance scorer with a configurable stopword set.
* Default ``is_available()`` that does **not** issue a live search —
  subclasses override only when a real probe is cheap.

The retry helper accepts an httpx client passed in by the subclass so we
don't presume a particular client lifetime; subclasses keep ownership of
their pool. Transient failures (timeout, network, 5xx) are retried with
exponential backoff capped at ``RETRY_BASE_DELAY * 2^(MAX_RETRY_ATTEMPTS-1)``.
"""

from __future__ import annotations

import asyncio
import re
from abc import ABC, abstractmethod
from datetime import date, datetime
from typing import Any, AsyncIterator, List, Optional, Set

import httpx
import structlog

from ..models.unified import (
    DataSource,
    OmicsType,
    SOURCE_REGISTRY,
    SourceInfo,
    UnifiedDataset,
    UnifiedSearchStatus,
)


# === SHARED CONSTANTS ===

# Common stopwords for relevance scoring (shared across all adapters).
COMMON_STOPWORDS: Set[str] = frozenset({
    "the", "a", "an", "and", "or", "but", "in", "on", "at", "to", "for",
    "of", "with", "by", "from", "as", "is", "was", "are", "were", "been",
    "be", "have", "has", "had", "do", "does", "did", "will", "would",
    "could", "should", "may", "might", "can", "its", "all", "any", "some",
    "more", "most", "than", "into", "data", "dataset", "datasets",
    "analysis", "study", "studies", "samples", "using", "based",
    "large", "small", "high", "low", "show", "shows", "showing",
    "find", "finds", "finding", "look", "search", "searching",
    "get", "want", "need", "looking", "research", "experiment", "experiments",
})

# Default timeout for HTTP requests (seconds).
DEFAULT_HTTP_TIMEOUT = 30.0

# Maximum retry attempts for transient failures.
MAX_RETRY_ATTEMPTS = 3

# Base delay between retries; actual delay is RETRY_BASE_DELAY * 2**attempt.
RETRY_BASE_DELAY = 1.0

# Default per-adapter httpx connection limits. Each adapter holds its own
# client; without this cap, 14 adapters × httpx default 100 = 1400 sockets.
# Override at construction time if a specific source needs more headroom.
DEFAULT_KEEPALIVE_CONNECTIONS = 10
DEFAULT_MAX_CONNECTIONS = 20


def build_default_limits(
    keepalive: int = DEFAULT_KEEPALIVE_CONNECTIONS,
    max_connections: int = DEFAULT_MAX_CONNECTIONS,
) -> httpx.Limits:
    """Return an :class:`httpx.Limits` for adapter HTTP clients.

    Centralizing the limits lets the operator tune them once, and forces a
    sensible cap on every adapter that adopts it.
    """
    return httpx.Limits(
        max_keepalive_connections=keepalive,
        max_connections=max_connections,
    )


class BaseSourceAdapter(ABC):
    """Abstract base class for all omics data source adapters.

    Each adapter is responsible for:

    * Connecting to a specific data source API.
    * Translating search queries to source-specific formats.
    * Normalizing results to :class:`UnifiedDataset`.
    * Handling pagination, rate limiting, and retries.

    Use :meth:`_fetch_with_retry` for the retry loop instead of writing
    your own — adapters that don't need source-specific logic should not
    own a copy of the same exponential-backoff code.
    """

    def __init__(self):
        self.logger = structlog.get_logger(self.__class__.__name__)

    # ----------------------------------------------------- abstract surface

    @property
    @abstractmethod
    def source(self) -> DataSource:
        """The data source this adapter handles."""

    @property
    @abstractmethod
    def supported_omics(self) -> List[OmicsType]:
        """List of omics types this source supports."""

    @abstractmethod
    async def search(
        self,
        query: str,
        max_results: int = 50,
        organism: Optional[str] = None,
        **kwargs,
    ) -> List[UnifiedDataset]:
        """Search the data source and return normalized results."""

    @abstractmethod
    async def get_dataset(self, accession: str) -> Optional[UnifiedDataset]:
        """Fetch a single dataset by its native accession."""

    # ----------------------------------------------------- info / utilities

    @property
    def source_info(self) -> Optional[SourceInfo]:
        """Source info from the registry, if registered."""
        return SOURCE_REGISTRY.get(self.source)

    @property
    def source_name(self) -> str:
        """Human-readable name of the source."""
        info = self.source_info
        return info.name if info else self.source.value.upper()

    @property
    def source_icon(self) -> str:
        """Emoji icon for the source."""
        info = self.source_info
        return info.icon if info else "🧬"

    def build_unified_id(self, accession: str) -> str:
        """Build the unified ID in format ``{source}:{accession}``."""
        return f"{self.source.value}:{accession}"

    def parse_accession_from_unified_id(self, unified_id: str) -> str:
        """Extract the native accession from a unified ID."""
        if ":" in unified_id:
            return unified_id.split(":", 1)[1]
        return unified_id

    def _create_source_url(self, accession: str) -> str:
        """Create the direct URL to the dataset page (override per source)."""
        info = self.source_info
        return info.url if info else ""

    # ------------------------------------------------------- streaming hook

    async def search_streaming(
        self,
        query: str,
        max_results: int = 50,
        organism: Optional[str] = None,
        **kwargs,
    ) -> AsyncIterator[UnifiedSearchStatus | UnifiedDataset]:
        """Default streaming implementation.

        This is a *fake stream*: the underlying ``search()`` call blocks
        until completion before any results are yielded. Override in
        subclasses whose source actually streams. Yields a "searching"
        status, the results, and a "complete" status.
        """
        yield UnifiedSearchStatus(
            stage="searching",
            message=f"Searching {self.source_name}...",
            progress=0.0,
            current_source=self.source,
        )
        results = await self.search(query, max_results, organism, **kwargs)
        yield UnifiedSearchStatus(
            stage="complete",
            message=f"Found {len(results)} results from {self.source_name}",
            progress=1.0,
            current_source=self.source,
            total_found=len(results),
        )
        for result in results:
            yield result

    async def is_available(self) -> bool:
        """Cheap availability check.

        Default returns ``True`` — issuing a real ``search("test", ...)`` on
        every probe (the previous behavior) burns a network round-trip and
        leaves "test" queries in upstream provider analytics. Subclasses
        should override only when the source exposes a fast HEAD/ping
        endpoint.
        """
        return True

    async def close(self) -> None:
        """Release any underlying connection pool. Default no-op."""
        return None

    # ---------------------------------------------------------- retry helper

    async def _fetch_with_retry(
        self,
        client: httpx.AsyncClient,
        url: str,
        *,
        method: str = "GET",
        params: Optional[dict] = None,
        json: Optional[dict] = None,
        max_attempts: int = MAX_RETRY_ATTEMPTS,
        base_delay: float = RETRY_BASE_DELAY,
    ) -> httpx.Response:
        """Issue an HTTP request with exponential-backoff retries.

        Retries on:

        * :class:`httpx.TimeoutException`
        * :class:`httpx.NetworkError` (connection drops, DNS failures, ...)
        * Any 5xx status code.

        4xx responses are returned to the caller as-is — those represent a
        client error and shouldn't be retried.

        Args:
            client: The HTTP client to use; the caller owns its lifetime.
            url: Target URL.
            method: HTTP method (``GET`` or ``POST``).
            params: Optional query parameters (GET).
            json: Optional JSON body (POST).
            max_attempts: How many times to attempt before giving up.
            base_delay: Seconds; actual delay is ``base_delay * 2**attempt``.

        Returns:
            The httpx Response object on the last attempt (success or 4xx).

        Raises:
            The last :class:`httpx.TimeoutException` or
            :class:`httpx.NetworkError` if all attempts fail.
        """
        last_exc: Optional[Exception] = None
        source_label = self.source_name

        for attempt in range(max_attempts):
            try:
                if method.upper() == "POST":
                    response = await client.post(url, json=json, params=params)
                else:
                    response = await client.get(url, params=params)

                if response.status_code >= 500 and attempt < max_attempts - 1:
                    delay = base_delay * (2 ** attempt)
                    self.logger.warning(
                        "Source returned 5xx, retrying",
                        source=source_label,
                        status=response.status_code,
                        attempt=attempt + 1,
                        delay=delay,
                    )
                    await asyncio.sleep(delay)
                    continue
                return response
            except (httpx.TimeoutException, httpx.NetworkError) as exc:
                last_exc = exc
                if attempt < max_attempts - 1:
                    delay = base_delay * (2 ** attempt)
                    self.logger.warning(
                        "Source request failed, retrying",
                        source=source_label,
                        error=str(exc),
                        error_type=type(exc).__name__,
                        attempt=attempt + 1,
                        delay=delay,
                    )
                    await asyncio.sleep(delay)
                    continue
                break

        # All attempts exhausted.
        if last_exc is not None:
            raise last_exc
        raise httpx.RequestError(f"{source_label} request failed after {max_attempts} attempts")

    # ---------------------------------------------------- field normalizers

    def _normalize_organism(self, organism: Any) -> List[str]:
        """Normalize organism field to list of strings."""
        if not organism:
            return []
        if isinstance(organism, str):
            return [organism]
        if isinstance(organism, list):
            return [str(o) for o in organism if o]
        return []

    def _normalize_list(self, value: Any) -> List[str]:
        """Normalize any value to a list of non-empty strings."""
        if not value:
            return []
        if isinstance(value, str):
            return [v.strip() for v in value.split(",") if v.strip()]
        if isinstance(value, list):
            return [str(v) for v in value if v]
        return [str(value)]

    def _normalize_date(self, value: Any) -> Optional[str]:
        """Normalize a date field to an ISO ``YYYY-MM-DD`` string.

        Accepts :class:`datetime`, :class:`date`, ``str``. Numeric (Unix
        epoch) inputs are not supported — they're rare in the upstream
        APIs we consume and silently producing wrong results would be
        worse than returning ``None``.
        """
        if value is None or isinstance(value, bool):
            return None
        if isinstance(value, datetime):
            return value.strftime("%Y-%m-%d")
        if isinstance(value, date):
            return value.strftime("%Y-%m-%d")
        if isinstance(value, str):
            stripped = value.strip()
            return stripped if stripped else None
        return None

    # -------------------------------------------------------- relevance score

    def _compute_keyword_relevance(
        self,
        query: str,
        title: str,
        description: str,
        additional_stopwords: Optional[Set[str]] = None,
    ) -> float:
        """Compute keyword-overlap relevance for a result.

        Returns a value in ``[0.3, 1.0]`` (results that came back at all
        get at least 0.3 — they matched something upstream). Title matches
        weigh 60%, description matches 40%, with bonuses for exact phrase
        hits. Used as a fallback when embedding-based scoring isn't
        available; the orchestrator may overwrite this score later.
        """
        if not query:
            return 0.5

        stopwords = COMMON_STOPWORDS
        if additional_stopwords:
            stopwords = stopwords | additional_stopwords

        query_lower = query.lower()
        query_words = [w.strip() for w in re.split(r"\W+", query_lower) if len(w.strip()) > 2]
        keywords = [w for w in query_words if w not in stopwords]
        if not keywords:
            keywords = query_words[:3]
        if not keywords:
            return 0.3

        title_lower = (title or "").lower()
        desc_lower = (description or "").lower()
        title_matches = sum(1 for kw in keywords if kw in title_lower)
        desc_matches = sum(1 for kw in keywords if kw in desc_lower)

        total = len(keywords)
        title_score = title_matches / total
        desc_score = desc_matches / total
        base_score = (title_score * 0.6) + (desc_score * 0.4)

        # Bonuses for exact phrase match — these are honest signals.
        if query_lower in title_lower:
            base_score = min(1.0, base_score + 0.3)
        elif query_lower in desc_lower:
            base_score = min(1.0, base_score + 0.15)

        # Floor of 0.3: results returned by the upstream API matched
        # *something*, so don't show 0.0 even when our keywords miss.
        return max(0.3, min(1.0, base_score))

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}(source={self.source.value})>"
