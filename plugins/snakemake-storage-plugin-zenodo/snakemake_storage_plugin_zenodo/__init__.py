# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import asyncio
import os
import shutil
import time
from collections.abc import Iterable
from contextlib import asynccontextmanager
from dataclasses import dataclass, field
from functools import cached_property
from pathlib import Path
from typing import Any
from urllib.parse import urlparse

import aiohttp
import platformdirs
import requests
import snakemake_storage_plugin_http as http_base
from snakemake_interface_common.exceptions import WorkflowError
from snakemake_interface_common.logging import get_logger
from snakemake_interface_storage_plugins.common import Operation
from snakemake_interface_storage_plugins.io import IOCacheStorageInterface, Mtime
from snakemake_interface_storage_plugins.settings import StorageProviderSettingsBase
from snakemake_interface_storage_plugins.storage_object import StorageObjectRead
from snakemake_interface_storage_plugins.storage_provider import (
    ExampleQuery,
    QueryType,
    StorageProviderBase,
    StorageQueryValidationResult,
)
from tqdm import tqdm

logger = get_logger()


def is_zenodo_url(url):
    parsed = urlparse(url)
    return parsed.netloc.endswith("zenodo.org") and parsed.scheme in ("http", "https")


# Patch the original HTTP StorageProvider off zenodo urls, so that there is no conflict
orig_valid_query = http_base.StorageProvider.is_valid_query
http_base.StorageProvider.is_valid_query = classmethod(
    lambda c, q: (
        StorageQueryValidationResult(
            query=q, valid=False, reason="Deactivated in favour of zenodo"
        )
        if is_zenodo_url(q)
        else orig_valid_query(q)
    )
)


# Define settings for the Zenodo storage plugin
@dataclass
class StorageProviderSettings(StorageProviderSettingsBase):
    cache: Path = field(
        default_factory=lambda: Path(
            platformdirs.user_cache_dir("snakemake-pypsa-eur", ensure_exists=False)
        ),
        metadata={
            "help": "Cache directory for downloaded Zenodo files (default: platform-dependent user cache dir)",
            "env_var": False,
        },
    )
    max_concurrent_downloads: int | None = field(
        default=3,
        metadata={
            "help": "Maximum number of concurrent Zenodo downloads",
            "env_var": False,
        },
    )


# Implementation of storage provider
class StorageProvider(StorageProviderBase):
    def __post_init__(self):
        super().__post_init__()

        # Set up cache directory
        self.cache_dir = Path(self.settings.cache)
        self.cache_dir.mkdir(exist_ok=True, parents=True)

        # Create semaphore for limiting concurrent downloads
        self._download_semaphore: asyncio.Semaphore = asyncio.Semaphore(
            int(self.settings.max_concurrent_downloads)
        )

        # Initialize session tracking
        self._session: aiohttp.ClientSession | None = None
        self._session_refcount: int = 0

    def use_rate_limiter(self) -> bool:
        """Return False if no rate limiting is needed for this provider."""
        return False

    def rate_limiter_key(self, query: str, operation: Operation) -> Any:
        raise NotImplementedError()

    def default_max_requests_per_second(self) -> float:
        raise NotImplementedError()

    @classmethod
    def example_queries(cls) -> list[ExampleQuery]:
        """Return an example query with description for this storage provider."""
        return [
            ExampleQuery(
                query="https://zenodo.org/records/17249457/files/ARDECO-SNPTD.2021.table.csv",
                description="A zenodo file URL",
                type=QueryType.INPUT,
            )
        ]

    @classmethod
    def is_valid_query(cls, query: str) -> StorageQueryValidationResult:
        """Only handle zenodo.org URLs"""
        if is_zenodo_url(query):
            return StorageQueryValidationResult(query=query, valid=True)

        return StorageQueryValidationResult(
            query=query,
            valid=False,
            reason="Not a Zenodo URL (only zenodo.org URLs are handled by this plugin)",
        )

    @classmethod
    def get_storage_object_cls(cls):
        return StorageObject

    def list_objects(self, query: Any) -> Iterable[str]:
        raise NotImplementedError()

    @asynccontextmanager
    async def session(self):
        """
        Reentrant async context manager for aiohttp.ClientSession.

        Creates a session on first entry and reuses it for nested calls.
        The session is closed only when all context managers have exited.

        Usage:
            async with (
                provider.session() as session,
                session.get(url) as response
            ):
                    ...
        """
        # Increment reference count
        self._session_refcount += 1

        # Create session on first entry
        if self._session is None:
            self._session = aiohttp.ClientSession()

        try:
            yield self._session
        finally:
            # Decrement reference count
            self._session_refcount -= 1
            # Close session when no more references
            if self._session_refcount == 0:
                await self._session.close()
                self._session = None


@dataclass
class ObjectState:
    exists: bool
    mtime: float
    size: int


# Implementation of storage object
class StorageObject(StorageObjectRead):
    def __post_init__(self):
        super().__post_init__()
        self.query_path: Path = self.provider.cache_dir / self.local_suffix()
        self.query_path.parent.mkdir(exist_ok=True, parents=True)

    def local_suffix(self):
        parsed = urlparse(self.query)
        return f"{parsed.netloc}{parsed.path}"

    def get_inventory_parent(self) -> str | None:
        """Return the parent directory of this object."""
        # this is optional and can be left as is
        return None

    def _stat(self, follow_symlinks: bool = True):
        """Get file stats, bypassing Path.stat() cache"""
        return os.stat(self.query_path, follow_symlinks=follow_symlinks)

    @cached_property
    def state(self) -> ObjectState:
        exists = self.query_path.exists()
        if exists:
            return ObjectState(exists=True, mtime=0, size=self._stat().st_size)

        # Perform HEAD request with rate limit checking
        try:
            response = requests.head(str(self.query), allow_redirects=True, timeout=30)
        except requests.RequestException as e:
            raise WorkflowError(f"Failed to query {self.query} from Zenodo", e)

        # Check rate limits
        wait_time = self._get_rate_limit_wait_time(response.headers)
        if wait_time is not None:
            time.sleep(wait_time)
            raise WorkflowError("Rate limit exceeded, retrying after wait")

        exists = response.status_code == 200
        size = int(response.headers.get("content-length", 0)) if exists else 0
        return ObjectState(exists=exists, mtime=0, size=size)

    def exists(self) -> bool:
        # return True if the object exists
        return self.state.exists

    def mtime(self) -> float:
        # return the modification time
        return self.state.mtime

    def size(self) -> int:
        # return the size in bytes
        return self.state.size

    def _get_rate_limit_wait_time(self, headers) -> float | None:
        """
        Calculate wait time based on rate limit headers.

        Returns:
            float | None: Wait time in seconds if rate limited, None otherwise
        """
        remaining = int(headers.get("X-RateLimit-Remaining", 100))
        reset_time = int(headers.get("X-RateLimit-Reset", 0))

        if remaining >= 1:
            return None

        wait_seconds = max(0, reset_time - time.time() + 1)
        logger.info(
            f"Zenodo rate limit exceeded. Waiting {wait_seconds:.0f}s until reset..."
        )
        return wait_seconds

    async def inventory(self, cache: IOCacheStorageInterface):
        """
        Gather file metadata (existence, size) from cache or remote.
        Checks local cache first, then queries remote if needed.
        """
        key = self.cache_key()

        if key in cache.exists_in_storage:
            # Already inventorized
            return

        state = self.state
        cache.exists_in_storage[key] = state.exists
        cache.mtime[key] = Mtime(storage=state.mtime)
        cache.size[key] = state.size

    def cleanup(self):
        """Nothing to cleanup"""
        pass

    def retrieve_object(self):
        return NotImplementedError()

    async def managed_retrieve(self):
        """Async download with concurrency control and progress bar"""
        local_path = self.local_path()

        local_path.parent.mkdir(parents=True, exist_ok=True)

        # If already in cache, copy from cache
        if self.query_path.exists():
            # Extract record ID and filename from URL
            parsed = urlparse(str(self.query))
            path_parts = parsed.path.strip("/").split("/")
            record_id = path_parts[1] if len(path_parts) > 1 else "unknown"
            filename = self.query_path.name

            logger.info(f"Retrieved {filename} of zenodo record {record_id} from cache")
            shutil.copy2(self.query_path, local_path)
            return

        # Download from Zenodo using a get request, rate limit errors are handled and trigger a retry
        try:
            async with (
                self.provider._download_semaphore,
                self.provider.session() as session,
                session.get(str(self.query), allow_redirects=True) as response,
            ):
                wait_time = self._get_rate_limit_wait_time(response.headers)
                if wait_time is not None:
                    await asyncio.sleep(wait_time)
                    raise aiohttp.ClientError(
                        "Rate limit exceeded, retrying after wait"
                    )

                if response.status != 200:
                    raise WorkflowError(
                        f"Failed to download from Zenodo: HTTP {response.status}"
                    )

                total_size = int(response.headers.get("content-length", 0))

                # Download to local path with progress bar
                with local_path.open(mode="wb") as f:
                    with tqdm(
                        total=total_size,
                        unit="B",
                        unit_scale=True,
                        desc=self.query_path.name,
                        position=None,
                        leave=True,
                    ) as pbar:
                        async for chunk in response.content.iter_chunked(8192):
                            f.write(chunk)
                            pbar.update(len(chunk))

                self.query_path.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(local_path, self.query_path)
                logger.info(
                    f"Cached {self.query_path.name} to {self.provider.cache_dir}"
                )
        except (aiohttp.ClientError, asyncio.TimeoutError, OSError) as e:
            # Clean up partial downloads on network or file system errors
            if local_path.exists():
                local_path.unlink()
            raise WorkflowError(f"Failed to retrieve {self.query} from Zenodo", e)
