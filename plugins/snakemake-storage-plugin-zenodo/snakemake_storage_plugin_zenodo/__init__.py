# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import asyncio
import hashlib
import os
import shutil
import time
from collections.abc import Iterable
from contextlib import asynccontextmanager
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any
from urllib.parse import urlparse

import aiohttp
import platformdirs
import snakemake_storage_plugin_http as http_base
from reretry import retry
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
from tqdm_loggable.auto import tqdm

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


@dataclass
class ZenodoFileMetadata:
    """Metadata for a file in a Zenodo record."""

    md5sum: str
    size: int


class WrongChecksum(Exception):
    def __init__(self, observed: str, expected: str):
        self.observed = observed
        self.expected = expected
        super().__init__(f"Checksum mismatch: expected {expected}, got {observed}")


retry_decorator = retry(
    exceptions=(aiohttp.ClientError, asyncio.TimeoutError, OSError, WrongChecksum),
    tries=3,
    delay=3,
    backoff=2,
    logger=get_logger(),
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

        # Cache for Zenodo record metadata to avoid repeated API calls
        self._record_cache: dict[str, dict[str, ZenodoFileMetadata]] = {}

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

    @retry_decorator
    async def get_metadata(
        self, record_id: str, netloc: str
    ) -> dict[str, ZenodoFileMetadata]:
        """
        Retrieve and cache file metadata for a Zenodo record.

        Args:
            record_id: The Zenodo record ID
            netloc: Network location (e.g., "zenodo.org")

        Returns:
            Dictionary mapping filename to ZenodoFileMetadata
        """
        # Check cache first
        if record_id in self._record_cache:
            return self._record_cache[record_id]

        # Fetch from API
        api_url = f"https://{netloc}/api/records/{record_id}"

        async with self.httpr("get", api_url) as response:
            if response.status != 200:
                raise WorkflowError(
                    f"Failed to fetch Zenodo record metadata: HTTP {response.status}"
                )

            data = await response.json()

            # Parse files array and build metadata dict
            files_metadata = {}
            files = data.get("files", [])
            for file_info in files:
                filename = file_info.get("key")
                checksum_str = file_info.get("checksum", "")
                size = file_info.get("size", 0)

                if not filename:
                    continue

                # Parse MD5 checksum
                if checksum_str.startswith("md5:"):
                    md5sum = checksum_str[4:].lower()
                    files_metadata[filename] = ZenodoFileMetadata(
                        md5sum=md5sum, size=size
                    )

            # Store in cache
            self._record_cache[record_id] = files_metadata

            return files_metadata

    @asynccontextmanager
    async def httpr(self, method, url):
        async with (
            self._download_semaphore,
            self.session() as session,
            session.request(method, url, allow_redirects=True) as response,
        ):
            wait_time = self._get_rate_limit_wait_time(response.headers)
            if wait_time is not None:
                await asyncio.sleep(wait_time)
                raise aiohttp.ClientError("Rate limit exceeded, retrying after wait")

            yield response


# Implementation of storage object
class StorageObject(StorageObjectRead):
    def __post_init__(self):
        super().__post_init__()
        self.query_path: Path = self.provider.cache_dir / self.local_suffix()
        self.query_path.parent.mkdir(exist_ok=True, parents=True)

        # Parse URL to extract record ID and filename
        # URL format: https://zenodo.org/records/{record_id}/files/{filename}
        parsed = urlparse(str(self.query))
        _records, record_id, _files, filename = parsed.path.strip("/").split(
            "/", maxsplit=3
        )

        if _records != "records" or _files != "files":
            raise WorkflowError(
                f"Invalid Zenodo URL format: {self.query}. "
                "Expected format: https://zenodo.org/records/{{record_id}}/files/{{filename}}"
            )

        self.record_id = record_id
        self.filename = filename
        self.netloc = parsed.netloc

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

    async def managed_exists(self) -> bool:
        exists = self.query_path.exists()
        if exists:
            return True

        metadata = await self.provider.get_metadata(self.record_id, self.netloc)
        return self.filename in metadata

    async def managed_mtime(self) -> float:
        return 0

    async def managed_size(self) -> int:
        exists = self.query_path.exists()
        if exists:
            return self.query_path.stat().st_size

        metadata = await self.provider.get_metadata(self.record_id, self.netloc)
        return metadata[self.filename].size if self.filename in metadata else 0

    async def inventory(self, cache: IOCacheStorageInterface):
        """
        Gather file metadata (existence, size) from cache or remote.
        Checks local cache first, then queries remote if needed.
        """
        key = self.cache_key()
        if key in cache.exists_in_storage:
            # Already inventorized
            return

        exists = self.query_path.exists()
        if exists:
            cache.exists_in_storage[key] = exists
            cache.mtime[key] = 0
            cache.size[key] = self.query_path.stat().st_size
            return

        metadata = await self.provider.get_metadata(self.record_id, self.netloc)
        exists = self.filename in metadata
        cache.exists_in_storage[key] = exists
        cache.mtime[key] = Mtime(storage=0)
        cache.size[key] = metadata[self.filename].size if exists else 0

    def cleanup(self):
        """Nothing to cleanup"""
        pass

    def exists(self):
        raise NotImplementedError()

    def size(self):
        raise NotImplementedError()

    def mtime(self):
        raise NotImplementedError()

    def retrieve_object(self):
        return NotImplementedError()

    async def verify_checksum(self, path: Path):
        """
        Fetch the MD5 checksum for this file from the Zenodo API.

        Uses the provider's caching mechanism to avoid repeated API calls
        for files from the same record.

        Raises:
            WrongChecksum
        """
        # Get cached or fetch record metadata
        metadata = await self.provider.get_metadata(self.record_id, self.netloc)
        if self.filename not in metadata:
            raise WorkflowError(
                f"File {self.filename} not found in Zenodo record {self.record_id}"
            )

        checksum_expected = metadata[self.filename].md5sum
        with open(path, "rb") as f:
            checksum_observed = hashlib.file_digest(f, "md5").hexdigest().lower()

        if checksum_expected != checksum_observed:
            raise WrongChecksum(observed=checksum_observed, expected=checksum_expected)

    @retry_decorator
    async def managed_retrieve(self):
        """Async download with concurrency control and progress bar"""
        local_path = self.local_path()
        local_path.parent.mkdir(parents=True, exist_ok=True)

        # If already in cache, verify checksum and copy
        if self.query_path.exists():
            # Verify cached file checksum
            logger.info(
                f"Retrieved {self.filename} of zenodo record {self.record_id} from cache"
            )
            shutil.copy2(self.query_path, local_path)
            return

        try:
            # Download from Zenodo using a get request, rate limit errors are detected and
            # raise WorkflowError to trigger a retry
            async with self.provider.httpr("get", str(self.query)) as response:
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

            await self.verify_checksum(local_path)

            # Copy to cache after successful verification
            self.query_path.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(local_path, self.query_path)
            logger.info(f"Cached {self.query_path.name} to {self.provider.cache_dir}")

        except (aiohttp.ClientError, asyncio.TimeoutError, OSError, WrongChecksum):
            # Clean up partial downloads on network or file system errors
            if local_path.exists():
                local_path.unlink()
            raise
