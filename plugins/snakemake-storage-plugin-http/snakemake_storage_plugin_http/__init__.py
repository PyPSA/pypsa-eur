__author__ = "Christopher Tomkins-Tinch, Johannes Köster"
__copyright__ = "Copyright 2023, Christopher Tomkins-Tinch, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

from contextlib import contextmanager
from dataclasses import dataclass, field
import email
from functools import partial
import os
import re
import shutil
from typing import Any, Iterable, List, Optional
import requests
from requests.auth import AuthBase, HTTPBasicAuth, HTTPDigestAuth
from requests_oauthlib import OAuth1
from urllib.parse import urlparse

from snakemake_interface_common.logging import get_logger
from snakemake_interface_storage_plugins.settings import StorageProviderSettingsBase
from snakemake_interface_storage_plugins.storage_provider import (
    StorageProviderBase,
    StorageQueryValidationResult,
    ExampleQuery,
    QueryType,
)
from snakemake_interface_storage_plugins.storage_object import StorageObjectRead
from snakemake_interface_storage_plugins.io import IOCacheStorageInterface, Mtime
from snakemake_interface_common.exceptions import WorkflowError
from snakemake_interface_storage_plugins.common import Operation


AUTH_METAVAR = "AUTH_TYPE=ARG1,ARG2,..."
logger = get_logger()


def parse_auth(arg):
    matches = re.match(r"^(?P<type>\w+)=(?P<arg>(\w+,?)+)$", arg)
    if not matches:
        raise WorkflowError(
            f"Authentication requires a string of the form {AUTH_METAVAR}"
        )
    auth_type = matches.group("type")
    if auth_type == "HTTPBasicAuth":
        return HTTPBasicAuth(*matches.group("arg").split(","))
    elif auth_type == "HTTPDigestAuth":
        return HTTPDigestAuth(*matches.group("arg").split(","))
    elif auth_type == "OAuth1":
        return OAuth1(*matches.group("arg").split(","))
    else:
        raise WorkflowError(
            f"Authentication type {auth_type} not supported. "
            f"Please choose one of HTTPBasicAuth, HTTPDigestAuth, or OAuth1."
        )


# Define settings for your storage plugin (e.g. host url, credentials).
# They will occur in the Snakemake CLI as --storage-<storage-plugin-name>-<param-name>
# Make sure that all defined fields are 'Optional' and specify a default value
# of None or anything else that makes sense in your case.
# Note that we allow storage plugin settings to be tagged by the user. That means,
# that each of them can be specified multiple times (an implicit nargs=+), and
# the user can add a tag in front of each value (e.g. tagname1:value1 tagname2:value2).
# This way, a storage plugin can be used multiple times within a workflow with different
# settings.
@dataclass
class StorageProviderSettings(StorageProviderSettingsBase):
    auth: Optional[AuthBase] = field(
        default=None,
        metadata={
            "help": "HTTP(S) authentication. AUTH_TYPE is the class name of "
            "requests.auth (e.g. HTTPBasicAuth), ARG1,ARG2,... are the arguments "
            "required by the specified type.",
            "metavar": AUTH_METAVAR,
            "parse_func": parse_auth,
            "env_var": True,
        },
    )
    allow_redirects: Optional[bool] = field(
        default=True,
        metadata={
            "help": "Allow redirects when retrieving files.",
        },
    )
    supports_head: Optional[bool] = field(
        default=True,
        metadata={
            "help": "Whether the storage provider supports HTTP HEAD requests.",
        },
    )


# Required:
# Implementation of your storage provider
class StorageProvider(StorageProviderBase):
    def rate_limiter_key(self, query: str, operation: Operation) -> Any:
        """Return a key for identifying a rate limiter given a query and an operation.

        This is used to identify a rate limiter for the query.
        E.g. for a storage provider like http that would be the host name.
        For s3 it might be just the endpoint URL.
        """
        parsed = urlparse(query)
        return parsed.netloc

    @classmethod
    def example_queries(cls) -> List[ExampleQuery]:
        """Return an example query with description for this storage provider."""
        return [
            ExampleQuery(
                query="https://example.com/file.txt",
                description="A file URL",
                type=QueryType.INPUT,
            )
        ]

    def default_max_requests_per_second(self) -> float:
        """Return the default maximum number of requests per second for this storage
        provider."""
        return 10.0

    def use_rate_limiter(self) -> bool:
        """Return False if no rate limiting is needed for this provider."""
        return True

    @classmethod
    def is_valid_query(cls, query: str) -> StorageQueryValidationResult:
        try:
            parsed = urlparse(query)
        except Exception as e:
            return StorageQueryValidationResult(
                query=query,
                valid=False,
                reason=f"cannot be parsed as URL ({e})",
            )
        if not (parsed.scheme == "http" or parsed.scheme == "https"):
            return StorageQueryValidationResult(
                query=query,
                valid=False,
                reason="scheme must be http or https",
            )
        return StorageQueryValidationResult(
            query=query,
            valid=True,
        )

    def list_objects(self, query: Any) -> Iterable[str]:
        raise NotImplementedError()


# Required:
# Implementation of storage object (also check out
# snakemake_interface_storage_plugins.storage_object for more base class options)
class StorageObject(StorageObjectRead):
    def local_suffix(self):
        parsed = urlparse(self.query)
        return f"{parsed.netloc}{parsed.path}"

    async def inventory(self, cache: IOCacheStorageInterface):
        """From this file, try to find as much existence and modification date
        information as possible. Only retrieve that information that comes for free
        given the current object.
        """
        name = self.local_path()
        with self.httpr(verb="HEAD") as httpr:
            res = ResponseHandler(httpr)
            cache.mtime[name] = Mtime(storage=res.mtime())
            cache.exists_in_storage[name] = res.exists()
            cache.size[name] = res.size()

    def get_inventory_parent(self) -> Optional[str]:
        """Return the parent directory of this object."""
        # this is optional and can be left as is
        return None

    def cleanup(self):
        # nothing to be done here
        pass

    def exists(self) -> bool:
        with self.httpr(verb="HEAD") as httpr:
            return ResponseHandler(httpr).exists()

    def mtime(self) -> float:
        with self.httpr(verb="HEAD") as httpr:
            return ResponseHandler(httpr).mtime()

    def size(self) -> int:
        with self.httpr(verb="HEAD") as httpr:
            return ResponseHandler(httpr).size()

    def retrieve_object(self):
        with self.httpr(stream=True) as httpr:
            # Find out if the source file is gzip compressed in order to keep
            # compression intact after the download.
            # Per default requests decompresses .gz files.
            # More details can be found here:
            # https://stackoverflow.com/questions/25749345/how-to-download-gz-files-with-requests-in-python-without-decoding-it?noredirect=1&lq=1  # noqa E501
            # Since data transferred with HTTP compression need to be decompressed
            # automatically, check the header and decode if the content is encoded.
            if (
                not self.query.endswith(".gz")
                and httpr.headers.get("Content-Encoding") == "gzip"
            ):
                # Decode non-gzipped sourcefiles automatically.
                # This is needed to decompress uncompressed files that are compressed
                # for the transfer by HTTP compression.
                httpr.raw.decode_content = True
            # if the destination path does not exist
            os.makedirs(os.path.dirname(self.local_path()), exist_ok=True)
            with open(self.local_path(), "wb") as f:
                shutil.copyfileobj(httpr.raw, f)

    @contextmanager  # makes this a context manager. after 'yield' is __exit__()
    def httpr(self, verb="GET", stream=False):
        r = None
        try:

            if verb.upper() == "GET":
                request = partial(requests.get, stream=stream)
            if verb.upper() == "HEAD":
                if self.provider.settings.supports_head:
                    request = requests.head
                else:
                    request = requests.get

            r = request(
                self.query,
                stream=stream,
                auth=self.provider.settings.auth,
                allow_redirects=self.provider.settings.allow_redirects,
            )

            logger.debug("X-RateLimit-Remaining: %s", res.headers.get("X-RateLimit-Remaining"))
            yield r
        finally:
            if r is not None:
                r.close()


@dataclass
class ResponseHandler:
    response: requests.Response

    def exists(self):
        if self.response.status_code in range(300, 308):
            raise WorkflowError(
                "The file specified appears to have been moved (HTTP "
                f"{self.response.status_code}), check the URL or enable "
                "redirects in the http storage provider "
                "(--storage-http-allow-redirects true): "
                f"{self.response.url}"
            )
        return self.response.status_code == requests.codes.ok

    def mtime(self):
        file_mtime = get_header_item(self.response, "last-modified", default=None)
        logger.debug(f"HTTP last-modified: {file_mtime}")

        epochTime = 0

        if file_mtime is not None:
            modified_tuple = email.utils.parsedate_tz(file_mtime)
            if modified_tuple is None:
                logger.debug(
                    "HTTP last-modified not in RFC2822 format: `{}`".format(file_mtime)
                )
            else:
                epochTime = email.utils.mktime_tz(modified_tuple)

        return epochTime

    def size(self):
        content_size = int(get_header_item(self.response, "content-size", default=0))

        return content_size


def get_header_item(httpr, header_name, default):
    """
    Since HTTP header capitalization may differ, this returns
    a header value regardless of case
    """

    header_value = default
    for k, v in httpr.headers.items():
        if k.lower() == header_name:
            header_value = v
    return header_value
