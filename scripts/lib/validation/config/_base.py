# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""Base classes for config validation models."""

from collections.abc import Iterator
from typing import Any

from pydantic import BaseModel


class ConfigModel(BaseModel):
    """Base model for all config classes with dict-like access for Snakemake compatibility."""

    def __getitem__(self, key: str) -> Any:
        """Enable: config['key']."""
        return getattr(self, key)

    def __contains__(self, key: str) -> bool:
        """Enable: 'key' in config."""
        return hasattr(self, key)

    def get(self, key: str, default: Any = None) -> Any:
        """Enable: config.get('key', default)."""
        return getattr(self, key, default)

    def keys(self) -> Iterator[str]:
        """Enable: config.keys()."""
        return iter(self.model_fields.keys())

    def values(self) -> Iterator[Any]:
        """Enable: config.values()."""
        return (getattr(self, k) for k in self.model_fields.keys())

    def items(self) -> Iterator[tuple[str, Any]]:
        """Enable: config.items()."""
        return ((k, getattr(self, k)) for k in self.model_fields.keys())
