# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""Base classes for config validation models."""

from abc import ABC, abstractmethod
from collections.abc import Iterator
from typing import TYPE_CHECKING, Any, TypeVar, overload

from pydantic import BaseModel, create_model

if TYPE_CHECKING:
    from scripts.lib.validation.config import ConfigSchema

_registry: list[type["ConfigUpdater"]] = []
T = TypeVar("T", bound=type[BaseModel])


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


class ConfigUpdater(ABC):
    """Abstract base class for updating the PyPSA-Eur base config"""

    def __init__(self, base_config: type["ConfigSchema"]):
        self.base_config = base_config

    def __init_subclass__(cls):
        """Register any subclasses so they are automatically captured for use later."""
        super().__init_subclass__()
        _registry.append(cls)

    @property
    @abstractmethod
    def name(self) -> str:
        """
        Name of custom config.
        If not an empty string, this name will be appended to the default config filename when the config file is saved, i.e.:

        `config.default.<prior-update-names>.{self.name}.<post-update-names>.yaml`
        """

    @abstractmethod
    def update(self) -> type["ConfigSchema"]:
        """Function in which the custom config schema is created and returned."""

    @overload
    def _apply_updates(
        self,
        __base__: None = None,
        __doc__: str | None = None,
        **updates: tuple[Any, Any],
    ) -> type["ConfigSchema"]: ...

    @overload
    def _apply_updates(
        self, __base__: T, __doc__: str | None = None, **updates: tuple[Any, Any]
    ) -> T: ...

    def _apply_updates(
        self,
        __base__: T | None = None,
        __doc__: str | None = None,
        **updates: tuple[Any, Any],
    ) -> T:
        """
        Helper function to apply updates to a pydantic model.

        Updates should come in the form of a type (can be another pydantic model) and a pydantic.Field with at least a default value.

        Parameters
        ----------
        __base__ : type[pydantic.BaseModel] | None, optional
            Optional base model to use instead of the class-level base config.
        __doc__ : str, optional
            Optional docstring for the new config model.
            If not provided, the docstring of the base model will be used.
        **updates : tuple[Any, Any]
            The actual updates to apply, in the form of field_name=(type, Field(...)).

        Returns
        -------
        Updated config model with the new fields included.
        """
        base: T | type[ConfigSchema]
        if __base__ is None:
            base = self.base_config
        else:
            base = __base__
        return create_model(
            base.__name__,
            __doc__=__doc__ or base.__doc__,
            __base__=base,
            **updates,
        )
