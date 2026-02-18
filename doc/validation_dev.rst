.. SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
..
.. SPDX-License-Identifier: CC-BY-4.0

##########################################
Validation
##########################################

PyPSA-Eur uses `Pydantic <https://docs.pydantic.dev/>`_ models for validation.
This system provides type checking, default values, and documentation in a single place.

Configuration
=============

The configuration validation system consists of:

- **Pydantic models** in ``scripts/lib/validation/config/`` that define all options and validates the snakemake config.
- **Auto-generated files**: ``config/config.default.yaml`` and ``config/schema.default.json``.

Adding a New Config Option
--------------------------

To add a new option to an existing config section, edit the corresponding module in
``scripts/lib/validation/config/``. Each field uses Pydantic's ``Field()`` function.

For example, the ``logging`` section in ``scripts/lib/validation/config/__init__.py``:

.. code-block:: python

   from typing import Literal
   from pydantic import Field
   from scripts.lib.validation.config._base import ConfigModel

   class LoggingConfig(ConfigModel):
       """Configuration for top level `logging` settings."""

       # ... existing fields ...

       # An option with default 0.5, float type and between 0 and 1. If anything else is passed,
       # the validation will fail
       new_option: float = Field(
           0.5,  # default value
           description="Threshold for the new feature.",  # shown in docs and IDE
           ge=0,  # greater than or equal
           le=1,  # less than or equal
           examples=[0.3, 0.7],  # example values for docs
       )

For more Field parameters, see the `Pydantic Field documentation <https://docs.pydantic.dev/latest/concepts/fields/>`__.

.. note::
    If you are making config changes in a fork of PyPSA-Eur to meet your project-specific needs,
    you should instead use the config updater class described in :ref:`soft_fork_ext`.

Adding a New Config Section
---------------------------

To add a nested config section, define a helper class and add it to an existing config.
For example, adding a ``file`` section to ``LoggingConfig``:

.. code-block:: python

   class _LoggingFileConfig(ConfigModel):
       """Configuration for logging to file."""

       enabled: bool = Field(False, description="Enable file logging.")
       path: str = Field("logs/pypsa.log", description="Log file path.")
       format: str | None = Field(None, description="Custom log format for that file.")


   class LoggingConfig(ConfigModel):
       """Configuration for top level `logging` settings."""

       # ... existing fields ...

       file: _LoggingFileConfig = Field(
           default_factory=_LoggingFileConfig,
           description="File logging configuration.",
       )


There is one python module for each top level configuration. Helper classes for nested
keys usee underscore prefix (e.g., ``_LoggingFileConfig``) by convention.

.. note::
    If you are making config changes in a fork of PyPSA-Eur to meet your project-specific needs,
    you should instead use the config updater class described in :ref:`soft_fork_ext`.

Regenerating Config Files
-------------------------

Snakemake will only read from the ``config/config.default.yaml``, which needs to be generated
after making changes to the Pydantic model. To regenerate the default config and JSON
schema:

.. code-block:: console

   $ pixi run generate-config

This updates ``config/config.default.yaml`` and ``config/schema.default.json``.
For example, the two examples above would now generate:

.. code-block:: yaml

   logging:
     level: INFO
     format: "%(levelname)s:%(name)s:%(message)s"
     new_option: 0.5
     file:
       enabled: false
       path: logs/pypsa.log
       format: null

Always commit these regenerated files alongside your model changes.

Custom Validators
-----------------

For validation logic beyond simple type checks and constraints, Pydantic provides
``field_validator`` (for single fields) and ``model_validator`` (for cross-field validation).

**Field Validator**: Validate a single field's value. For example, ensuring the log level
is uppercase:

.. code-block:: python

   from pydantic import Field, field_validator
   from scripts.lib.validation.config._base import ConfigModel

   class LoggingConfig(ConfigModel):
       """Configuration for top level `logging` settings."""

       level: str = Field("INFO", description="Logging level.")

       @field_validator("level")
       @classmethod
       def validate_level(cls, v):
           if v.upper() != v:
               raise ValueError("Logging level must be uppercase (e.g., 'INFO', 'DEBUG').")
           return v

**Model Validator**: Validate relationships between multiple fields. For example,
ensuring the file path is set when file logging is enabled:

.. code-block:: python

   from pydantic import Field, model_validator
   from scripts.lib.validation.config._base import ConfigModel

   class LoggingConfig(ConfigModel):
       """Configuration for top level `logging` settings."""

       file_enabled: bool = Field(False, description="Enable file logging.")
       file_path: str | None = Field(None, description="Log file path.")

       @model_validator
       def check_file_path_required(self):
           if self.file_enabled and not self.file_path:
               raise ValueError("file_path is required when file_enabled is True.")
           return self

Again, find more information in the Pydantic documentation on
`Field Validators <https://docs.pydantic.dev/latest/concepts/validators/#field-validators>`_
and `Model Validators <https://docs.pydantic.dev/latest/concepts/validators/#model-validators>`_.

.. _soft_fork_ext:

Extending for Soft-Forks
------------------------

If you maintain a soft-fork of PyPSA-Eur with custom config options, you have two approaches:

**Allow extra fields**: The ``ConfigSchema`` uses ``extra="allow"`` by default, so
unrecognized config keys won't cause validation errors. Your custom options will pass
through without type checking. Only if you changed existing config settings, you will
need to adjust the schema. But you will lose the sync of Pydantic model and defaults
YAML, which is currently enforced via an upstream CI job.

**Extend the schema**: It is better to add full validation of your additional
configuration.
The cleanest way to do this is to use the config updater base class that we make available.
You impose config changes in subclasses of the base class and by importing those into `scripts.lib.validation.config_updates.py` they will be used to automatically overwrite the configuration.

In the below example, two updates are made to the default config.

.. code-block:: python

    from typing import Literal

    from pydantic import BaseModel, Field

    from scripts.lib.validation.config._base import ConfigUpdater
    from scripts.lib.validation.config._schema import ConfigSchema


    class ClusteringConfigUpdater(ConfigUpdater):
        name: str = "update_clustering"

        def update(self) -> type[ConfigSchema]:
            # To update and existing config item, we need it's most recent state, as defined in `self.base_config`
            clustering_config = self.base_config().clustering.__class__
            mode_config = clustering_config.model_fields["mode"]

            current_description = mode_config.description or ""
            new_description = current_description + " (extra) foobar: new item."
            new_list = Literal[mode_config.annotation, "foobar"]

            clustering_schema = self._apply_updates(
                __base__=clustering_config,
                mode=(new_list, Field(mode_config.default, description=new_description)),
            )
            new_schema = self._apply_updates(
                clustering=(clustering_schema, Field(default_factory=clustering_schema))
            )

            return new_schema


    class MyNewConfigSection(BaseModel):
        my_new_field: str = Field("foo")


    class NewConfigItem(ConfigUpdater):
        name: str = "new_section"

        def update(self) -> type[ConfigSchema]:
            new_schema = self._apply_updates(
                new_section=(MyNewConfigSection, Field(default_factory=MyNewConfigSection)),
            )
            return new_schema

If this code were stored in the script ``scripts/_my_config_updates.py`` then ``scripts.lib.validation.config_updates.py`` would now include:

.. code-block:: python

    import scripts._my_config_updates

This is sufficient for both updates to be imported.

.. admonition:: Config filename

    When generating the config files with the above example (``pixi run generate-config``),
    you would now generate ``config/config.default.update_clustering.new_section.yaml``.
    To override the base ``config/config.default.yaml``, you can set the ``name`` property of your updater classes to empty strings: ``""``.

.. admonition:: Chaining updates

    Several separate update scripts can exist and be used to create chained updates of the schema.
    They will be used to update the schema in the order they appear in ``scripts.lib.validation.config_updates.py``.
    this means that you can update the same config item multiple times.
    If you are importing config changes from a submodule and you want to catch cases where you are both updating the same config item, you can add a check in your ``update`` method, such as:

    .. code:: python

        if clustering_config != ClusteringConfig:
            raise ValueError(
                "You are trying to update the clustering config item after it has already been updated by another config updater."
                " This could have unexpected consequences."
            )

.. autoclass:: lib.validation.config._base::ConfigUpdater
    :members: name, update, _apply_updates