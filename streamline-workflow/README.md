<!--
SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
SPDX-License-Identifier: CC-BY-4.0
-->

# Streamline Workflow Development

This folder contains development work for the PyPSA-EUR workflow streamlining project (GitHub Discussion #1529).

## Files

- `STREAMLINE_WORKFLOW_IMPLEMENTATION_PLAN.md` - Detailed implementation plan
- `test_horizon_wildcard_validated.smk` - **Working test Snakefile (validated)**
- `test_streamlined.yaml` - Test configuration file

## Testing

To test the streamlined workflow:

```bash
snakemake -s rules/test_streamlined.smk test_overnight -n --configfile config/test/config.overnight.yaml
```

## Key Findings

✅ **Validated Concepts**:
- `{horizon}` wildcard works correctly
- 4-step workflow: `base → clustered → composed_{horizon} → solved_{horizon}`
- Multi-run support with `{run}` wildcard
- Config-driven approach replacing wildcard expansion

## Status

- Basic workflow structure: ✅ Tested and working
- Multi-horizon processing: ✅ Validated
- Sequential dependencies: ⏳ Next to implement
- Full script integration: ⏳ Next to implement

This folder can be deleted once the streamlined workflow is fully implemented and integrated into the main codebase.