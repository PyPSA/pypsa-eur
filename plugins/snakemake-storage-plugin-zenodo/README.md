<!--
SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
SPDX-License-Identifier: CC-BY-4.0
-->

# Snakemake Storage Plugin: Zenodo

A Snakemake storage plugin for downloading files from Zenodo with local caching, checksum verification, and adaptive rate limiting.

## Features

- **Local caching**: Downloads are cached to avoid redundant transfers (can be disabled)
- **Checksum verification**: Automatically verifies MD5 checksums from Zenodo API
- **Rate limit handling**: Automatically respects Zenodo's rate limits using `X-RateLimit-*` headers with exponential backoff retry
- **Concurrent download control**: Limits simultaneous downloads to prevent overwhelming Zenodo
- **Progress bars**: Shows download progress with tqdm
- **Immutable URLs**: Returns mtime=0 since Zenodo URLs are persistent
- **Environment variable support**: Configure via environment variables for CI/CD workflows

## Installation

From the pypsa-eur repository root:

```bash
pip install -e plugins/snakemake-storage-plugin-zenodo
```

## Configuration

The Zenodo storage plugin works alongside other storage providers (like HTTP). Snakemake automatically routes URLs to the correct provider based on the URL pattern.

Register additional settings in your Snakefile if you want to customize the defaults:

```python
# Optional: Configure Zenodo storage with custom settings
# This extends the existing storage configuration (e.g., for HTTP)
storage zenodo:
    provider="zenodo",
    cache="~/.cache/snakemake-pypsa-eur",  # Default location
    max_concurrent_downloads=3,  # Download max 3 files at once
```

If you don't explicitly configure it, the plugin will use default settings automatically.

### Settings

- **cache** (optional): Cache directory for downloaded files
  - Default: Platform-dependent user cache directory (via `platformdirs.user_cache_dir("snakemake-pypsa-eur")`)
  - Set to `""` (empty string) to disable caching
  - Files are cached here to avoid re-downloading
  - Environment variable: `SNAKEMAKE_STORAGE_ZENODO_CACHE`

- **skip_remote_checks** (optional): Skip metadata checking with Zenodo API
  - Default: `False` (perform checks)
  - Set to `True` or `"1"` to skip remote existence/size checks (useful for CI/CD)
  - Environment variable: `SNAKEMAKE_STORAGE_ZENODO_SKIP_REMOTE_CHECKS`

- **max_concurrent_downloads** (optional): Maximum concurrent downloads
  - Default: `3`
  - Controls how many Zenodo files can be downloaded simultaneously
  - No environment variable support

## Usage

Use Zenodo URLs directly in your rules. Snakemake automatically detects zenodo.org URLs and routes them to this plugin:

```python
rule download_data:
    input:
        storage("https://zenodo.org/records/3520874/files/natura.tiff"),
    output:
        "resources/natura.tiff"
    shell:
        "cp {input} {output}"
```

Or if you configured a tagged storage entity:

```python
rule download_data:
    input:
        storage.zenodo(
            "https://zenodo.org/records/3520874/files/natura.tiff"
        ),
    output:
        "resources/natura.tiff"
    shell:
        "cp {input} {output}"
```

The plugin will:
1. Check if the file exists in the cache (if caching is enabled)
2. If cached, copy from cache (fast)
3. If not cached, download from Zenodo with:
   - Progress bar showing download status
   - Automatic rate limit handling with exponential backoff retry
   - Concurrent download limiting
   - MD5 checksum verification against Zenodo API metadata
4. Store in cache for future use (if caching is enabled)

### Example: CI/CD Configuration

For continuous integration environments where you want to skip caching and remote checks:

```yaml
# GitHub Actions example
- name: Run snakemake workflows
  env:
    SNAKEMAKE_STORAGE_ZENODO_CACHE: ""
    SNAKEMAKE_STORAGE_ZENODO_SKIP_REMOTE_CHECKS: "1"
  run: |
    snakemake --cores all
```

## Rate Limiting and Retry

Zenodo API limits:
- **Guest users**: 60 requests/minute
- **Authenticated users**: 100 requests/minute

The plugin automatically:
- Monitors `X-RateLimit-Remaining` header
- Waits when rate limit is reached
- Uses `X-RateLimit-Reset` to calculate wait time
- Retries failed requests with exponential backoff (up to 5 attempts)
- Handles transient errors: HTTP errors, timeouts, checksum mismatches, and network issues

## URL Handling

- Only handles URLs containing `zenodo.org`
- Other HTTP(S) URLs are handled by the standard `snakemake-storage-plugin-http`
- Both plugins can coexist in the same workflow

### Plugin Priority

When using `storage()` without specifying a plugin name, Snakemake checks all installed plugins:
- **Zenodo plugin**: Only accepts zenodo.org URLs (`is_valid_query` returns True only for zenodo.org)
- **HTTP plugin**: Accepts all HTTP/HTTPS URLs (including zenodo.org)

If both plugins are installed, zenodo.org URLs are ambiguous - both plugins accept them.
Typically snakemake would raise an error: **"Multiple suitable storage providers found"** if you try to use `storage()` without specifying which plugin to use, ie. one needs to explicitly call the Zenodo provider for zenodo.org URLs using `storage.zenodo(url)` instead of `storage(url)`,
but we monkey-patch the http plugin to refuse zenodo.org urls.

## License

MIT License
