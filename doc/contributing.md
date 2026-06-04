<!-- SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur> -->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# Contributing

We welcome anyone interested in contributing to this project, be it with new
ideas, suggestions, by filing bug reports or contributing code to our [GitHub
repository](https://github.com/PyPSA/PyPSA-Eur).

## Where to start

* If you already have some code changes, you can submit them directly as a [pull request](https://github.com/PyPSA/pypsa-eur/pulls).
* If you are wondering where we would greatly appreciate your efforts, check out the `help wanted` tag in the [issues list](https://github.com/PyPSA/pypsa-eur/issues) and initiate a discussion there.
* If you start working on a feature in the code, let us know by opening an issue or a draft pull request.
  This helps all of us to keep an overview on what is being done and helps to avoid a situation where we
  are doing the same work twice in parallel.

## Setting up the development environment

For linting, formatting and checking your code contributions
against our guidelines (e.g. we use [Black](https://github.com/psf/black) as code style
use [pre-commit](https://pre-commit.com/index.html):

1. Install [pixi](https://pixi.sh/latest/).
1. Usage:
    * To automatically activate `pre-commit` on every `git commit`: Run `pre-commit install`
    * To manually run it: `pre-commit run --all`

!!! note
    Note that installing `pre-commit` locally is not strictly necessary. If you create a Pull Request the `pre-commit CI` will be triggered automatically and take care of the checks.

For all code contributions we follow the four eyes principle (two person principle), i.e. all suggested code
including our own are reviewed by a second person before they are incorporated into our repository.

If you are unfamiliar with pull requests, the GitHub help pages have a nice [guide](https://help.github.com/en/articles/about-pull-requests).

To **discuss** with other PyPSA users, organise projects, share news, and get in touch with the community you can use the [Discord server](https://discord.gg/AnuJBk23FU).

## Contributing to the documentation

We strive to keep documentation useful and up to date for all PyPSA users. If you encounter an area where documentation is not available or insufficient, we very much welcome your contribution. Here is How To:

1. Install [pixi](https://pixi.sh/latest/).
2. Make your changes in the corresponding `.md` file under `pypsa-eur/doc`.
3. Compile your changes by running `pixi run build-docs site` from the project root.
   You can also preview live with `mkdocs serve`. HTML files to review can be found under `site/`.
4. Contribute your documentation in a pull request ([here is a guide](https://help.github.com/en/articles/about-pull-requests)).
