..
  SPDX-FileCopyrightText: 2019-2023 The PyPSA-Eur Authors

  SPDX-License-Identifier: CC-BY-4.0

#######################
Contributing
#######################

We welcome anyone interested in contributing to this project, be it with new
ideas, suggestions, by filing bug reports or contributing code to our `GitHub
repository <https://github.com/PyPSA/PyPSA-Eur>`_.

* If you already have some code changes, you can submit them directly as a `pull request <https://github.com/PyPSA/pypsa-eur/pulls>`_.
* If you are wondering where we would greatly appreciate your efforts, check out the ``help wanted`` tag in the `issues list <https://github.com/PyPSA/pypsa-eur/issues>`_ and initiate a discussion there.
* If you start working on a feature in the code, let us know by opening an issue or a draft pull request.
  This helps all of us to keep an overview on what is being done and helps to avoid a situation where we
  are doing the same work twice in parallel.

For linting, formatting and checking your code contributions
against our guidelines (e.g. we use `Black <https://github.com/psf/black>`_ as code style
use `pre-commit <https://pre-commit.com/index.html>`_:

1. Installation ``conda install -c conda-forge pre-commit`` or ``pip install pre-commit``
2. Usage:
    * To automatically activate ``pre-commit`` on every ``git commit``: Run ``pre-commit install``
    * To manually run it: ``pre-commit run --all``

Note that installing `pre-commit` locally is not strictly necessary. If you create a Pull Request the `pre-commit CI` will be triggered automatically and take care of the checks.

For all code contributions we follow the four eyes principle (two person principle), i.e. all suggested code
including our own are reviewed by a second person before they are incorporated into our repository.

If you are unfamiliar with pull requests, the GitHub help pages have a nice `guide <https://help.github.com/en/articles/about-pull-requests>`_.

To ask and answer general usage questions, join the `PyPSA and PyPSA-Eur mailing list <https://groups.google.com/forum/#!forum/pypsa>`_.
