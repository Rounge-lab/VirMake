=============
Contributing
=============

This page is for developers and contributors working on **VirMake** itself.
To use the pipeline, see :doc:`install` and :doc:`usage`.

Development installation
------------------------

Development happens on the ``dev`` branch. The development environment is
managed with `pixi <https://pixi.sh>`_, which installs **VirMake** as an
editable package together with the development tools:

.. code-block:: console

    git clone https://github.com/Rounge-lab/VirMake.git
    cd VirMake
    git switch dev
    pixi install

Pixi tasks:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Task
     - Description
   * - ``pixi run lint`` / ``pixi run lint-fix``
     - Ruff lint and format checks, with and without automatic fixes.
   * - ``pixi run pytest``
     - Unit tests with a coverage report.
   * - ``pixi run test``
     - Lint plus unit tests — the same suite CI runs.
   * - ``pixi run smoke``
     - Prepares a run directory from the bundled test data and performs a
       dryrun of the whole workflow.
   * - ``pixi run typecheck``
     - Static type checking with Pyright.
   * - ``pixi run build``
     - Builds the sdist and wheel into ``dist/``.

Install the pre-commit hooks once; they run the lint checks on every
commit and the test suite on every push to ensure compliance with CI:

.. code-block:: console

    pixi run pre-commit-install

The documentation builds in its own environment:

.. code-block:: console

    pixi run -e docs docs

Testing
-------

The test suite has three layers:

* **Unit tests** (``tests/``) cover the command-line interface, run
  directory preparation and configuration generation. They run on any
  platform, including Windows and macOS development machines.
* **Integration tests** execute parts of the actual Snakemake workflow.
  They are Linux-only and skip on other platforms.
* The **smoke** task — preparation of a test-data run followed by a dryrun
  of the full workflow — catches broken rule dependencies, missing input
  files and configuration mistakes that unit tests cannot see.

Continuous integration
----------------------

GitHub Actions runs on every push to ``dev`` and on every pull request:

* ``checks (ubuntu-latest)`` — lint, tests and the pre-commit hygiene
  checks
* ``checks (windows-latest)`` — the same test suite on Windows
* ``smoke (ubuntu-latest)`` — the smoke dryrun

Pull requests into ``main`` must pass all three legs before merging.

Branches and pull requests
--------------------------

* Day-to-day development happens on ``dev``, via direct commits or
  short-lived feature branches with pull requests into ``dev``.
* ``main`` only receives changes through pull requests, protected by the
  checks above.
* Releases are cut from ``main`` only.

Releasing
---------

Releases are created by the *Release* workflow. Tags and GitHub releases
are never made by hand.

1. Make sure ``dev`` is green in CI.

2. Bump ``__version__`` in ``src/virmake/__init__.py``.

3. Merge ``dev`` into ``main`` via pull request and wait for the required
   checks.

4. On GitHub, open *Actions → Release → Run workflow*, select the ``main``
   branch and enter the version number. The workflow refuses to run from
   any other branch.

5. The workflow then:

   * runs the test suite and the smoke dryrun,
   * verifies that the requested version matches the package and that
     the tag does not exist yet,
   * builds the sdist and wheel,
   * creates the GitHub release with both artifacts attached, the
     sdist's sha256 in the release notes and an automatic changelog,
   * if a bioconda recipe exists, updates it in the ``Rounge-lab``
     bioconda-recipes fork with the new version and sha256, pushes a
     branch, and prints the link for opening the pull request into
     bioconda.

6. Open the bioconda pull request from the printed link and submit it
   for review.
