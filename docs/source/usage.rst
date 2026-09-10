=====
Usage
=====

Basic Usage
-----------

.. note::

    When installing from source, make sure to run all commands from the root
    directory of the repository.

1. Activate the **VirMake** conda environment:

.. code-block:: console

    conda activate ./venv

2. Run the **VirMake** pipeline:

Command-line interface
^^^^^^^^^^^^^^^^^^^^^^

Synopsis
""""""""

.. code-block:: console

    virmake -h

usage: virmake [OPTIONS] COMMAND [ARGS]...

Options
"""""""

.. option:: -h, --help

    Show the help message and exit.

Commands
""""""""

.. option:: clean

    Clean VirMake directory.

.. option:: db

    Downloads databases for steps defined in config. Databases will be
    downloaded as needed if not downloaded using this option.

.. option:: prep

    Downloads and creates all environments needed to run the workflow offline.

.. option:: run

    Runs the main workflow.

Environment preparation
^^^^^^^^^^^^^^^^^^^^^^^

Running VirMake requires the preparation of software and databases, which may
take a substantial amount of time. To prepare a run by downloading and setting
up environments, use the :program:`virmake prep` command:

.. program:: virmake prep

Synopsis
""""""""

.. code-block:: console

    virmake prep [OPTIONS]

Options
"""""""

.. option:: -h, --help

    Show the help message and exit.

.. option:: -c, --threads
   **Type:** `int`
   **Default:** 24

    Number of threads to use per multi-threaded job.

Database download
^^^^^^^^^^^^^^^^^

To prepare a run by downloading and setting up databases, use the :program:`virmake db`
command:

.. program:: virmake db

Synopsis
""""""""

.. code-block:: console

    virmake db [OPTIONS]

Options
"""""""

.. option:: -h, --help

    Show the help message and exit.

.. option:: -n, --dryrun

    Test execution of the command.

.. option:: -c, --threads
   **Type:** `int`
   **Default:** 24

    Maximum number of threads to use per multi-threaded job.

These commands will set up all requirements for running the **VirMake**
pipeline, including the prerequisites for any steps specified in the workflow
configuration file ``config/params.yaml`` under :confval:`rule_inclusion` (See
:ref:`config-workflow`).

Running the workflow
^^^^^^^^^^^^^^^^^^^^

.. program:: virmake run

Synopsis
""""""""

.. code-block:: console

    virmake run [OPTIONS] <WORKFLOW>

Positional arguments
""""""""""""""""""""

.. describe:: <WORKFLOW>

    - all (default)
    - qc
    - assembly
    - identification
    - mapping
    - taxonomy
    - function
    - stats

    Specify which part of the workflow to execute.
    For more details on the workflow, see :ref:`workflow`.

Options
"""""""

.. option:: -h, --help

    Show the help message and exit.

.. option:: -n, --dryrun

    Test execution of the command.

.. option:: -s, --slurm

    Use `Slurm <https://slurm.schedmd.com/overview.html>`_ cluster to run
    parallel jobs.

.. option:: -c, --threads
    **Type:** `int`
    **Default:** 24

    Maximum number of threads to use per multi-threaded job.

..  NOTE! Is this supposed to be config/config.yaml ?
    NOTE! Why is this parsed as a string while other file options are parsed as
    paths?
.. option:: -p, --profile
    **Type:** `str`
    **Default:** ``./config``

    Snakemake profile *e.g.* for cluster execution.

.. option:: -d, --workflow-dir
    **Type:** `path`
    **Default:** ``./workflow``

    Location to run virmake.

..  NOTE! Should the help message mention that this is the Workflow parameters
   config file?
.. option:: -C, --config-file
    **Type:** `path`
    **Default:** ``./config/params.yaml``

    Config file generated during virmake setup (See :ref:`config-workflow`).

.. option:: -T, --jobs_at_once
    **Type:** `int`
    **Default:** 3

    Number of jobs to add to queue at once.

Please note that the :option:`-c` (or :option:`--threads`) option controls both
the number of jobs and the number of threads in each job when running without
using *Slurm*. Multiple jobs are used for assembly and other steps that can be
run separately for each sample. For example, running ``./virmake run -c 8``
will start 8 jobs using 8 threads each (for a total of 64 threads) during
assembly.

