=====
Usage
=====

Basic Usage
-----------

.. note::

    When installing from source, make sure to run all commands from the root
    directory of the repository.

#. Activate the **VirMake** conda environment:

.. code-block:: console

    conda activate ./venv

#. Run the **VirMake** the pipeline:

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
up environments, use the :option:`prep` command:


Synopsis
""""""""

.. code-block:: console

    virmake prep -h

Usage: virmake prep [OPTIONS]

Options
"""""""

.. option:: -h, --help

    Show the help message and exit.

.. option:: -c, --threads INTEGER

    Number of threads to use per multi-threaded job

To prepare a run by downloading and setting up databases, use the db command:

Synopsis
""""""""

.. code-block:: console

    virmake db -h

Usage: virmake db [OPTIONS]

Options
"""""""

.. option:: -h, --help

    Show the help message and exit.

.. option:: -c, --threads INTEGER

    Number of threads to use per multi-threaded job

.. option:: -n, --dryrun

    Test execution of the command

These commands will set up all requirements for running the **VirMake**
pipeline, including the prerequisites for any steps specified in the workflow
configuration file ``config/params.yaml`` under ``rule_inclusion`` (See
:ref:`config-workflow`).

Running the workflow
^^^^^^^^^^^^^^^^^^^^

Synopsis
""""""""

.. code-block:: console

    virmake run -h

Usage: virmake run [OPTIONS]

Options
"""""""

.. option:: -h, --help

    Show the help message and exit.

.. option:: -c, --threads INTEGER

    Number of threads to use per multi-threaded job.

.. option:: -n, --dryrun

    Test execution of the command.

.. option:: -p, --profile TEXT

    Snakemake profile e.g. for cluster execution.

.. option:: -d, --workflow-dir PATH

    Location to run virmake.

.. option:: -C, --config-file PATH

    Config file generated during virmake setup.

.. option:: -s, --slurm

    Use `*Slurm* <https://slurm.schedmd.com/overview.html>`_ cluster to run
    parallel jobs.

.. option:: -T, --jobs_at_once INTEGER

    Number of jobs to add to queue at once.

Please note that the :option:`-c` (or :option:`--threads`) option controls both
the number of jobs and the number of threads in each job when running without
using *Slurm*. Multiple jobs are used for assembly and other steps that can be
run separately for each sample. For example, running ``./virmake run -c 8``
will start 8 jobs using 8 threads each (for a total of 64 threads) during
assembly.

