HPC Configuration
=================

This page describes how to run **VirMake** on a *Slurm* cluster, and which
parts of the configuration can be adjusted.

Profiles
--------

:program:`virmake prep` generates a Snakemake profile for each run
directory:

.. code-block:: text

    myrun/config/
    |-- config.yaml           # workflow profile (always used)
    `-- slurm/
        `-- config.yaml       # Slurm profile (used with --slurm)

Both files are plain YAML and can be edited freely. :program:`virmake prep`
with ``--force`` regenerates them, overwriting edits.

The workflow profile
--------------------

``config/config.yaml`` configures the workflow regardless of where it
runs. As generated it looks like this (``conda-prefix`` holds the
environment directory chosen during :program:`virmake prep`):

.. code-block:: yaml

    software-deployment-method: conda
    conda-prefix: /cluster/shared/virmake/envs
    show-failed-logs: true
    rerun-incomplete: true
    keep-going: true
    printshellcmds: true
    jobname: '{rule}.{jobid}'
    max-jobs-per-second: 1

.. list-table::
   :header-rows: 1
   :widths: 32 68

   * - Setting
     - Description
   * - ``software-deployment-method``
     - Every rule runs inside its own conda environment.
   * - ``conda-prefix``
     - Directory the environments are created in (the ``--envs-dir``
       chosen at preparation).
   * - ``show-failed-logs``
     - Logs of failed jobs are printed at the end of a run.
   * - ``rerun-incomplete``
     - Jobs interrupted mid-run are re-run automatically on the next
       invocation.
   * - ``keep-going``
     - When a job fails, independent jobs continue.
   * - ``printshellcmds``
     - Every job prints the command it runs.
   * - ``jobname``
     - Cluster job names show rule and job id (e.g. ``fastp.12`` in
       ``squeue``).
   * - ``max-jobs-per-second``
     - Job submissions are throttled to one per second, which is easy
       on the scheduler.

Slurm accounts
--------------

Many clusters require an account on every job submission. Accounts are set
with ``virmake run --slurm --account`` (see :ref:`Running on a cluster`);
the value is stored in the Slurm profile ``config/slurm/config.yaml`` under
``default-resources.slurm_account`` and can also be edited there by hand.

The Slurm profile
-----------------

``config/slurm/config.yaml`` is used whenever the workflow runs with
``--slurm``. As generated it looks like this:

.. code-block:: yaml

    executor: slurm
    default-resources:
      slurm_account: ''
      runtime: 120 * attempt
      mem_mb: 7500 * attempt
    max-status-checks-per-second: 0.1
    latency-wait: 60
    restart-times: 2
    jobs: 200

.. list-table::
   :header-rows: 1
   :widths: 38 62

   * - Setting
     - Description
   * - ``default-resources.slurm_account``
     - The Slurm account (see above). Once ``--account`` has been given,
       the account appears here.
   * - ``default-resources.runtime``
     - Runtime per job in minutes, scaled by retry attempt: a job first
       tried with 120 minutes gets 240 on its first retry, 360 on its
       second.
   * - ``default-resources.mem_mb``
     - Memory per job in MB, scaled the same way.
   * - ``max-status-checks-per-second``
     - Job status is polled once every 10 seconds.
   * - ``latency-wait``
     - After a job finishes, up to 60 seconds are allowed for its output
       files to appear on the shared filesystem.
   * - ``restart-times``
     - Failed jobs are retried up to twice, each time with the larger memory and
       runtime that follow from the scaling above.
   * - ``jobs``
     - At most 200 jobs run at the same time.

Jobs are submitted to the ``normal`` partition, the *Slurm* executor's
default. To use a different partition, add it under ``default-resources``:

.. code-block:: yaml

    default-resources:
      slurm_partition: long

Individual workflow steps request their own resources where needed (the
assembly step, for example, asks for substantially more memory); the
values above apply to jobs without explicit requests, and on retries. The
per-step requests are configured in the workflow configuration file (see
:doc:`config-params`).

Storage of environments and databases
-------------------------------------

The database and software environment directories are chosen when the run
is prepared and are shared between runs; they can become large (see the
``--db-dir`` and ``--envs-dir`` options in :doc:`usage`). On clusters, it may
be preferable to store them at shared project storage locations.

Network access
--------------

Software environments are created and databases downloaded on the
machines that run the jobs. On clusters where compute nodes lack
internet access, prefetch both from the login node first (see
:ref:`Environments and databases` in :doc:`usage`).

Resuming and re-running
-----------------------

**VirMake** records finished steps as it goes: an interrupted run is
simply started again with the same command, and only the remaining jobs
run (``rerun-incomplete`` and ``restart-times`` are already set in the
profiles). If a run crashed hard while holding the lock, pass Snakemake's
unlock option:

.. code-block:: console

    virmake run myrun all --slurm --unlock

Any other Snakemake option can be passed through the same way (see
:doc:`usage`).
