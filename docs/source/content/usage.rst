=====
Usage
=====

Basic Usage
-----------

The examples below assume the ``virmake`` command is on your path (see
:doc:`install`). If you installed **VirMake** into a conda environment,
activate it first:

.. code-block:: console

    conda activate virmake

**VirMake** organises every analysis in its own *run directory*, created by
:program:`virmake prep`. The run directory records the input samples, the
workflow configuration and the locations of the shared databases and
software environments. All further commands operate on this directory.

Trying **VirMake** on the bundled test data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. Prepare a run directory with a small synthetic dataset:

.. code-block:: console

    virmake prep testrun --test-data

2. Review the generated configuration in ``testrun/config/``.

3. Show what the workflow would do, without running anything:

.. code-block:: console

    virmake run testrun all --dryrun

4. Run the whole workflow on it:

.. code-block:: console

    virmake run testrun all --cores 8

.. note::

   The synthetic reads are tiny, but a full ``all`` run still downloads all
   databases and creates the software environments the first time they are
   needed. Their location can be specified using the ``--db-dir`` and
   ``--envs-dir`` options with the ``virmake prep`` command.


Running your own data
^^^^^^^^^^^^^^^^^^^^^

Input samples are recognised by file name, and multiple samples may be
present in the same folder:

* Raw reads must be paired-end and gzip compressed, named
  ``{sample}_1.fastq.gz`` and ``{sample}_2.fastq.gz``.
* Reads that have already passed quality control must be paired-end and
  named ``{sample}_1.fastq.gz`` and ``{sample}_2.fastq.gz`` (or
  uncompressed ``.fastq``).
* Pre-assembled contigs are named ``{sample}.fasta``.

.. code-block:: console

    virmake prep myrun --reads myreads

The input options can be combined: for example, reads given with ``--reads`` are
used for quality control and mapping, while contigs given with ``--contigs`` are
used instead of assembling.

Command-line interface
----------------------

.. program:: virmake

Synopsis
^^^^^^^^

.. code-block:: console

    virmake [OPTIONS] COMMAND [ARGS]...

Options
^^^^^^^

.. option:: -h, --help

    Show the help message and exit.

Commands
^^^^^^^^

.. option:: prep

    Prepare a run directory: samples table, configuration and results directory.

.. option:: run

    Run the workflow on a prepared run directory.

.. option:: envs

    Prefetch the software environments needed by the configured workflow steps.

.. option:: db

    Prefetch the databases needed by the configured workflow steps.

Preparing a run
---------------

.. program:: virmake prep

Synopsis
^^^^^^^^

.. code-block:: console

    virmake prep [OPTIONS] RUNDIR

Positional arguments
^^^^^^^^^^^^^^^^^^^^

.. describe:: RUNDIR

    Path of the run directory to prepare. It is created if missing.
    **VirMake** refuses to overwrite a non-empty run directory unless
    ``--force`` is given.

Options
^^^^^^^

.. option:: -h, --help

    Show the help message and exit.

.. option:: --reads

    :Type: ``path``

    Directory with paired-end raw reads (gzip compressed, see `Basic Usage`_).

.. option:: --qc-reads

    :Type: ``path``

    Directory with reads that have already passed quality control. The trimming
    step is skipped for these samples.

.. option:: --contigs

    :Type: ``path``

    Directory with pre-assembled contigs. The assembly step is skipped for these
    samples.

.. option:: --test-data

    Download a small `synthetic dataset
    <https://doi.org/10.5281/zenodo.23388915>`_ into the run directory and
    use it as quality-controlled reads. This option cannot be combined
    with the input options above.

.. option:: --db-dir

    :Type: ``path``
    :Default: ``~/.local/share/virmake/databases``

    Directory for the workflow databases. Shared between runs.

.. option:: --envs-dir

    :Type: ``path``
    :Default: ``~/.local/share/virmake/envs``

    Directory for the workflow's software environments. Shared between runs.

.. option:: --yes

    Use the default database and environment locations without asking.

.. option:: --force

    Regenerate the samples table and the configuration in an existing run
    directory. ``results/`` is never touched.

Unless ``--yes`` is given, **prep** shows the database and environment locations
it intends to use and asks for confirmation. Both locations are created once and
shared between runs. On clusters, the default locations in the ``$HOME``
directory may not be preferred. You can specify their locations as follows:

.. code-block:: console

    virmake prep myrun --reads myreads \
        --db-dir /cluster/shared/virmake/databases \
        --envs-dir /cluster/shared/virmake/envs

Running the workflow
--------------------

.. program:: virmake run

Synopsis
^^^^^^^^

.. code-block:: console

    virmake run [OPTIONS] RUNDIR [WORKFLOW]

Positional arguments
^^^^^^^^^^^^^^^^^^^^

.. describe:: RUNDIR

    A run directory prepared by :program:`virmake prep`.

.. describe:: WORKFLOW

    :Type: ``str``
    :Default: ``all``
    :Alternatives:

      * ``qc``
      * ``assembly``
      * ``identification``
      * ``mapping``
      * ``taxonomy``
      * ``function``
      * ``stats``

    Which part of the workflow to execute, together with everything it depends
    on. See :doc:`workflow` for what each part does.

Options
^^^^^^^

.. option:: -h, --help

    Show the help message and exit.

.. option:: -n, --dryrun

    Test execution: show the jobs that would run, without executing anything.
    Use this to check a prepared run before committing to the real thing.

.. option:: --cores

    :Type: ``int``
    :Default: ``1``

    Number of cores **VirMake** may use at once. When running locally, as many
    jobs as fit within this budget are started in parallel.

.. option:: --slurm

    Submit the workflow's jobs to a SLURM cluster (see :ref:`Running on a
    cluster` and :doc:`config-hpc`).

.. option:: --account

    SLURM account to submit jobs under. The account is recorded in the run's
    SLURM configuration, so it needs to be given only once per run directory.

Further options
^^^^^^^^^^^^^^^

Any further options are passed through to the underlying Snakemake invocation;
for a full overview, refer to the official `Snakemake documentation
<https://snakemake.readthedocs.io/en/stable/executing/cli.html>`_.
Common recovery scenarios are shown in :doc:`config-hpc`.

Checks at the start of a run
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Before starting, **VirMake** checks that the databases and software environments
needed by the configured steps are in place, and warns about anything missing.
Missing pieces are then created automatically when the workflow reaches the step
that needs them. :program:`virmake envs` and :program:`virmake db` allow you to
fetch required environments and databases ahead of time (see :ref:`Environments
and databases`).

Where the results end up
^^^^^^^^^^^^^^^^^^^^^^^^

Results are written to ``results/output/``, with one subdirectory per step.
The main output tables most users are interested in are:

``results/output/mapping/rel_abundance_table.tsv``
   relative abundance of each vOTU per sample

``results/output/statistics/sample_stats.tsv``
   per-sample summary

``results/output/statistics/vOTU_stats.tsv``
   per-vOTU summary

For the full inventory of output files, see :doc:`workflow`.

.. note::

   Temporary files are cleaned up after a successful full run.

.. _Environments and databases:

Environments and databases
--------------------------

The tools used by the workflow (fastp, metaSPAdes, VirSorter2, ...) are
installed into per-tool conda environments, and the viral databases (CheckV,
DRAM-v, vConTACT3, ...) are downloaded once into the database directory chosen
during preparation. Nothing needs to be done manually: environments are created
on first use, and missing databases are downloaded by the workflow itself.

To prefetch everything ahead of time:

.. code-block:: console

    virmake envs myrun --threads 4
    virmake db myrun --threads 4

Both commands accept ``--threads`` to limit how many environments or downloads
are prepared at once; :program:`virmake db` additionally accepts ``-n`` to only
show what would be downloaded.

.. _Running on a cluster:

Running on a cluster
--------------------

On a SLURM cluster, submit jobs instead of running them locally:

.. code-block:: console

    virmake run myrun all --slurm --account myproject

After the first run, the account is recorded in the run directory and later
invocations only need ``--slurm``. Job parameters (partition, memory, runtime,
number of concurrent jobs) are set by the SLURM profile that :program:`virmake
prep` generated; see :doc:`config-hpc` for how to adjust it.

The run directory
-----------------

A prepared run directory looks like this:

.. code-block:: text

    myrun/
    |-- samples.tsv
    |-- test_data/            # only with --test-data
    |-- config/
    |   |-- params.yaml       # workflow configuration
    |   |-- config.yaml       # Snakemake profile
    |   `-- slurm/
    |       `-- config.yaml   # Snakemake SLURM profile
    `-- results/
        |-- output/
        |-- log/
        |-- benchmark/
        `-- temp/

The samples table
^^^^^^^^^^^^^^^^^

``samples.tsv`` is a tab-separated table with one row per sample and the columns
``sample_id``, ``r1``, ``r2``, ``qc_r1``, ``qc_r2`` and ``contigs``. It is
generated by :program:`virmake prep` and records the input files as absolute
paths. Normally there is no need to edit it; it can be edited by hand to mix
inputs per sample (for example reads for some samples and pre-assembled contigs
for others). Running :program:`virmake prep` with ``--force`` regenerates the
table and overwrites such edits.

Workflow configuration
^^^^^^^^^^^^^^^^^^^^^^

``config/params.yaml`` holds the workflow parameters (see :doc:`config-params`).
It is generated by :program:`virmake prep` and can be edited freely; ``--force``
regenerates it. A run always uses the configuration inside its own directory.
