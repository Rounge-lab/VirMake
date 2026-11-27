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

2. Review ``testrun/samples.tsv`` and the configuration in
   ``testrun/config/``.

3. Preview the QC jobs:

.. code-block:: console

    virmake run testrun qc --dryrun

4. Run QC locally:

.. code-block:: console

    virmake run testrun qc --cores 4

5. Inspect the QC artifacts in ``testrun/results/output/qc/`` and the
   completion manifest at
   ``testrun/results/output/qc_manifest.json``.

.. note::

   The test dataset is supplied as quality-controlled reads. It exercises
   FastQC and artifact verification, but does not exercise raw-read trimming
   with fastp.

   QC does not require the viral reference databases. Software environments
   are created on first use, which can take longer than processing the small
   test dataset.

To try the complete analysis, first review the enabled stages and database
locations, then preview and execute the full workflow:

.. code-block:: console

    virmake run testrun all --dryrun
    virmake run testrun all --cores 8

.. warning::

   A full run can require substantial database downloads, disk space,
   memory and computation even when the input reads are tiny.


Running your own data
^^^^^^^^^^^^^^^^^^^^^

:program:`virmake prep` discovers samples from input filenames. Multiple
samples may be present in each input directory.

For paired-end reads, use either of these naming conventions:

.. code-block:: text

    SampleA_1.fastq.gz
    SampleA_2.fastq.gz

    SampleB_R1.fastq.gz
    SampleB_R2.fastq.gz

Uncompressed ``.fastq`` files are also accepted. Pre-assembled contigs use
the name ``{sample}.fasta``.

Choose the input option according to the processing already performed:

``--reads``
   Raw reads. VirMake runs fastp before assessing the processed reads with
   FastQC.

``--qc-reads``
   Reads that have already been preprocessed. VirMake skips fastp but still
   runs FastQC.

``--contigs``
   Pre-assembled contigs. VirMake uses these instead of assembling that
   sample's reads.

For example:

.. code-block:: console

    virmake prep myrun --reads myreads

Input options can be combined:

.. code-block:: console

    virmake prep myrun --qc-reads processed_reads --contigs assemblies

A sample may have both reads and contigs: the contigs are used for
identification, while the reads remain available for QC and mapping.
Contigs-only samples do not produce QC artifacts.

Review the generated ``samples.tsv`` before running the workflow. See
`The samples table`_ for explicit declarations, including single-end reads.

.. important::

   Do not declare the same sample as both raw and supplied-QC reads.
   Choose the read stage that describes the files you intend to analyse.


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

.. option:: --version

    Show the **VirMake** version and exit.

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

    Directory with paired-end raw reads. Accepted mate suffixes are
    ``_1/_2`` and ``_R1/_R2``, with ``.fastq`` or ``.fastq.gz`` extensions.
    These reads are processed with fastp.

.. option:: --qc-reads

    :Type: ``path``

    Directory with preprocessed paired-end reads, using the same naming
    conventions as ``--reads``. VirMake skips fastp for these samples
    but still assesses the supplied reads with FastQC.

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

    Regenerate the samples table and configuration in an existing run
    directory. ``results/`` is never touched.

.. warning::

   ``--force`` overwrites manual edits to ``samples.tsv``, workflow
   parameters and generated profiles. Back up those files first if you
   want to preserve custom settings.

   It does not clear previous results or guarantee a fresh analysis.
   Use a new run directory for an independent analysis.

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
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Before invoking Snakemake, VirMake checks that the samples table exists
and warns if it appears to contain no samples. It also performs basic
checks of the recorded environment prefix and selected database locations.

These checks are advisory, not a complete validation of every tool or
database. Snakemake constructs the requested workflow and determines its
actual dependencies.

Software environments are created when needed. Database setup jobs are
scheduled where required by the workflow. Use :program:`virmake envs`
and :program:`virmake db` to prepare resources ahead of time.

A header-only sample sheet can be useful while preparing a run, but QC
requires at least one read-bearing sample.

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

   Intermediate files marked as temporary can be removed as soon as their
   consumers finish, including during a QC-only run. A successful ``all``
   run also removes the workflow's ``results/temp/`` directory.

   Generated processed reads are retained only when ``keep_qc_reads`` is
   enabled. Supplied input reads are not deleted by this setting.

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

Both commands accept ``--threads`` to set Snakemake's core budget.
:program:`virmake db` additionally accepts ``-n`` to preview the jobs.

.. warning::

   The ``--clean`` options remove shared resources, not just resources
   belonging to one run. Other run directories may depend on the same
   environments or databases. Do not clean them while another analysis
   is using them.

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

The Snakemake controller remains in the session where the command was
started; individual workflow jobs are submitted to SLURM. Follow your
cluster's policy on whether the controller should run on a login node or
inside an allocation.

Environment solving and installation can happen in the controller before
jobs are submitted. SLURM submission does not necessarily move this work
onto a compute node. Consider preparing environments ahead of time in an
allocation if required by cluster policy.

The run directory, input files, databases and environment prefix must be
accessible from the compute nodes. Do not use a node-local temporary
directory as a shared environment prefix.

Keep the controller session alive until the workflow finishes, using
your cluster's recommended session or allocation mechanism.

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

``samples.tsv`` is a UTF-8, tab-separated declaration table with one row
per sample. Its header must contain exactly these seven columns:

.. code-block:: text

    sample_id	read_type	read_layout	read_stage	r1	r2	contigs

``sample_id``
   Unique sample identifier.

``read_type``
   Read technology category: ``short`` or ``long`` for read-bearing
   samples. The current version supports short reads only.

``read_layout``
   ``single`` for one read file, or ``paired`` for two mate files.

``read_stage``
   ``raw`` for reads requiring preprocessing, or ``qc`` for reads that
   have already been preprocessed.

``r1``
   Path to the first mate or the sole single-end read file.

``r2``
   Path to the second mate. Leave blank for single-end samples.

``contigs``
   Optional path to pre-assembled contigs. A sample can have contigs
   together with reads, or contigs alone.

For contigs-only samples, leave the read-related fields blank.

For example, the following declarations represent raw paired reads,
supplied-QC single reads, and a contigs-only sample. ``-`` means an empty
cell in this illustration; do not write literal hyphens into the sheet.

.. list-table::
   :header-rows: 1

   * - sample_id
     - read_type
     - read_layout
     - read_stage
     - r1
     - r2
     - contigs
   * - RawPair
     - short
     - paired
     - raw
     - /data/RawPair_1.fastq.gz
     - /data/RawPair_2.fastq.gz
     - -
   * - ProcessedSingle
     - short
     - single
     - qc
     - /data/ProcessedSingle.fastq.gz
     - -
     - -
   * - AssemblyOnly
     - -
     - -
     - -
     - -
     - -
     - /data/AssemblyOnly.fasta

:program:`virmake prep` writes absolute input paths. Prefer absolute
paths when editing declarations by hand.

.. important::

   The legacy ``qc_r1`` and ``qc_r2`` columns are no longer accepted.
   Put the supplied read paths in ``r1`` and ``r2`` and declare
   ``read_stage`` as ``qc``.

   Sample identifiers must be unique. QC artifact identities are
   normalized, so identifiers that become identical after normalization
   can cause verification-spec errors. Avoid names that differ only by
   case, such as ``SampleA`` and ``samplea``.

Changing ``read_stage`` changes how VirMake handles the files; it does
not preprocess them. Declare ``qc`` only when the supplied reads have
already undergone the preprocessing you require.

Running :program:`virmake prep` with ``--force`` regenerates the table
and overwrites manual edits.

Workflow configuration
^^^^^^^^^^^^^^^^^^^^^^

``config/params.yaml`` holds the workflow parameters (see :doc:`config-params`).
It is generated by :program:`virmake prep` and can be edited freely; ``--force``
regenerates it. A run always uses the configuration inside its own directory.


Resuming and recovering a run
-----------------------------

To resume an interrupted analysis, run the same command against the
same prepared directory:

.. code-block:: console

    virmake run myrun qc

Snakemake uses existing outputs and its execution metadata to determine
which jobs need to run. A normal rerun does not necessarily execute every
job again.

If Snakemake reports incomplete outputs, request their regeneration:

.. code-block:: console

    virmake run myrun qc --rerun-incomplete

If the run directory is locked, first confirm that no Snakemake
controller or submitted job is still writing to it. Only then remove a
stale lock:

.. code-block:: console

    virmake run myrun qc --unlock

.. warning::

   Do not use ``--nolock`` to bypass an active workflow. Two controllers
   writing to the same run directory can damage its state.
