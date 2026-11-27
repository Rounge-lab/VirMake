========
Workflow
========

.. _header-workflow:

Overview
--------

.. only:: html

   .. raw:: html
      :file: ../_static/workflow-overview.svg

.. only:: not html

   .. figure:: /_static/workflow-overview.svg
      :alt: Overview of the VirMake workflow
      :width: 80%
      :align: center

      Overview of the **VirMake** workflow.

Pre-processing & QC
-------------------

QC uses the read declarations in ``samples.tsv`` rather than inferring
processing state from output files or input-directory settings.

* Raw short reads are processed with fastp, then assessed with FastQC.
* Supplied-QC short reads skip fastp and are assessed directly with FastQC.
* Contigs-only samples have no QC artifacts.

QC supports single-end and paired-end short-read declarations. This does
not imply that every downstream analysis tool supports both layouts.

To run this stage on its own:

.. code-block:: console

    virmake run myrun qc --dryrun
    virmake run myrun qc --cores 4

Raw-read preprocessing
^^^^^^^^^^^^^^^^^^^^^^

`fastp <https://github.com/OpenGene/fastp>`_ processes raw reads using
the defaults of the installed tool version, apart from the workflow's
input, output, reporting and thread settings.

Single-end samples use one input file. Paired-end samples use both mates
in one fastp invocation.

Generated reads are written to:

.. code-block:: text

    results/output/reads/reads_qc/{sample}_1.fastq.gz
    results/output/reads/reads_qc/{sample}_2.fastq.gz # paired-end only

The fastp JSON report is retained at:

.. code-block:: text

    results/output/qc/fastp/{sample}.json

.. note::

   ``keep_qc_reads`` controls retention of generated processed reads.
   When disabled, those reads are temporary and may be removed after
   their consumers finish. Their filenames do not change.

   Reads supplied through ``--qc-reads`` remain at their declared input
   locations and are not removed by this setting.

FastQC assessment
^^^^^^^^^^^^^^^^^

`FastQC <https://github.com/s-andrews/FastQC>`_ assesses each declared
short-read file separately: processed reads for raw samples, or the
original supplied files for samples declared as ``qc``.

Persistent FastQC data files are stored at:

.. code-block:: text

    results/output/qc/fastqc/{sample}_1_fastqc_data.txt
    results/output/qc/fastqc/{sample}_2_fastqc_data.txt # paired-end only

Artifact verification and completion
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After producing the reports, VirMake verifies them with ``workflow-coda``.

``results/output/qc_spec.yaml`` declares the expected artifact families:

* One fastp report per raw short-read sample.
* One FastQC data file per short-read sample's first mate.
* One FastQC data file per paired short-read sample's second mate.

Successful verification writes:

.. code-block:: text

    results/output/qc_manifest.json

The manifest records artifact identities, paths, completion status,
checksums and verification records. Invalid artifacts fail the job rather than
being accepted as successful QC outputs.

.. important::

   Artifact verification checks that reports satisfy their declared
   format and verification requirements. It does not mean that every
   FastQC module passed or that the reads are suitable for every
   biological analysis. Review the quality metrics before interpreting
   downstream results.

Assembly
--------

After QC, reads are assembled. The assembler is selected with ``assembler``.
Assembly runs with default parameters, and `MetaQUAST
<https://quast.sourceforge.net/metaquast>`_ evaluates the result.

Assemblies are stored in ``results/output/metaSpades/{sample}/contigs.fasta``,
the assembly quality report in ``results/output/metaQUAST/report.html``.
Pre-assembled contigs provided with :option:`virmake prep --contigs` replace
this step for those samples.

Viral Identification
--------------------

Viral sequences are identified from contigs using `VirSorter2
<https://github.com/jiarong/VirSorter2>`_ or `geNomad
<https://github.com/apcamargo/genomad/>`_, selected with ``identifier``.  For
VirSorter2, the cutoffs are configured through ``virsorter2.id.*`` (see
:doc:`config-params`); databases are located via ``path.database.virsorter2``
and ``path.database.genomad``.

Identified sequences are quality checked with `CheckV
<https://bitbucket.org/berkeleylab/CheckV>`_ and filtered to the
``quality_threshold`` category. The filtered sequences are stored in
``results/output/virus_identification/{sample}/predicted_viruses.fasta``.

Dereplication
-------------

Viral genomes are clustered into vOTUs with `Galah
<https://github.com/wwood/galah>`_, using the thresholds in
``dereplication.ani``, ``dereplication.precluster_ani`` and
``dereplication.min_aligned_fraction``. Clusters are stored in
``results/output/dereplication/galah_clusters.tsv``, representative sequences in
``results/output/dereplication/repr_viral_seqs.fasta``, and
``results/output/dereplication/old_to_new_ids.tsv`` maps each sequence to its
representative vOTU. vOTU identifiers are configurable through
``dereplication.vOTU_prefix`` and related settings.

Mapping & Microdiversity
------------------------

If reads were provided, quality-controlled reads are mapped to the
representative vOTUs with `Bowtie2 <https://github.com/BenLangmead/bowtie2>`_
and coverage computed with `BBTools pileup
<https://github.com/bbushnell/BBTools>`_. A genome is considered present if its
coverage exceeds ``min_coverage``. The relative abundance table is stored at
``results/output/mapping/rel_abundance_table.tsv``.

Strain-level diversity is analysed with `InStrain
<https://github.com/MrOlm/instrain>`_; results are stored in
``results/output/instrain/comparison/output/comparison_comparisonsTable.tsv``.

.. tip::

   InStrain can lead to high resource usage. Disable it by setting
   ``rule_inclusion.all.instrain`` to ``false``.

Taxonomy
--------

Dereplicated vOTUs are classified with `vConTACT3
<https://bitbucket.org/MAVERICLab/vcontact3/src>`_ using the default
``prokaryotes`` ``--db-domain`` option. The database version and export formats
are set through ``vcontact3.*``. Assignments are stored in
``results/output/vcontact3/vOTU_assignments.csv``.

Functional Annotation
---------------------

Before annotation, vOTUs are re-run through VirSorter2 with the
``--prep-for-dramv`` option and more stringent cutoffs
(``virsorter2.for_dramv.min_length`` and ``virsorter2.for_dramv.min_score``),
and checked with CheckV. Annotation is done with `DRAM-v
<https://github.com/BortonWrightonLabs/DRAM>`_ using its default databases; the
minimum contig size is set by ``min_contig_size``. The output is stored in
``results/output/DRAMv/distilled/vMAG_stats.tsv``.

.. tip::

   Required database construction is resource-intensive. Disable functional
   annotation by setting ``rule_inclusion.all.function`` to ``false``.

Summary Statistics
------------------

Sample-wise and vOTU-wise summary tables are generated in
``results/output/statistics/sample_stats.tsv`` and
``results/output/statistics/vOTU_stats.tsv``.

Output Files
^^^^^^^^^^^^

The ``rel_abundance_table.tsv`` table shows the relative abundance of each
dereplicated vOTU per sample.

.. table:: Relative abundance table example
    :width: 75%
    :widths: auto
    :align: center

    +--------+-----------+-----------+-----------+----+
    |vOTU    |Sample1    |Sample2    |Sample3    |... |
    +========+===========+===========+===========+====+
    |vOTU001 |0.4        |0.1        |0          |    |
    +--------+-----------+-----------+-----------+----+
    |vOTU002 |0.6        |0.3        |0.1        |    |
    +--------+-----------+-----------+-----------+----+
    |vOTU003 |0          |0.2        |0          |    |
    +--------+-----------+-----------+-----------+----+
    |vOTU004 |0          |0.1        |0.2        |    |
    +--------+-----------+-----------+-----------+----+
    |vOTU005 |0          |0.3        |0          |    |
    +--------+-----------+-----------+-----------+----+
    |...     |           |           |           |    |
    +--------+-----------+-----------+-----------+----+


The ``sample_stats.tsv`` table shows genome, contig, and read-level statistics
per sample. The following columns are included:

.. describe:: Genome

    * ``n_genomes``
    * ``genome_median_length``
    * ``genome_min_length``
    * ``genome_max_length``
    * ``genomes_N50``
    * ``n_provirus``
    * ``checkv_complete``
    * ``checkv_high_quality``
    * ``checkv_medium_quality``
    * ``checkv_low_quality``

.. describe:: Contig

    * ``# contigs``
    * ``Largest contig``
    * ``Total length``
    * ``N50``
    * ``N90``
    * ``n_present``

.. describe:: Read

    * ``total_reads``
    * ``mapped_reads``
    * ``mapped_percent``

The ``vOTU_stats.tsv`` table shows statistics per dereplicated vOTU. The
following columns are included:

    * ``n_genomes``
    * ``vOTU_length``
    * ``vOTU_provirus``
    * ``pangenome_median_length``
    * ``pangenome_min_length``
    * ``pangenome_max_length``
    * ``pangenome_n_provirus``
    * ``Closer``
    * ``Accession``
    * ``Status``
    * ``VC``
    * ``Level``
    * ``Host``
    * ``BaltimoreGroup``
    * ``Realm``
    * ``Kingdom``
    * ``Phylum``
    * ``Class``
    * ``Order``
    * ``Family``
    * ``Subfamily``
    * ``Genus``
    * ``n_present``
    * ``mean_present``
    * ``Gene count``
    * ``Strand switches``
    * ``potential AMG count``
    * ``Transposase present``
    * ``Possible Non-Viral Contig``
