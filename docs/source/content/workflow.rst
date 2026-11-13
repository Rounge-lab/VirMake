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

1) Trimming
^^^^^^^^^^^

If raw paired-end reads were provided during preparation using :option:`virmake
prep --reads`, `FASTP <https://github.com/OpenGene/fastp>`_ is used for trimming
with default options:

* ``--qualified_quality_phred`` >= Q15
* ``--unqualified_percent_limit`` 40
* ``--length_limit`` 0
* ``--disable_adapter_trimming`` False

The input directory is recorded in ``path.input_reads`` by :program:`virmake
prep`.

Trimmed reads are stored in ``results/output/fastp_pe/`` as
``{sample}_1.fastq.gz`` and ``{sample}_2.fastq.gz``. QC reports are stored in
``{sample}.html`` and ``{sample}.json``.

2) Quality summary
^^^^^^^^^^^^^^^^^^

After quality trimming, or alternatively if quality-controlled reads were
provided using :option:`virmake prep --qc-reads`, `FastQC
<https://github.com/s-andrews/FastQC>`_ performs a quality assessment,
summarised with `MultiQC <https://github.com/MultiQC/MultiQC>`_. The quality
report is stored at ``results/output/multiqc/multiqc.html``.

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

Report Generation
-----------------

Sample-wise and vOTU-wise summary reports are generated in
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
