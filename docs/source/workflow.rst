========
Workflow
========

.. warning::

    This page is under construction and may contain incorrect and/or incomplete
    information.

.. _header-workflow:

Overview
--------

.. figure:: /_static/workflow-overview.svg
    :alt: Workflow diagram
    :width: 80%
    :align: center

    Overview of the **VirMake** workflow.

Pre-processing & QC
-------------------

1) Trimming
^^^^^^^^^^^
..  NOTE! fastp has 1 thread hardcoded

If raw paired-end reads were specified during setup using
:option:`setup.py --reads`, `FASTP <https://github.com/OpenGene/fastp>`_ is
used for trimming using default options:

:--qualified_quality_phred: >=Q15
:--unqualified_percent_limit: 40
:--length_limit: 0
:--disable_adapter_trimming: False

Input file location is specified in <:confval:`path.input_reads`> (for more
info, see :doc:`config-params`).

Output files are stored in <:confval:`path.output`>/fastp_pe/.  For each
sample, paired trimmed reads are stored in ``{sample}_1.fastq`` and
``{sample}_2.fastq``. QC reports are stored in ``html`` and ``json`` format in
``{sample}.html`` and ``{sample}.json`` respectively.

2) Quality summary
^^^^^^^^^^^^^^^^^^

After quality trimming, or alternatively, if quality-trimmed reads were
provided using the :option:`setup.py --qc-reads` option,
`FastQC <https://github.com/s-andrews/FastQC>`_ is used to perform a
quality assessment. Finally, a summary is made using
`MultiQC <https://github.com/MultiQC/MultiQC>`_.

The output quality report is stored at
<:confval:`path.output`>/multiqc/multiqc.html.

Assembly
--------

After QC, the reads are assembled. The choice of assembler is specified in
:confval:`assembler`, though currently, only `MetaSpades
<https://github.com/ablab/spades>`_ is supported.  Assembly is performed with
default parameters and `MetaQuast <https://quast.sourceforge.net/metaquast>`_
is subsequently used with the ``--max-ref-number`` option set to ``0`` to
evaluate the assembly quality and generate the quality report.

Assemblies for each sample are stored in
<:confval:`path.output`>/metaSpades/``{sample}``/contigs.fasta. The assembly
quality report is kept at <:confval:`path.output`>/metaQUAST/report.html.

Viral Identification
--------------------

Viral sequences are then identified from the resulting contigs. Alternatively,
pre-assembled contigs can be provided to **VirMake** during setup with
:option:`setup.py --contigs`. **VirMake** supports two methods for
identification of viral sequences; `VirSorter2
<https://github.com/jiarong/VirSorter2>`_, or `geNomad
<https://github.com/apcamargo/genomad/>`_, as specified in
:confval:`identifier`.

..  NOTE! The min_contig_length parameters will be updated
   to 3000, then DRAM-v will use 1000 since viruses can be shorter than
   the contig they were identified in

Databases can be pre-downloaded using :option:`virmake db` and their location
is specified in :confval:`path.database.virsorter2` and
:confval:`path.database.genomad` for VirSorter2 and geNomad respectively.

When using ``VirSorter2`` for identification, output files will be stored in
<:confval:`path.output`>/virsorter/``{sample}``. This includes the predicted
viral sequences ``viruses.fasta``, the virus table ``virus_table.tsv``, and
viral boundries ``final-viral-combined.fa`` and ``final-viral-boundary.tsv``.

Similarly, when running ``geNomad``, output files will be stored in
<:confval:`path.output`>/virsorter/``{sample}``, including the predicted
viral sequences ``viruses.fasta``, the virus table ``virus_table.tsv``. Summary
files will be stored in
<:confval:`path.output`>/virsorter/``{sample}``/``{sample}_summary``/, including
``{sample}_virus.fna`` and ``{sample}_virus_summary.tsv``.

Identified viral sequences are then quality checked using `CheckV
<https://bitbucket.org/berkeleylab/CheckV>`_ and quality filtered to a minimum
of :confval:`quality_threshold` and quality filtered viral sequences are
stored at
<:confval:`path.output`>/virus_identification/``{sample}``/predicted_viruses.fasta".

Dereplication
-------------

..  NOTE! The virsorter_for_dram and checkv_vOTU_virsorter2 rules do not
   seem to run in the pipeline unless I missed a rule that calls them.

After viral identification, extracted viral genomes are clustered into vOTUs
using `Galah <https://github.com/wwood/galah>`_ with at average nucleotide
identity (ANI) :confval:`dereplication.ani`, precluster ANI
:confval:`dereplication.precluster_ani`, and with a minimal aligned fraction
:confval:`dereplication.min_aligned_fraction`. The clusters are stored in
<:confval:`path.output`>/dereplication/galah_clusters.tsv, representative
sequences are stored in
<:confval:`path.output`>/dereplication/repr_viral_seqs.fasta, and the
<:confval:`path.output`>/dereplication/old_to_new_ids.tsv table details which
sequences belong to which representative vOTU.

Mapping & Microdiversity
------------------------

.. NOTE! InStrain runs if :confval:`rule_inclusion.all.instrain` set to true
   NOTE! true by default!

.. NOTE! samtools view is run with a hardcoded 4 threads

If reads were provided, either as raw reads using :option:`setup.py --reads` or
as quality filtered reads using :option:`setup.py --qc-reads`, reads passing
the quality control are mapped to the representative viral vOTUs using
`Bowtie2 <https://github.com/BenLangmead/bowtie2>`_ and coverage is computed
using `BBTools pileup <https://github.com/bbushnell/BBTools>`_. A genome is
considered present if the coverage exceeds the specified :confval:`min_coverage`
threshold. The resulting relative abundance table is stored at
<:confval:`path.output`>/mapping/rel_abundance_table.tsv.

.. tip::
   InStrain can lead to high resource usage. Disable InStrain by setting
   :confval:`rule_inclusion.all.instrain` to false.

Strain-level diversity analysis is then conducted using `InStrain
<https://github.com/MrOlm/instrain>`_. The resulting diversity comparison is
stored at
<:confval:`path.output`>/instrain/comparison/output/comparison_comparisonsTable.tsv
