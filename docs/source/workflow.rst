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
default parameters and `Metaquast <https://quast.sourceforge.net/metaquast>`_
is subsequently used with the ``--max-ref-number`` option set to ``0`` to
evaluate the assembly quality and generate the quality report.

Assemblies for each sample are stored in
<:confval:`path.output`>/metaSpades/``{sample}``/contigs.fasta. The assembly
quality report is kept at <:confval:`path.output`>/metaQUAST/report.html.

Viral Identification
--------------------

Viral sequences are then identified from the resulting contigs. Alternatively,
pre-assembled contigs can be provided to **VirMake** during setup with
:option:`setup.py --contigs`.

