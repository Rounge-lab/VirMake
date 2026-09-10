========
Workflow
========

.. _header-workflow:

Overview
--------

.. figure:: /_static/workflow-overview.svg
    :alt: Workflow diagram
    :width: 80%
    :align: center

    Overview of the **VirMake** workflow.

QC
--

Raw read pre-processing and quality control
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1) Trimming
"""""""""""
..  NOTE! fastp has 1 thread hardcoded

If raw paired-end reads were specified during setup using
:option:`setup.py --reads`, `FASTP <https://github.com/OpenGene/fastp>`_ is
used for trimming using default options:

:--qualified_quality_phred: >=Q15
:--unqualified_percent_limit: 40
:--length_limit: 0
:--disable_adapter_trimming: False

Input file location is specified in <:ref:`params.yaml input_reads
<params-input-reads>`> (for more info, see :ref:`config-workflow`).

Output files are stored in <:ref:`params.yaml output
<params-output>`>/fastp_pe/.  For each sample, paired trimmed reads are
stored in ``{sample}_1.fastq`` and ``{sample}_2.fastq``. QC reports are
stored in ``html`` and ``json`` format in ``{sample}.html`` and
``{sample}.json`` respectively.

2) Quality summary
""""""""""""""""""

After quality trimming, or alternatively, if quality-trimmed reads were
provided using the :option:`setup.py --qc-reads` option,
`FastQC <https://github.com/s-andrews/FastQC>`_ is used to perform a
quality assessment. Finally, a summary is made using
`MultiQC <https://github.com/MultiQC/MultiQC>`_.

The output quality report is stored at <:ref:`params.yaml output
<params-output>`>/multiqc/multiqc.html.
