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

..  NOTE! fastp has 1 thread hardcoded

If raw paired-end reads were specified during setup using
:option:`setup.py --reads`, `FASTP <https://github.com/OpenGene/fastp>`_ is used
for trimming and QC using default options

:`--qualified_quality_phred`: >=Q15
:`--unqualified_percent_limit`: 40
:`--length_limit`: 0
:`--disable_adapter_trimming`: False

Input file location is specified in :option:`params.yaml --input_reads` (for
more info, see :ref:`config-workflow`).

Output files are stored in :option:`params.yaml --output`/fastp_pe/. For each
sample, paired trimmed reads are stored in ``{sample}_1.fastq`` and
``{sample}_2.fastq``. QC reports are stored in ``html`` and ``json`` format in
``{sample}.html`` and ``{sample}.json`` respectively.
