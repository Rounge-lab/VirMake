============
Installation
============

.. _prerequisites:


Prerequisites
-------------

Dependencies
^^^^^^^^^^^^

* `Git <https://git-scm.com/>`_
* Conda (either `Miniconda
  <https://continuumio-docs.readthedocs-hosted.com/miniconda/install/>`_ or
  `Anaconda <https://www.anaconda.com/docs/getting-started/installation>`_)

System Requirements
^^^^^^^^^^^^^^^^^^^

**VirMake** is designed for the *Linux* operating system.

The following OS versions have been explicitly tested:

* x86_64 Red Hat Enterprise Linux 9.8 (Plow)

Hardware
""""""""

* RAM: minimum 64 GB
* Storage: minimum 200 GB of free disk space
* Additional memory and disk space will be needed depending on the data to be
  analysed

.. _installation:

Installation
------------

.. _source_installation:

From Source
^^^^^^^^^^^

1. Clone the git repository

.. code-block:: console

    git clone https://github.com/Rounge-lab/VirMake.git

2. Run the **VirMake** setup script

.. program:: setup.py

Synopsis
""""""""

.. code-block:: console

    python setup.py [-h] [--work-dir WORK_DIR] [--reads READS] [--qc-reads
    QC_READS] [--contigs CONTIGS]

Options
"""""""

.. option:: -h, --help

    Show the help message and exit.

.. option:: --work-dir

    :Type: ``path``
    :Default: ``./results```

    Use this option to specify the output directory.

.. option:: --reads

    :Type: ``path``

    Path where reads can be found.

.. option:: --qc-reads

    :Type: ``path``

    Path where QC reads can be found.

.. option:: --contigs

    :Type: ``path``

    Path where contigs can be found.

Use :option:`--reads`, :option:`--qc-reads`, or :option:`--contigs` to
specify the directories where input files are found.

Reads specified with the :option:`--reads` option or the
:option:`--qc-reads` option must be paired-end reads and found in the
specified folder in files named either
``{sample}_1.fastq`` and ``{sample}_2.fastq`` (for uncompressed reads), or
``{sample}_1.fastq.gz`` and ``{sample}_2.fastq.gz`` (for gzip compressed
reads), where ``{sample}`` is the sample name. Multiple samples may be present.

If contigs are specified with the :option:`--contigs` option, they must be
found in the specified folder in files named ``{sample}.fasta``, where
``{sample}`` is the sample name. Multiple samples may be present.
