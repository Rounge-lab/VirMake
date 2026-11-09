============
Installation
============

.. _prerequisites:

Prerequisites
-------------

Dependencies
^^^^^^^^^^^^

* Conda (either `Miniconda
  <https://docs.anaconda.com/miniconda/>`_ or `Anaconda
  <https://www.anaconda.com/docs/getting-started/installation>`_)

System Requirements
^^^^^^^^^^^^^^^^^^^

**VirMake** is designed for the *Linux* operating system.
The following OS versions have been explicitly tested:

* x86_64 Red Hat Enterprise Linux 9.8 (Plow)

Hardware
""""""""

* RAM: minimum 64 GB
* Storage: minimum 200 GB of free disk space
* Additional memory and disk space will be needed depending on the data to
  be analysed

.. _installation:

Installation
------------

**VirMake** can be installed one of in three ways. The workflow's tool
environments (fastp, metaSPAdes, VirSorter2, ...) are downloaded and created
through conda automatically the first time **VirMake** is run (see
:doc:`usage`).

.. _bioconda_installation:

From Bioconda
^^^^^^^^^^^^^

This is the recommended way to install **VirMake**:

.. code-block:: console

    conda create -n virmake -c conda-forge -c bioconda virmake
    conda activate virmake
    virmake --help

The package installs **VirMake** together with snakemake and its Python
dependencies.

.. note::

   The bioconda package becomes available with the **VirMake** v1.1.0
   release. Until then, use one of the two methods below.

.. _pip_methods:


From a GitHub Release
^^^^^^^^^^^^^^^^^^^^^

.. note::

   Requires ``Python>=3.12``

Each release publishes a source archive on the `GitHub releases page
<https://github.com/Rounge-lab/VirMake/releases>`_.

1. Download the source archive (``virmake-<version>.tar.gz``) from the
   latest release.

2. Install it with pip:

.. code-block:: console

    python -m pip install virmake-<version>.tar.gz
    virmake --help


.. _source_installation:

From Source
^^^^^^^^^^^

Install the newest development state directly from the git repository:

.. note::

   Requires ``Git``

1. Clone the git repository

.. code-block:: console

    git clone https://github.com/Rounge-lab/VirMake.git

2. Install **VirMake**:

.. code-block:: console

    conda create -n virmake python=3.12
    conda activate virmake
    python -m pip install ./VirMake
    virmake --help


The source installation tracks whatever branch you cloned. Contributors and
developers who want a live checkout should follow the development installation
in :doc:`contributing`.
