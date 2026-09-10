=============
Configuration
=============
.. ./config/params.yaml

.. _config-workflow:

Workflow parameters
-------------------

.. program:: params.yaml

path:
^^^^^

.. option:: benchmark

    :Type: ``dir``
    :Default: ``./results/benchmark/``


.. option:: envs

   :Type: ``dir``
   :Default: ``./workflow/envs/``

    Folder containing the environment.yaml files for each tool included in the
    **VirMake** workflow.

.. option:: input_contigs

    :Type: ``path``

    This parameter is set when contigs were specified during setup using
    :option:`setup.py --contigs`.


.. option:: input_reads

    :Type: ``path``

    This parameter is set when paired-end reads were specified during setup
    using :option:`setup.py --reads`.

.. option:: log

    :Type: ``dir``
    :Default: ``./results/log/``

.. option:: output

    :Type: ``dir``
    :Default: ``./results/output/``

.. option:: samples

    :Type: ``path``
    :Default: ``./results/samples.tsv``

    This is a tab separated file containing the following columns:
    - sample_id
    - r1
    - r2
    - qc_r1
    - qc_r2
    - contigs

.. option:: scripts

    :Type: ``dir``
    :Default: ``./workflow/scripts/``

    This directory contains scripts needed to run the **VirMake** workflow.

.. option:: temp

    :Type: ``dir``
    :Default: ``./results/temp/``

    This directory is used to store temporary files when running the
    **VirMake** workflow.

.. option:: virmake

    :Type: ``dir``
    :Default: ``./``

    The **VirMake** repository root directory.

database:
"""""""""

.. option:: DRAM

    :Type: ``path``
    :Default: ``./resources/databases/DRAM/DRAM_data``

.. option:: INPHARED

    :Type: ``path``
    :Default: ``./resources/databases/INPHARED``

.. option:: RefSeq

    :Type: ``path``
    :Default: ``./resources/databases/RefSeq/viral.q.q.genomic.fna``

.. option:: checkv

    :Type: ``path``
    :Default: ``./resources/databases/checkv``

.. option:: genomad

    :Type: ``path``
    :Default: ``./resources/databases/genomad``

.. option:: vcontact2

    :Type: ``path``
    :Default: ``./resources/databases/vcontact2``

.. option:: vibrant

    :Type: ``path``
    :Default: ``./resources/databases/vibrant/vibrant-1.2.1``

.. option:: virsorter2

    :Type: ``path``
    :Default: ``./resources/databases/virsorter2``


rule_inclusion:
^^^^^^^^^^^^^^^

all:
""""

.. confval:: rule_inclusion.all.assembly
   :type: bool
   :default: true

.. confval:: rule_inclusion.all.function
   :type: bool
   :default: true

.. confval:: rule_inclusion.all.identification
   :type: bool
   :default: true

.. confval:: rule_inclusion.all.instrain
   :type: bool
   :default: true

.. confval:: rule_inclusion.all.mapping
   :type: bool
   :default: true

.. confval:: rule_inclusion.all.metaquast
   :type: bool
   :default: true

.. confval:: rule_inclusion.all.qc
   :type: bool
   :default: true

.. confval:: rule_inclusion.all.stats
   :type: bool
   :default: true

.. confval:: rule_inclusion.all.taxonomy
   :type: bool
   :default: true

.. confval:: rule_inclusion.stats.dramv
   :type: bool
   :default: true

.. confval:: rule_inclusion.stats.instrain
   :type: bool
   :default: false

.. confval:: rule_inclusion.stats.mapping
   :type: bool
   :default: true

.. confval:: rule_inclusion.stats.metaquast
   :type: bool
   :default: true


.. _config-hpc:

HPC configuration
-----------------
