=============
Configuration
=============
.. ./config/params.yaml

.. _config-workflow:

Workflow parameters
-------------------

These workflow parameters are stored in ``config/params.yaml`` unless another
parameters file is provided using the :option:`virmake run --config-file`
option.

..  NOTE! Is there actually another supported assembler choice?

.. _params-assembler:

.. option:: assembler

    :Type: ``str``
    :Default: metaspades

.. _params-identifier:

.. option:: identifier

    :Type: ``str``
    :Default: virsorter2
    :Alternative: genomad

.. option:: min_contig_size

    :Type: ``int``
    :Default: 1000

.. option:: min_coverage

    :Type: ``int``
    :Default: 75

.. option:: quality_threshold

    :Type: ``str``
    :Default: medium

.. option:: slurm_account

    :Type: ``str``
    :Default: default

.. option:: threads

    :Type: ``int``
    :Default: 24

.. option:: trim_percentage

    :Type: ``float``
    :Default: 0.05

Dereplication:
^^^^^^^^^^^^^^

.. option:: ani

    :Type: ``int``
    :Default: 97

.. option:: min_aligned_fraction

    :Type: ``int``
    :Default: 70

.. option:: precluster_ani

    :Type: ``int``
    :Default: 95

.. option:: vOTU_num_len

    :Type: ``int``
    :Default: 5

.. option:: vOTU_num_start

    :Type: ``int``
    :Default: 1

.. option:: vOTU_prefix

    :Type: ``str``
    :Default: vOTU

.. option:: vOTU_suffix

    :Type: ``str``
    :Default: ''

job_type:
^^^^^^^^^

.. option:: big

    :Type: ``str``
    :Default: bigmem

.. option:: normal

    :Type: ``str``
    :Default: normal

.. option:: small

    :Type: ``str``
    :Default: normal

memory:
^^^^^^^
..  big: 32000$
    metaquast: 63000$
    metaspades: 63000$
    normal: 16000$
    small: 8000$
    tiny: 1000$
    vcontact2: 63000$

.. option:: big:

    :Type: ``str``
    :Default: 32000

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

.. _params-input-reads:

.. option:: input_reads

    :Type: ``path``

    This parameter is set when paired-end reads were specified during setup
    using :option:`setup.py --reads`.

.. option:: log

    :Type: ``dir``
    :Default: ``./results/log/``

.. _params-output:

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

path.database:
""""""""""""""

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

rule_inclusion.all:
"""""""""""""""""""

.. option:: assembly

   :Type: ``bool``
   :Default: True

.. option:: function

   :Type: ``bool``
   :Default: True

.. option:: identification

   :Type: ``bool``
   :Default: True

.. option:: instrain

   :Type: ``bool``
   :Default: True

.. option:: mapping

   :Type: ``bool``
   :Default: True

.. option:: metaquast

   :Type: ``bool``
   :Default: True

.. option:: qc

   :Type: ``bool``
   :Default: True

.. option:: stats

   :Type: ``bool``
   :Default: True

.. option:: taxonomy

   :Type: ``bool``
   :Default: True

rule_inclusion.stats:
"""""""""""""""""""""

.. option:: dramv

   :Type: ``bool``
   :Default: True

.. option:: instrain

   :Type: ``bool``
   :Default: false

.. option:: mapping

   :Type: ``bool``
   :Default: True

.. option:: metaquast

   :Type: ``bool``
   :Default: True


.. _config-hpc:

HPC configuration
-----------------

