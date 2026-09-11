======================
Workflow Configuration
======================

.. _config-workflow:

Workflow parameters
-------------------

These workflow parameters are stored in ``config/params.yaml`` unless another
parameters file is provided using the :option:`virmake run --config-file`
option.

.. yaml-config:: _data/params.yaml

..  NOTE! Is there actually another supported assembler choice?

    envs
    Folder containing the environment.yaml files for each tool included in the
    **VirMake** workflow.

    input_contigs
    This parameter is set when contigs were specified during setup using
    :option:`setup.py --contigs`.

    input_reads
    This parameter is set when paired-end reads were specified during setup
    using :option:`setup.py --reads`.

    samples
    This is a tab separated file containing the following columns:
    - sample_id
    - r1
    - r2
    - qc_r1
    - qc_r2
    - contigs

    scripts
    This directory contains scripts needed to run the **VirMake** workflow.

    temp
    This directory is used to store temporary files when running the
    **VirMake** workflow.

    virmake
    The **VirMake** repository root directory.

