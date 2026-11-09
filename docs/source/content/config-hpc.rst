HPC Configuration
=================

.. warning::

    This page is under construction and may contain incorrect and/or incomplete
    information.


..  NOTE! havin this be --profile and the workflow params option being called
   --config-file seems a little confusing.
   NOTE! why is the default value for --profile the directory ./config/ when
   the default value for --config-file is the file ./config/params.yaml?

These HPC parameters are stored in ``config/config.yaml`` unless another
parameters file is provided using the :option:`virmake run --profile`
option.

.. yaml-config:: _data/config.yaml
