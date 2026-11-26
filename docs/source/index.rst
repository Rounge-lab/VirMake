=====================
VirMake Documentation
=====================


VirMake: a Snakemake pipeline for viral metagenomic data analysis
-----------------------------------------------------------------

**VirMake** is a `Snakemake <https://snakemake.readthedocs.io/en/stable/>`_
based pipeline that offers viral metagenomic data analysis on paired-end data.
It offers taxonomic and functional annotation, supports offline running and has
support for HPC cluster execution. It is made for Linux based systems and has
been tested on x86_64-based Linux.

License
-------

**VirMake** is provided under the GNU `General Public License v3.0
<https://www.gnu.org/licenses/gpl-3.0.en.html>`_ (GPLv3). It is free software,
and you are welcome to redistribute and/or modify it under the conditions of
the GPLv3 terms. This software is provided 'As is' and comes with ABSOLUTELY
NO WARRANTY.

Citation
--------

When using **VirMake**, please cite our paper:

Birkeland, E.E., Istvan, P., Avershina, E. et al. *VirMake: a flexible and
user-friendly pipeline for viral taxonomic and functional characterization from
shotgun metagenomic sequencing data.* **Microbiome** (2026).
`https://doi.org/10.1186/s40168-026-02582-8 <https://doi.org/10.1186/s40168-026-02582-8>`_

.. note::

   This project is under active development.

Contents
--------

.. toctree::
   :caption: User Documentation
   :maxdepth: 2

   content/install
   content/usage
   content/workflow
   content/config-params
   content/config-hpc

.. toctree::
   :caption: For Contributors
   :maxdepth: 2

   content/contributing
