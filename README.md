<h1 align="center">
  <picture>
    <source
      media="(prefers-color-scheme: dark)"
      srcset="docs/source/_static/virmake-logo-dark.svg">
    <img
      align="center"
      width="160"
      src="docs/source/_static/virmake-logo-light.svg"
      alt="VirMake">
  </picture>
  <br>
<\h1>

<h4 align="center">A <a href="https://snakemake.readthedocs.io/en/stable/">Snakemake</a> pipeline for viral metagenomic data analysis</h4>

<p align="center">
  <a href="https://github.com/Rounge-lab/VirMake/releases">
    <img src="https://img.shields.io/github/v/release/Rounge-lab/VirMake"
         alt="Latest release version">
  </a>
  <a href="https://virmake.readthedocs.io/en/latest/#">
    <img src="https://app.readthedocs.org/projects/virmake/badge/?version=latest"
         alt="Documentation status">
  </a>
  <a href="./LICENCE">
    <img src="https://img.shields.io/badge/license-GPL--v3.0-blue.svg"
         alt="Licence">
  </a>
</p>

<p align="center">
  <a href="#about-virmake">About VirMake</a> •
  <a href="#documentation">Documentation</a> •
  <a href="#quick-start">Quick start</a> •
  <a href="#license">License</a> •
  <a href="#citation">Citation</a> •
</p>


## About VirMake

**VirMake** is a [Snakemake](https://snakemake.readthedocs.io/en/stable/)-based
pipeline that offers viral metagenomic data analysis on paired-end data.  It
offers taxonomic and functional annotation, supports offline running and has
support for HPC cluster execution. It is made for Linux based systems and has
been tested on x86_64-based Linux.

## Documentation

Detailed documentation can be found on the
[VirMake Wiki page](https://rtd-virmake-test.readthedocs.io/en/latest/install.html).
Here you will find prerequisites, installation instructions, detailed usage
documentation, workflow description, and detailed configuration options. Please
refer to the Wiki page before opening issues.

Quick start
-----------

> [!NOTE]
> **Dependencies:** `git` `conda`

1. Clone the git repository

    ```bash
    git clone https://github.com/Rounge-lab/VirMake.git
    ```

2. Run the **VirMake** setup script with your input files (raw reads/QC
   reads/contigs)

    ```bash
    python setup.py [-h] [--work-dir WORK_DIR] [--reads READS] [--qc-reads
    QC_READS] [--contigs CONTIGS]
    ```

3. Set up the environments and databases

    ```bash
    virmake prep [OPTIONS]
    virmake db [OPTIONS]
    ```

4. Run the workflow

    ```bash
    virmake run [OPTIONS] <WORKFLOW>
    ```

## License

**VirMake** is provided under the GNU [General Public License
v3.0](https://www.gnu.org/licenses/gpl-3.0.en.html) (GPLv3). It is free
software, and you are welcome to redistribute and/or modify it under the
conditions of the GPLv3 terms. This software is provided 'As is' and comes with
ABSOLUTELY NO WARRANTY.

## Citation

When using **VirMake**, please cite our paper:

Birkeland E, Istvan P, Avershina E, Rognes T, Rounge TB.
*VirMake: a flexible and user-friendly pipeline for viral taxonomic and
functional characterisation from shotgun metagenomic sequencing data*.
bioRxiv. 2025. doi:
[https://doi.org/10.1101/2025.02.07.637044](https://doi.org/10.1101/2025.02.07.637044)
