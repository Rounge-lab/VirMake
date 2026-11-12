<h1 align="center">
  <picture>
    <source
      media="(prefers-color-scheme: dark)"
      srcset="docs/source/_static/virmake_logo_dark.svg">
    <img
      align="center"
      width="160"
      src="docs/source/_static/virmake_logo_light.svg"
      alt="VirMake">
  </picture>
  <br>
</h1>

<h4 align="center">A <a href="https://snakemake.readthedocs.io/en/stable/">Snakemake</a> pipeline for viral metagenomic data analysis</h4>

<p align="center">
  <a href="https://github.com/Rounge-lab/VirMake/releases">
    <img src="https://img.shields.io/github/v/release/Rounge-lab/VirMake"
         alt="Latest release version">
  </a>
  <a href="https://github.com/Rounge-lab/VirMake/actions/workflows/ci.yml">
    <img src="https://github.com/Rounge-lab/VirMake/actions/workflows/checks.yml/badge.svg"
         alt="CI status">
  </a>
  <a href="https://virmake.readthedocs.io/en/latest/">
    <img src="https://readthedocs.org/projects/virmake/badge/?version=latest"
         alt="Documentation status">
  </a>
  <a href="./LICENSE">
    <img src="https://img.shields.io/badge/license-GPL--v3.0-blue.svg"
         alt="Licence">
  </a>
</p>

<p align="center">
  <a href="#about-virmake">About VirMake</a> •
  <a href="#documentation">Documentation</a> •
  <a href="#quick-start">Quick start</a> •
  <a href="#license">License</a> •
  <a href="#citation">Citation</a>
</p>

## About VirMake

**VirMake** is a [Snakemake](https://snakemake.readthedocs.io/en/stable/)-based
pipeline that offers viral metagenomic data analysis on paired-end data.  It
offers taxonomic and functional annotation, supports offline running and has
support for HPC cluster execution. It is made for Linux based systems and has
been tested on x86\_64-based Linux; workflow execution is not supported on
Windows or macOS.

## Documentation

Detailed documentation is hosted at
[virmake.readthedocs.io](https://virmake.readthedocs.io/en/latest/).
Here you will find prerequisites, installation instructions, detailed usage
documentation, workflow description, and detailed configuration options.
Please refer to the documentation before opening issues.

## Quick start

> [!NOTE]
> **Dependencies:** conda (either
> [Miniconda](https://docs.anaconda.com/miniconda/) or
> [Anaconda](https://www.anaconda.com/docs/getting-started/installation))

1. Install **VirMake**:

    ```bash
    conda create -n virmake -c conda-forge -c bioconda virmake
    conda activate virmake
    ```

    > The bioconda package becomes available with the v1.1.0 release.
    > Until then, install from the latest [GitHub release](https://github.com/Rounge-lab/VirMake/releases)
    > (see [virmake.readthedocs.io/install](https://virmake.readthedocs.io/en/latest/install.html)).

2. Try the bundled synthetic dataset:

    ```bash
    virmake prep testrun --test-data
    virmake run testrun all --cores 8
    ```

3. Run your own data (paired-end reads, recognised by file name. Please note
   the naming conventions in the documentation):

    ```bash
    virmake prep myrun --reads myreads
    virmake run myrun all --cores 8
    ```

4. On a Slurm cluster:

    ```bash
    virmake run myrun all --slurm --account myproject
    ```

## License

**VirMake** is provided under the GNU [General Public License
v3.0](https://www.gnu.org/licenses/gpl-3.0.en.html) (GPLv3). It is free
software, and you are welcome to redistribute it and/or modify it under the
conditions of the GPLv3 terms. This software is provided 'As is' and comes with
ABSOLUTELY NO WARRANTY.

## Citation

When using **VirMake**, please cite our paper:

Birkeland E, Istvan P, Avershina E, Rognes T, Rounge TB.
*VirMake: a flexible and user-friendly pipeline for viral taxonomic and
functional characterisation from shotgun metagenomic sequencing data*.
bioRxiv. 2025. doi:
[https://doi.org/10.1101/2025.02.07.637044](https://doi.org/10.1101/2025.02.07.637044).
