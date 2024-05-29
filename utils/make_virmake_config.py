import yaml
import pathlib
import click
import os


# CLI command tool to choose the viral identifier

def make_config(virmake_path,
                work_dir,
                input_reads,
                input_contigs):
    """Creates a default config structure."""
    config = {}
    db_path = virmake_path / "resources" / "databases"

    config["slurm_account"] = "default"
    config["assembler"] = "metaspades"
    config["identifier"] = "virsorter2" # may be changed to vibrant or genomad
    config["trim_percentage"] = 0.05
    config["min_coverage"] = 75
    config["min_contig_size"] = 1000

    config["vibrant"] = {
        "is_virome": "no",
        }  
           
    config["virsorter2"] = {
        "id": {
            "min_length": 3000,
            "min_score": 0.5,
            "viral_groups": "dsDNAphage,ssDNA,NCLDV,RNA,lavidaviridae",
        },
        "for_dramv": {
            "min_length": 100,
            "min_score": 0.01,
            "viral_groups": "dsDNAphage,ssDNA,NCLDV,RNA,lavidaviridae",
        },
    }
 

    config["path"] = {
        "virmake": str(virmake_path),
        "envs": str(virmake_path / "workflow" / "envs"),
        "input": str(virmake_path / work_dir / "input"),
        "output": str(virmake_path / work_dir / "output"),
        "log": str(virmake_path / work_dir / "log"),
        "benchmark": str(virmake_path / work_dir / "benchmark"),
        "temp": str(virmake_path / work_dir / "temp"),
        "scripts": str(virmake_path / "workflow" / "scripts"),
        "samples": str(virmake_path / work_dir / "samples.tsv"),
        "input_reads": str(os.path.normpath(input_reads)),
        "input_contigs": str(os.path.normpath(input_contigs)),
        "database": {
            "RefSeq": str(db_path / "RefSeq/viral.1.1.genomic.fna"),
            "virsorter2": str(db_path / "virsorter2"),
            "genomad": str(db_path / "genomad"),
            "vibrant": str(db_path / "vibrant" / "vibrant-1.2.1"),
            "checkv": str(db_path / "checkv"),
            "INPHARED": str(db_path / "INPHARED"),
            "vcontact2": str(db_path / "vcontact2"),
            "DRAM": str(db_path / "DRAM" / "DRAM_data"),
        },
    }

    config["dereplication"] = {
        "ani": 97,
        "precluster_ani": 95,
        "min_aligned_fraction": 70,
        "vOTU_num_start": 1,
        "vOTU_prefix": "vOTU",
        "vOTU_suffix": "",
        "vOTU_num_len": 5
    }
    config["quality_threshold"] = "medium"
    config["threads"] = 24
    config["job_type"] = {
        "small": "normal",
        "normal": "normal",
        "big": "bigmem",
    }
    config["memory"] = {
        "tiny": 1000,
        "small": 8000,
        "normal": 16000,
        "big": 32000,
        "vcontact2": 63000,
        "metaquast": 63000,
        "metaspades": 63000,
    }
    config["time"] = {
        "tiny": "30 m",
        "small": "1 h",
        "normal": "6 h",
        "big": "13 h",
        "vcontact2": "72 h",
        "metaquast": "24 h",
        "metaspades": "2 h",
    }
    config["include_tables"] = [
        "raw_coverage_table.tsv",
        "trimmed_mean_coverage.tsv",
        "amg_summary.tsv",
        "quality_summary.tsv",
        "quality_summary.tsv",
        "quality_summary.tsv",
        "contig_plasmid_summary",
        "results_vcontact2",
        "comparison_comparisonsTable.tsv",
        "transposed_report.tsv",
        "viral_genomes_combined.csv",
        "derep95_combined.fasta.clstr",
    ]
    config["rule_inclusion"] = {
        "all": {
            "qc": True,
            "assembly": True,
            "metaquast": True,
            "identification": True,
            "taxonomy": True,
            "mapping": True,
            "instrain": True,
            "function": True,
            "stats": True,
        },
        "stats": {
            "metaquast": True,
            "mapping": True,
            "instrain": False,
            "dramv": True,
        }
    }

    return config

@click.command()
@click.argument("virmake_path")
@click.option(
    "-w",
    "--work-dir",
    default="results",
    help="Path to working directory.",
)
@click.option(
    "-r",
    "--reads",
    default="",
    help="Location of read files",
)
@click.option(
    "-c",
    "--contigs",
    default="",
    help="Location of assembled contigs.",
)
def run_setup(virmake_path, work_dir, reads, contigs):
    """Saves the config to a file."""
    virmake_path = pathlib.Path(virmake_path)
    config_path = virmake_path / "config" / "params.yaml"
    config = make_config(virmake_path, 
                         work_dir=work_dir, 
                         input_reads=reads, 
                         input_contigs=contigs)
    if config_path.exists():
        print(f"Config file {config_path} already exists.")
        print("Skipping...")
    else:
        with open(config_path, "w+") as f:
            yaml.dump(config, f)


if __name__ == "__main__":
    run_setup()
