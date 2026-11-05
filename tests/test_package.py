"""Version and packaged workflow assets in place"""

import importlib.resources

import virmake


def test_version():
    assert virmake.__version__ == "1.1.0"


def test_workflow_assets_present():
    pkg = importlib.resources.files("virmake")
    assert (pkg / "workflow" / "Snakefile").is_file()
    assert (pkg / "workflow" / "envs" / "DRAMv.post-deploy.sh").is_file()
    assert (pkg / "resources" / "dummy" / "quality_summary_dummy.tsv").is_file()
    assert (pkg / "templates" / "params.yaml").is_file()
    assert (pkg / "templates" / "workflow-profile.yaml").is_file()
    assert (pkg / "templates" / "slurm-profile.yaml").is_file()
