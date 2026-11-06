"""VirMake assets are now bundled with the python package to make the install read-only
files are locations are resolved through :func:`importlib.resources.files`"""

from importlib.resources import files
from pathlib import Path


def package_root() -> Path:
    return Path(str(files("virmake")))


def workflow_dir() -> Path:
    """Dir with Snakefile and rule files"""
    return package_root() / "workflow"


def snakefile() -> Path:
    return workflow_dir() / "Snakefile"


def envs_dir() -> Path:
    """per-rule conda environments"""
    return workflow_dir() / "envs"


def scripts_dir() -> Path:
    """workflow scripts"""
    return workflow_dir() / "scripts"


def dummy_resources_dir() -> Path:
    """Dir for dummy TSV templates referenced by rules"""
    return package_root() / "resources" / "dummy"


def templates_dir() -> Path:
    """Dir with config templates for ``virmake prep``"""
    return package_root() / "templates"
