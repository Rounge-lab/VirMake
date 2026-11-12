# Configuration file for the Sphinx documentation builder.
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent / "_ext"))

from importlib.metadata import PackageNotFoundError, version as _pkg_version

# -- Project information

project = 'VirMake'
copyright = '2023, CRCbiome'

# Version from the installed package (single-sourced from pyproject).
try:
    release = _pkg_version("virmake")
except PackageNotFoundError:
    release = "0.0.0"
version = ".".join(release.split(".")[:2])

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinx_copybutton',
    'sphinx_toolbox.confval',
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'

# -- Options for EPUB output
epub_show_urls = 'footnote'

# -- Sphinx copybutton
copybutton_prompt_text = r">>> |\$ "
copybutton_prompt_is_regexp = True

# -- Styling
html_theme_options = {
        'logo_only': True,
}

html_static_path = ["_static"]

html_logo = "_static/virmake_logo_dark.svg"
html_favicon = "_static/virmake_favicon.svg"

html_context = {
    "display_github": True,
    "github_repo": "Rounge-lab/VirMake",
}

html_css_files = [
    "crcbiome-rtd.css",
]

html_js_files = [
    "dark-mode-toggle.js",
]
