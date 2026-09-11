# Configuration file for the Sphinx documentation builder.
import os
import sys

sys.path.insert(0, os.path.abspath("_ext"))

# -- Project information

project = 'VirMake'
copyright = '2023, CRCbiome'

release = '0.1'
version = '0.1.0'

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinx_copybutton',
    'sphinx_toolbox.confval',
    'yaml_config',
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
html_static_path = ["_static"]

html_css_files = [
    "virmake.css",
]

html_logo = "virmake_emblem_white.svg"
html_favicon = "virmake_favicon.svg"
