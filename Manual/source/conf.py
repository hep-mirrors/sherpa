# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#

import os
import sys
import re
sys.path.insert(0, os.path.abspath('.'))
sys.path.insert(0, os.path.abspath('./_extensions/'))

from shared_conf import *

# -- Project information -----------------------------------------------------

project = 'Sherpa Manual'

release = '[GIT]' # will be read from `configure.ac`

release = '3.0.0'
# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinxcontrib.bibtex',
    'gen_bash_completion'
]

# List of bibliography files, relative to the source directory
bibtex_bibfiles = ['manual/references.bib']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['manual/examples/*', 'manual/parameters/models/*']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

html_favicon = '_static/images/favicon.ico'

html_theme_options['extra_nav_links'] = {
    'Project Home': 'https://sherpa-team.gitlab.io',
    'Pre 3.0.0 Manuals': 'https://sherpa.hepforge.org/doc/',
}

# Remove logo name, because the logo already contains the string "SHERPA".
html_theme_options['logo_name'] = False

suppress_warnings = ['ref.option']

man_pages = [
    ('manpage', 'Sherpa', 'a Monte Carlo event generator for the Simulation of High-Energy Reactions of Particles ', 'Sherpa Team', 1)
]
