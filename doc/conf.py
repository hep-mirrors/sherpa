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

copyright = '2026, Sherpa Team'
author = 'Sherpa Team'

# -- Project information -----------------------------------------------------

project = 'Sherpa Manual'

release = '[GIT]' # will be defined from CMakeLists.txt

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinxcontrib.bibtex',
    'gen_bash_completion'
]

# List of bibliography files, relative to the source directory
bibtex_bibfiles = ['references.bib']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['examples/*', 'parameters/models/*']

templates_path = ["_templates"]

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'alabaster'
html_sidebars = {
    '**': [
        'about.html',
        'navigation.html',
        'relations.html',
        'searchbox.html',
        'donate.html',
        "singlemulti.html"
    ]
}

html_theme_options = {
    'logo': 'images/sherpa-logo.png',
    'logo_name': False, # Remove logo name, because the logo already contains the string "SHERPA"
    'logo_text_align': 'center',
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

html_favicon = '_static/images/favicon.ico'

suppress_warnings = ['ref.option']

man_pages = [
    ('manpage', 'Sherpa', 'a Monte Carlo event generator for the Simulation of High-Energy Reactions of Particles ', 'Sherpa Team', 1)
]
