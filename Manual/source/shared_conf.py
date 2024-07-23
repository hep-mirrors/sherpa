copyright = '2022, Sherpa Team'
author = 'Sherpa Team'

# The full version, including alpha/beta/rc tags
release = '3.0.0'


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
    ]
}

html_theme_options = {
    'logo': 'images/sherpa-logo.png',
    'logo_name': True,
    'logo_text_align': 'center',
}
