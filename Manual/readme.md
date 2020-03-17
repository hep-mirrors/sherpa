<!-- If you have questions or requests ask Valentin Boettcher <hiro@protagon.space> (Gitlab: vale9811) -->

# Sherpa Manual
This is the source for the sherpa manual, powered by the [sphinx
documentation generator](http://www.sphinx-doc.org/en/).

The online Version is available under: https://sherpa-team.gitlab.io/sherpa/

## Important Directories
 - `source/` the source of the
   - `source/index.rst` the root document
   - `source/manpage.rst` the manpage source
 - `source/manual` the sources of the manual contents
 - `source/man` the sources of the extra manpage contents

## Conventions
 - please use the first heading of the contents as filename:
   - spaces -> `-`
   - downcase
   - extension: `.rst`

## Build Dependencies
 - `python 3`
 - `sphinx >= 2.2.0`
 - `sphinxcontrib-bibtex`
 - optional:
   - `makeinfo` to build the info manual
   - a LaTeX distribution to build the pdf manual; usually a standard
     texlive installation should do (from texlive directly not the
     distribution packages)

See also `requirements.txt` for pinned versions of the python packages.


## Building the Docs
 - run `configure` with `--enable-manual`
 - run `make` in the Manual directory to build all targets
   - run `make sherpamanual_html` to build the html manual
   - run `make sherpamanual.pdf` to build the pdf manual
   - run `make sherpamanual.info` to build the info manual
   - run `make Sherpa.1` to build the manpage

Note that the manual is always built as if you specified `-j1` to
prevent race conditions.

## Caveats
If you see something like:
```rst
:ref:`text <text>`
```
it can be replaced by:
```rst
:ref:`text`
```
if the corresponding label refers to a heading!

## Completion Index
After building the html docs you can generate bash completion indices
by running: `make completion.index`.
This will create `completion.index` and `options.index`.

## Version/Release Number
The `release` option in `conf.py` may not be set, as it is automatically read from `../configure.ac`.

## Distribution
The tarball created by `make-dist` will include the rendered manual in
all formats if it was configured with `--enable-manual`.
