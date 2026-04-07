# Sherpa Manual
This is the source for the sherpa manual, powered by the [sphinx
documentation generator](http://www.sphinx-doc.org/en/).

The online Version is available under: https://sherpa-team.gitlab.io/sherpa/index.html

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

## Building the Docs
 - run `cmake` with `-DSHERPA_ENABLE_MANUAL=ON`
 - run `make` in the Manual directory to build all targets

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
The `release` option in `conf.py` may not be set, as it is automatically set in the sphinx call in `CMakeLists.txt`.

## Distribution
The manual will be automatically deployed to Gitlab pages through the `manual` and `manual-index` jobs.
