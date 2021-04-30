#!/bin/sh

###
# This script builds the manual for all revisions and in multiple formats.
# The final page will reside in $OUTDIR/_multi
#
# Author: Valentin Boettcher <hiro at protagon dot space>
##

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
OUTDIR=$SCRIPT_DIR/$1
MV_CONFDIR=Manual/.multi-version-config/

cd $SCRIPT_DIR/..

if [ "$#" -lt "1" ]; then
   echo "Usage: $0 [build directory]"
   exit 1
fi

## build versioned in multiple formats
sphinx-multiversion Manual/source/ $OUTDIR/_multi -- -c $MV_CONFDIR -A mode=html -A pdf=sherpamanual.pdf \
                    & sphinx-multiversion Manual/source/ $OUTDIR/_single -- -c $MV_CONFDIR -b singlehtml -A mode=singlehtml -A pdf=sherpamanual.pdf \
                    & sphinx-multiversion Manual/source/ $OUTDIR/_pdf -- -c $MV_CONFDIR -b latex

## rename single_page html files
cd $OUTDIR
find _single/ -iname "index.html" -exec rename .html _single.html '{}' \;

## copy single-page html files to the mult-page folder
rsync -mirv --include "*/"  --include="*.html" --exclude "*"  _single/* _multi/

## make the pdfs
find _pdf -iname "makefile" -exec sh -c 'cd $(dirname {}) && make'  \;

## copy single-page html files to the mult-page folder
rsync -mirv --include "*/"  --include="*.pdf" --exclude "*"  _pdf/* _multi/
