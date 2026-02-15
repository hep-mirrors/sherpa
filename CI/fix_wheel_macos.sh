#!/bin/bash
# Fix absolute Fortran runtime paths in macOS wheels.
#
# On macOS, gfortran bakes absolute paths to libgfortran, libquadmath,
# and libgcc_s into shared libraries. This script rewrites those to
# @rpath references and adds common Homebrew RPATH entries so the wheel
# is relocatable. Users must have gfortran installed on their system.
#
# Usage:
#   ./fix_wheel_macos.sh dist/*.whl
#
# Produces a fixed wheel in dist/fixed_wheels/

set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Usage: $0 <wheel_file>"
    exit 1
fi

WHEEL="$1"
WORKDIR=$(mktemp -d)
trap "rm -rf $WORKDIR" EXIT

echo "Unpacking wheel: $WHEEL"
unzip -q "$WHEEL" -d "$WORKDIR"

# Find all dylibs with absolute Fortran runtime references
for dylib in "$WORKDIR"/lib/SHERPA-MC/*.dylib; do
    [ -f "$dylib" ] || continue

    # Find absolute paths to Fortran runtime libraries
    deps=$(otool -L "$dylib" | grep -oE '/[^ ]+lib(gfortran|quadmath|gcc_s)[^ ]*\.dylib' || true)
    if [ -z "$deps" ]; then
        continue
    fi

    echo "Fixing: $(basename "$dylib")"
    for dep in $deps; do
        basename_dep=$(basename "$dep")
        echo "  $dep -> @rpath/$basename_dep"
        install_name_tool -change "$dep" "@rpath/$basename_dep" "$dylib"
    done

    # Add common Homebrew gfortran RPATH entries if not already present
    existing_rpaths=$(otool -l "$dylib" | grep -A1 LC_RPATH | grep path | awk '{print $2}' || true)
    for rpath in /opt/homebrew/lib/gcc/current /usr/local/lib/gcc/current; do
        if ! echo "$existing_rpaths" | grep -qF "$rpath"; then
            echo "  Adding RPATH: $rpath"
            install_name_tool -add_rpath "$rpath" "$dylib"
        fi
    done
done

# Repack the wheel
OUTDIR="$(dirname "$WHEEL")/fixed_wheels"
mkdir -p "$OUTDIR"
WHEEL_NAME=$(basename "$WHEEL")
echo "Repacking wheel: $OUTDIR/$WHEEL_NAME"
(cd "$WORKDIR" && zip -q -r "$OUTDIR/$WHEEL_NAME" .)

echo "Done. Fixed wheel: $OUTDIR/$WHEEL_NAME"
