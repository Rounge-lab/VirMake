#!/usr/bin/env bash
#
# Snakemake post-deploy script for DRAMv.yaml.
#
# DRAM 1.5.0 hardcodes two hosts that have since moved:
#   bcb.unl.edu (dbCAN)   -> now 302-redirects to an HTML landing page
#   fileshare.csb.univie.ac.at (VOGDB) -> redirects correctly today, but is
#                                          the same arrangement dbCAN had
# Paths are unchanged on both new hosts, so only the hostname is rewritten.
#
# This must be patched in the package, not worked around with flags:
# prepare_databases() builds its locs dict with
#     locs = {remove_suffix(i, '_loc'): j for i, j in locals().items()
#             if i.endswith('_loc') and j is not None}
# so --dbcan_fam_activities and --vog_annotations are silently ignored, and
# dbcan_subfam_ec is not exposed at all. The unvalidated HTML then reaches the
# description database as "UNIQUE constraint failed: dbcan_description.id".
#
set -euo pipefail

PATCHES="bcb.unl.edu|pro.unl.edu fileshare.csb.univie.ac.at|fileshare.lisc.univie.ac.at"

prefix="${CONDA_PREFIX:?CONDA_PREFIX is not set - run via snakemake, or set it manually}"

pkg_dir=$(find "$prefix/lib" -maxdepth 3 -type d -name mag_annotator | head -n 1)
if [ -z "$pkg_dir" ]; then
    echo "ERROR: mag_annotator package not found under $prefix" >&2
    exit 1
fi

for pair in $PATCHES; do
    old=${pair%%|*}
    new=${pair##*|}
    files=$(grep -rl "$old" "$pkg_dir" --include='*.py' || true)
    if [ -z "$files" ]; then
        echo "post-deploy: no references to $old - already patched?"
        continue
    fi
    for f in $files; do
        cp -n "$f" "$f.orig"
        sed -i "s|$old|$new|g" "$f"
        echo "post-deploy: $f: $old -> $new"
    done
    if grep -rq "$old" "$pkg_dir" --include='*.py'; then
        echo "ERROR: $old still referenced after patching" >&2
        exit 1
    fi
done

echo "post-deploy: DRAM download hosts patched"
