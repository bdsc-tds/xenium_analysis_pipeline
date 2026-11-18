#!/usr/bin/env bash
# Replace the copy of morphology.ome.tif that xeniumranger writes into outs/ with a relative symlink to
# the image of the input Xenium bundle, when both files are byte-identical. Saves one full copy of the
# image per segmentation output. Any other case keeps the copy and only prints why.
#
# Usage: link_morphology.sh <xenium_bundle_dir> <outs_dir>
set -euo pipefail

bundle_img="$1/morphology.ome.tif"
copy="$2/morphology.ome.tif"

# Nothing to do without both files, or when the copy is already a link
if [ ! -f "$bundle_img" ] || [ ! -f "$copy" ] || [ -L "$copy" ]; then
    echo "link_morphology: skipped, missing file or already a link: $copy"
    exit 0
fi

# Keep the copy unless it is byte-identical to the bundle image
if ! cmp -s "$bundle_img" "$copy"; then
    echo "link_morphology: kept the copy, it differs from $bundle_img"
    exit 0
fi

# Relative target between the physical folders, so links stay valid when the output path is itself a symlink
outs_dir=$(realpath "$(dirname "$copy")")
rel=$(realpath --relative-to="$outs_dir" "$(realpath "$bundle_img")")

# Create the link next to the copy and check it resolves to the bundle image before swapping it in
tmp="$outs_dir/.morphology.ome.tif.tmp_link"
rm -f "$tmp"
ln -s "$rel" "$tmp"
if [ ! "$tmp" -ef "$bundle_img" ]; then
    rm -f "$tmp"
    echo "link_morphology: link does not resolve to $bundle_img, kept the copy"
    exit 0
fi
mv -T "$tmp" "$outs_dir/morphology.ome.tif"
echo "link_morphology: $copy -> $rel"
