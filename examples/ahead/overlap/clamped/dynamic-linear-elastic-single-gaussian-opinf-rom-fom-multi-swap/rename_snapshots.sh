#!/bin/bash
# Renumbers snapshot_0000.png, snapshot_0010.png, snapshot_0020.png, ...
# to snapshot_0000.png, snapshot_0001.png, snapshot_0002.png, ...
#
# i.e. divides each file's numeric suffix by 10 (assumes every existing
# index is in fact a multiple of 10; any that aren't are skipped with a
# warning rather than silently mis-renamed).
#
# Usage:
#   ./rename_snapshots.sh [directory]
# If no directory is given, the current directory is used.

set -euo pipefail

dir="${1:-.}"

if [ ! -d "$dir" ]; then
    echo "Error: directory '$dir' does not exist." >&2
    exit 1
fi

# Rename into a temporary suffix first, then drop the suffix, so that a
# file's new name never collides with another file's not-yet-renamed old
# name mid-loop (e.g. avoids snapshot_0010.png -> snapshot_0001.png
# clobbering an original snapshot_0001.png if one ever existed).
shopt -s nullglob
files=("$dir"/snapshot_[0-9][0-9][0-9][0-9].png)

if [ ${#files[@]} -eq 0 ]; then
    echo "No snapshot_NNNN.png files found in '$dir'."
    exit 0
fi

echo "Found ${#files[@]} file(s). Renaming..."

tmp_suffix=".renametmp"

for f in "${files[@]}"; do
    base="$(basename "$f")"
    num="${base#snapshot_}"
    num="${num%.png}"

    # Strip leading zeros for arithmetic (avoid octal interpretation), then
    # divide by 10.
    num_noleading="$((10#$num))"

    if (( num_noleading % 10 != 0 )); then
        echo "Skipping '$base': index $num_noleading is not a multiple of 10." >&2
        continue
    fi

    new_index=$(( num_noleading / 10 ))
    new_name=$(printf "snapshot_%04d.png" "$new_index")

    mv -- "$f" "$dir/$new_name$tmp_suffix"
done

for f in "$dir"/snapshot_[0-9][0-9][0-9][0-9].png"$tmp_suffix"; do
    [ -e "$f" ] || continue
    mv -- "$f" "${f%$tmp_suffix}"
    echo "  -> $(basename "${f%$tmp_suffix}")"
done

echo "Done."
