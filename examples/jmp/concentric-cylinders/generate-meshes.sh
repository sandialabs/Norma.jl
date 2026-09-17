#!/bin/bash
# Regenerate the five meshes of this example with Cubit at a chosen
# refinement level (default 1, the committed meshes; 2 halves the element
# size). Run from any directory; Cubit must be on PATH or in CUBIT.
#
#   ./generate-meshes.sh [refinement]
set -euo pipefail
refinement=${1:-1}
cubit=${CUBIT:-cubit}
command -v "$cubit" > /dev/null || { echo "Cubit not found; set CUBIT to the executable" >&2; exit 1; }
cd "$(dirname "$0")"
wrapper=$(mktemp ./refinement-XXXX.jou)
trap 'rm -f "$wrapper"' EXIT
printf '${refinement = %s}\nplayback "concentric-cylinders.jou"\n' "$refinement" > "$wrapper"
"$cubit" -batch -nographics -nojournal -noecho "$wrapper"
ls -la monolithic/*.g nonoverlap/*.g overlap/*.g
