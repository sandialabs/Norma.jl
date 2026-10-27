#!/bin/bash
# Regenerate the meshes of this example with Cubit at a chosen refinement
# level and for a chosen number of Schwarz subdomains:
#
#   --refine R   element sizes are the base sizes divided by R (default 1;
#                2 halves the element size, 0.5 doubles it)
#   --sd N       number of subdomains, 2, 3 or 4 (default 2):
#                  2   inner, outer
#                  3   core, filler_in, outer
#                  4   core, filler_in, filler_out, shell
#   --quarter    keep only the quarter x >= 0, y >= 0, with the symmetry
#                sets nsx-/ssx- on x = 0 and nsy-/ssy- on y = 0
#
# Meshes go to nonoverlap/ and overlap/, plus monolithic/cylinders.g in
# every case, all under dynamic/ for the full model and under
# quasistatic/quarter/ with --quarter. Either value option may
# be written --refine 2 or --refine=2. Run from any directory; Cubit must
# be on PATH or in CUBIT.
#
#   ./generate-meshes.sh [--refine R] [--sd N] [--quarter]
#
# The committed meshes of the dynamic example are the full meshes at
# refinement 1 from --sd 2, and those of the quasi-static example the
# quarter meshes at refinement 0.5: monolithic/cylinders.g and
# nonoverlap/inner.g and outer.g from --sd 2, nonoverlap/core.g, filler_in.g,
# filler_out.g and shell.g from --sd 4:
#
#   ./generate-meshes.sh --sd 2
#   ./generate-meshes.sh --refine 0.5 --sd 2 --quarter
#   ./generate-meshes.sh --refine 0.5 --sd 4 --quarter
set -euo pipefail

usage() {
  echo "usage: $0 [--refine R] [--sd 2|3|4] [--quarter]" >&2
  exit "${1:-1}"
}

refinement=1
subdomains=2
quarter=0
while (( $# > 0 )); do
  case $1 in
    --quarter)   quarter=1; shift ;;
    --refine)    (( $# > 1 )) || usage; refinement=$2; shift 2 ;;
    --refine=*)  refinement=${1#*=}; shift ;;
    --sd)        (( $# > 1 )) || usage; subdomains=$2; shift 2 ;;
    --sd=*)      subdomains=${1#*=}; shift ;;
    -h|--help)   usage 0 ;;
    *)           echo "unknown argument '$1'" >&2; usage ;;
  esac
done

[[ $refinement =~ ^([0-9]+\.?[0-9]*|\.[0-9]+)$ ]] && [[ ! $refinement =~ ^0*\.?0*$ ]] \
  || { echo "--refine must be a positive number, not '$refinement'" >&2; exit 1; }
case $subdomains in
  2) domains="inner outer" ;;
  3) domains="core filler_in outer" ;;
  4) domains="core filler_in filler_out shell" ;;
  *) echo "--sd must be 2, 3 or 4, not '$subdomains'" >&2; exit 1 ;;
esac

cubit=${CUBIT:-cubit}
command -v "$cubit" > /dev/null || { echo "Cubit not found; set CUBIT to the executable" >&2; exit 1; }
cd "$(dirname "$0")"
prefix="dynamic/"
(( quarter )) && prefix="quasistatic/quarter/"
mkdir -p "${prefix}monolithic" "${prefix}nonoverlap" "${prefix}overlap"
wrapper=$(mktemp ./refinement-XXXX.jou)
trap 'rm -f "$wrapper"' EXIT
printf '${refinement = %s}\n${subdomains = %s}\n${quarter = %s}\nplayback "concentric-cylinders.jou"\n' \
  "$refinement" "$subdomains" "$quarter" > "$wrapper"
"$cubit" -batch -nographics -nojournal -noecho "$wrapper"
meshes=("${prefix}monolithic/cylinders.g")
for dir in nonoverlap overlap; do
  for domain in $domains; do
    meshes+=("$prefix$dir/$domain.g")
  done
done
ls -la "${meshes[@]}"
