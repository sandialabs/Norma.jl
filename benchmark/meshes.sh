#!/bin/bash
# Generates the hexahedral cubes of the benchmarks with Cubit: cube24.g
# (13 824 elements) for the implicit cases and cube48.g (110 592 elements) for
# the explicit case. Pass the path of the Cubit executable as the first
# argument when it is not on the PATH.
set -e
cd "$(dirname "$0")"
cubit=${1:-cubit}
for n in 24 48; do
    sed "s/#{n = 24}/#{n = $n}/" cube.jou > cube-$n.jou
    "$cubit" -batch -nographics -nojournal -noecho -information off -input cube-$n.jou > cube-$n.log 2>&1
    rm -f cube-$n.jou
    ls -l cube$n.g
done
