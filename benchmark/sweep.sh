#!/bin/bash
# Thread sweep of the benchmark cases, one Julia process per (case, thread
# count), appended to results.txt as BENCH lines. Run with the machine
# otherwise idle: a busy core delays every static chunk of the loops.
#   ./sweep.sh                      # all cases, 1 2 4 8 12 threads
#   ./sweep.sh "1 12" explicit48    # chosen thread counts and cases
set -e
cd "$(dirname "$0")"
threads=${1:-"1 2 4 8 12"}
# The threads pinned to cores, as bin/norma sets it for threaded runs.
export JULIA_EXCLUSIVE=${JULIA_EXCLUSIVE:-1}
shift || true
cases=${@:-"explicit48 implicit24-cg implicit24-cg-jacobi implicit24-direct quasistatic24-cg quasistatic24-jacobi quasistatic24-direct"}
echo "# $(date -u +%Y-%m-%dT%H:%M:%SZ) $(git -C .. rev-parse --short HEAD) $(hostname)" >> results.txt
for case in $cases; do
    for n in $threads; do
        julia -t "$n" --project=.. run.jl "$case.yaml" 2>&1 | grep "^BENCH" | tee -a results.txt
    done
done
