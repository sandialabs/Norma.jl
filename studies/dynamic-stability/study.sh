#!/bin/bash
# Run the dynamic stability study end to end: generate the cases, run them,
# summarize them, and draw the figures.  Resumable: completed cases are kept
# and skipped, so the same command continues an interrupted study.
#
#   ./study.sh [--jobs N] [--threads T] [--final-time T] [--runs DIR] [TIER ...]
#
# With no tier every tier of matrix.jl runs in order.  Run it detached and
# follow study.log:
#
#   setsid nohup ./study.sh > study.log 2>&1 &
#   tail -f study.log
#
# Defaults: 8 cases at a time with 4 threads each (32 threads).  On rigel
# (336 CPUs, other users' jobs) that leaves the machine mostly free; the
# beam cases do not run faster with more than 4 threads.
set -e
cd "$(dirname "$0")"
JOBS=8; THREADS=4; FINAL=""; RUNS=runs; TIERS=()
while [ $# -gt 0 ]; do
  case "$1" in
    --jobs) JOBS=$2; shift 2 ;;
    --threads) THREADS=$2; shift 2 ;;
    --final-time) FINAL="--final-time $2"; shift 2 ;;
    --runs) RUNS=$2; shift 2 ;;
    *) TIERS+=("$1"); shift ;;
  esac
done
[ ${#TIERS[@]} -eq 0 ] && TIERS=(A B C D E F)
export OPENBLAS_NUM_THREADS=1
echo "study: tiers ${TIERS[*]}, $JOBS jobs x $THREADS threads, runs in $RUNS, Norma $(git -C ../.. log -1 --format=%h), $(date)"
for tier in "${TIERS[@]}"; do
  echo "== tier $tier: generate ($(date +%T))"
  julia --project=../.. generate.jl --runs "$RUNS" $FINAL "$tier"
  echo "== tier $tier: run ($(date +%T))"
  julia --project=../.. run.jl --runs "$RUNS" --jobs "$JOBS" --threads "$THREADS" "$tier"
  echo "== tier $tier: summarize and plot ($(date +%T))"
  julia --project=../.. collect.jl --runs "$RUNS" --output "$RUNS/../summary.csv"
  python3 plot.py --runs "$RUNS" --output "$RUNS/../figures"
  python3 plot.py --runs "$RUNS" --output "$RUNS/../figures" --log
done
echo "study: done ($(date))"
