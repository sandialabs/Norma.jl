#!/usr/bin/env python3
"""Energy histories of the dynamic stability study.

    python3 plot.py [--runs DIR] [--output DIR] [--log]

For every problem, refinement level, and (for the beam) mesh ratio found
under the runs directory, draws one figure with a panel per coupling: the
blended total energy E(t)/E(0) of each integrator pair, with the monolithic
reference of the same level in gray. A run that failed ends where its
history ends, marked with a cross. Figures go to the output directory
(default: figures) as PNG files. --log draws the energy ratio on a
logarithmic axis, which shows exponential growth as a straight line.

Needs matplotlib.
"""
import argparse
import csv
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

COUPLINGS = [
    ("ov-dir", "overlap, Dirichlet"),
    ("ov-imp", "overlap, impedance"),
    ("no-dn", "nonoverlap, Dirichlet-Neumann"),
    ("no-imp", "nonoverlap, paired impedance"),
]
PAIRS = ["II", "IE", "EI", "EE"]
COLORS = {"II": "#1f77b4", "IE": "#ff7f0e", "EI": "#2ca02c", "EE": "#d62728"}


def read_case(directory):
    info = {}
    with open(os.path.join(directory, "case.txt")) as f:
        for line in f:
            key, value = line.split(":", 1)
            info[key.strip()] = value.strip()
    status = "not run"
    status_file = os.path.join(directory, "status.txt")
    if os.path.exists(status_file):
        with open(status_file) as f:
            status = f.readline().split(":", 1)[1].strip()
    info["status"] = status
    times, ratios = [], []
    energy_file = os.path.join(directory, "run-energy.csv")
    if os.path.exists(energy_file):
        with open(energy_file) as f:
            rows = list(csv.DictReader(f))
        if rows:
            e0 = float(rows[0]["total_energy"])
            times = [1.0e3 * float(r["time"]) for r in rows]
            ratios = [float(r["total_energy"]) / e0 for r in rows]
    info["times"], info["ratios"] = times, ratios
    return info


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    parser = argparse.ArgumentParser()
    parser.add_argument("--runs", default=os.path.join(here, "runs"))
    parser.add_argument("--output", default=os.path.join(here, "figures"))
    parser.add_argument("--log", action="store_true")
    args = parser.parse_args()
    cases = []
    for name in sorted(os.listdir(args.runs)):
        directory = os.path.join(args.runs, name)
        if os.path.exists(os.path.join(directory, "case.txt")):
            cases.append(read_case(directory))
    os.makedirs(args.output, exist_ok=True)
    groups = sorted({(c["problem"], c["level"], c["ratio"]) for c in cases if c["coupling"] != "mono"})
    for problem, level, ratio in groups:
        references = [c for c in cases if c["problem"] == problem and c["level"] == level and c["coupling"] == "mono"]
        members = [c for c in cases if (c["problem"], c["level"], c["ratio"]) == (problem, level, ratio)]
        if not any(len(c["times"]) > 1 for c in members + references):
            continue
        figure, axes = plt.subplots(2, 2, figsize=(11, 7), sharex=True, sharey=True)
        for axis, (coupling, title) in zip(axes.flat, COUPLINGS):
            for reference in references:
                axis.plot(reference["times"], reference["ratios"], color="0.5", lw=1.0,
                          ls="--" if reference["pair"] == "I" else ":", label=f"monolithic {reference['pair']}")
            for pair in PAIRS:
                match = [c for c in cases if (c["problem"], c["level"], c["ratio"], c["coupling"], c["pair"])
                         == (problem, level, ratio, coupling, pair)]
                if not match or not match[0]["times"]:
                    continue
                c = match[0]
                axis.plot(c["times"], c["ratios"], color=COLORS[pair], lw=1.2, label=pair)
                if c["status"] == "failed":
                    axis.plot(c["times"][-1], c["ratios"][-1], "x", color=COLORS[pair], ms=8)
            axis.set_title(title, fontsize=10)
            axis.axhline(1.0, color="k", lw=0.5)
            axis.grid(alpha=0.3)
            if args.log:
                axis.set_yscale("log")
        # Limits from the data, with a minimum vertical span: an energy ratio
        # that barely moves otherwise defeats the automatic tick placement.
        drawn = [c for c in members + references if len(c["times"]) > 1]
        t_max = max(c["times"][-1] for c in drawn)
        low = min(min(c["ratios"]) for c in drawn)
        high = max(max(c["ratios"]) for c in drawn)
        if args.log:
            low, high = min(low, 0.5), max(high, 2.0)
        else:
            margin = max(0.05 * (high - low), 1.0e-3)
            low, high = low - margin, high + margin
        axes.flat[0].set_xlim(0.0, t_max)
        axes.flat[0].set_ylim(low, high)
        for axis in axes[1]:
            axis.set_xlabel("time (ms)")
        for axis in axes[:, 0]:
            axis.set_ylabel("E(t) / E(0)")
        axes.flat[0].legend(fontsize=8)
        label = ("beam" if problem == "beam" else "nested cylinders") + f", level {level}" + (f", mesh ratio {ratio}" if problem == "beam" else "")
        figure.suptitle(label)
        figure.tight_layout()
        suffix = f"_r{int(round(100 * float(ratio))):03d}" if problem == "beam" else ""
        path = os.path.join(args.output, f"{problem}_L{level}{suffix}{'_log' if args.log else ''}.png")
        figure.savefig(path, dpi=130)
        plt.close(figure)
        print(path)


if __name__ == "__main__":
    main()
