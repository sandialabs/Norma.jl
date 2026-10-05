#!/usr/bin/env python3
"""Energy histories of the dynamic stability study, one figure per factor.

    python3 plot.py [--runs DIR] [--output DIR] [--log] [--tier A ...]

Each figure shows the effect of one factor (one tier of matrix.jl) with the
other factors at their defaults: one panel per integrator pair, one curve
per level of the factor, the energy E(t)/E(0) against time, with the
monolithic reference of the same level dashed in gray. A run that failed
ends where its history ends, marked with a cross. Figures go to the output
directory (default: figures) as PNG files named after the factor. --log
draws the energy ratio on a logarithmic axis, which shows exponential
growth as a straight line.

The membership of a case in a tier is read from its case.txt, so cases
generated for one tier appear in every figure whose other factors they
match (the default case appears in all of them).

Needs matplotlib.
"""
import argparse
import csv
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

DEFAULT = {"level": "1", "ratio": "1.0", "tolerance": "1.0e-12", "steps": "1", "solver": "anderson"}

# factor -> (tier, file name, title, pairs shown, the levels in plotting
# order as (predicate on the case, label))
FACTORS = {
    "coupling": ("A", "coupling", "Coupling", ["II", "IE", "EI", "EE"], [
        ("no-dn", "Dirichlet-Neumann exchange (baseline)"),
        ("no-cd", "constrained exchange"),
    ]),
    "ratio": ("B", "mesh-ratio", "Mesh ratio (clamped part coarsened)", ["II", "EE"], [
        ("1.0", "conforming"), ("0.75", "1:0.75"), ("0.5", "1:0.5"),
    ]),
    "tolerance": ("C", "tolerance", "Schwarz tolerance", ["II", "EE"], [
        ("1.0e-6", "1e-6"), ("1.0e-8", "1e-8"), ("1.0e-10", "1e-10"), ("1.0e-12", "1e-12"),
    ]),
    "steps": ("D", "time-steps", "Time steps", ["II", "EE", "IE"], [
        ("1", "equal steps"),
        ("4coarse", "4:1, Dirichlet side: clamped part (coarse step)"),
        ("4fine", "4:1, Dirichlet side: free part (fine step)"),
    ]),
    "solver": ("E", "solver", "Solver of the interface problem", ["II", "EE"], [
        ("anderson", "Anderson acceleration"), ("aitken", "Aitken, recursive"), ("fixed", "fixed relaxation 0.5"),
        ("direct", "direct interface solve"),
    ]),
    "level": ("F", "refinement", "Refinement level", ["II", "EE"], [
        ("1", "level 1"), ("2", "level 2"), ("4", "level 4"),
    ]),
}
COLORS = ["#0075A9", "#D97A00", "#1F8A4C", "#8A2B78", "#C41D24"]


def read_case(directory):
    info = {}
    with open(os.path.join(directory, "case.txt")) as f:
        for line in f:
            if ":" not in line:
                continue
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
        rows = [r for r in rows if r["total_energy"] not in ("NaN", "")]
        if rows:
            e0 = float(rows[0]["total_energy"])
            times = [1.0e3 * float(r["time"]) for r in rows]
            ratios = [float(r["total_energy"]) / e0 for r in rows]
    info["times"], info["ratios"] = times, ratios
    return info


def level_of(case, factor):
    """The value of the factor for a case, in the keys of FACTORS."""
    if factor == "coupling":
        return case["coupling"]
    if factor == "steps":
        return "1" if case["steps"] == "1" else case["steps"] + case["role"]
    return case[factor]


def at_defaults(case, factor):
    """True when every factor other than `factor` is at its default."""
    if case["coupling"] != "no-cd" and factor != "coupling":
        return False
    for key, value in DEFAULT.items():
        if key == factor or (factor == "steps" and key == "steps"):
            continue
        if key == "tolerance":
            if float(case[key]) != float(value):
                return False
        elif case[key] != value:
            return False
    return True


def draw(factor, cases, references, output, log):
    tier, stem, title, pairs, levels = FACTORS[factor]
    members = [c for c in cases if c["coupling"] != "mono" and c["pair"] in pairs and at_defaults(c, factor)]
    if not any(len(c["times"]) > 1 for c in members):
        return None
    figure, axes = plt.subplots(1, len(pairs), figsize=(4.6 * len(pairs), 4.2), sharey=True, squeeze=False)
    for axis, pair in zip(axes.flat, pairs):
        kind = pair[0] if pair in ("II", "EE") else "I"
        for reference in references:
            if reference["pair"] != kind:
                continue
            if factor != "level" and reference["level"] != "1":
                continue
            axis.plot(reference["times"], reference["ratios"], color="0.45", lw=1.0, ls="--",
                      label="monolithic" + (f", level {reference['level']}" if factor == "level" else ""))
        for color, (key, label) in zip(COLORS, levels):
            match = [c for c in members if c["pair"] == pair and level_of(c, factor) == key]
            if not match or not match[0]["times"]:
                continue
            c = match[0]
            axis.plot(c["times"], c["ratios"], color=color, lw=1.3, label=label)
            if c["status"] == "failed":
                axis.plot(c["times"][-1], c["ratios"][-1], "x", color=color, ms=8)
        axis.set_title(pair, fontsize=11)
        axis.axhline(1.0, color="k", lw=0.5)
        axis.grid(alpha=0.3)
        axis.set_xlabel("time (ms)")
        if log:
            axis.set_yscale("log")
    drawn = [c for c in members + references if len(c["times"]) > 1]
    t_max = max(c["times"][-1] for c in drawn)
    low = min(min(c["ratios"]) for c in drawn)
    high = max(max(c["ratios"]) for c in drawn)
    if log:
        low, high = max(min(low, 0.5), 1.0e-2), min(max(high, 2.0), 1.0e3)
    else:
        # Runs that end by element inversion reach ratios in the hundreds;
        # the frame keeps the range in which the others are legible.
        low, high = max(low, 0.0), min(high, 2.0)
        margin = max(0.05 * (high - low), 1.0e-3)
        low, high = low - margin, high + margin
    axes.flat[0].set_xlim(0.0, t_max)
    axes.flat[0].set_ylim(low, high)
    axes.flat[0].set_ylabel("E(t) / E(0)")
    handles, labels = [], []
    for axis in axes.flat:
        for handle, text in zip(*axis.get_legend_handles_labels()):
            if text not in labels:
                handles.append(handle)
                labels.append(text)
    figure.legend(handles, labels, loc="lower center", ncol=min(len(labels), 4), fontsize=9, frameon=False)
    figure.suptitle(f"{title} (other factors at their defaults)", fontsize=12)
    figure.tight_layout(rect=(0, 0.12 if len(labels) > 4 else 0.08, 1, 0.95))
    path = os.path.join(output, f"{stem}{'_log' if log else ''}.png")
    figure.savefig(path, dpi=130)
    plt.close(figure)
    return path


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    parser = argparse.ArgumentParser()
    parser.add_argument("--runs", default=os.path.join(here, "runs"))
    parser.add_argument("--output", default=os.path.join(here, "figures"))
    parser.add_argument("--log", action="store_true")
    parser.add_argument("--tier", nargs="*", default=None, help="tiers to plot (default: all with data)")
    args = parser.parse_args()
    cases = []
    for name in sorted(os.listdir(args.runs)):
        directory = os.path.join(args.runs, name)
        if os.path.exists(os.path.join(directory, "case.txt")):
            cases.append(read_case(directory))
    os.makedirs(args.output, exist_ok=True)
    references = [c for c in cases if c["coupling"] == "mono"]
    for factor, (tier, *_) in FACTORS.items():
        if args.tier and tier not in args.tier:
            continue
        path = draw(factor, cases, references, args.output, args.log)
        if path:
            print(path)


if __name__ == "__main__":
    main()
