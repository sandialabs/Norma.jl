"""
Plot the computed solution (z-displacement, z-velocity, z-acceleration) for
each of the two overlapping subdomains in the
dynamic-linear-elastic-single-gaussian-opinf-rom-fom-multi-swap example,
snapshot by snapshot.

This example uses Norma's FOM<->ROM swapping capability on BOTH
subdomains, each swapping twice over the course of the run:

    Subdomain 1 (coarse, z in [-0.5, -0.2]): ROM -> FOM -> ROM
    Subdomain 2 (fine,   z in [-0.3,  0.5]): FOM -> ROM -> FOM

Every time a subsim swaps onto a model whose Exodus output file name would
collide with the one already on disk, Norma renames the new output (and
therefore every CSV prefix associated with it) by appending "-phaseN". For
this run that produces the following on-disk CSV prefixes (confirmed by
inspecting the actual files and clamped-swap.yaml / clamped-swap.log):

    Subdomain 1: clamped-rom-1   (snapshots 0000-0130, ROM)
                 clamped-fom-1-phase3 (snapshots 0140-0880, FOM)
                 clamped-rom-1-phase4 (snapshots 0890-1000, ROM)

    Subdomain 2: clamped-fom-2   (snapshots 0000-0350, FOM)
                 clamped-rom-2-phase2 (snapshots 0360-0970, ROM)
                 clamped-fom-2-phase3 (snapshots 0980-1000, FOM)

A subdomain's CSV files follow the naming convention:
    <prefix>-disp-XXXX.csv
    <prefix>-velo-XXXX.csv
    <prefix>-acce-XXXX.csv
    <prefix>-time-XXXX.csv
where XXXX is a zero-padded snapshot index (step number). This script
walks across the swap boundaries transparently -- for any given snapshot
index it looks up which prefix segment is "live" for each subdomain, reads
that segment's CSV files, and plots the result as one continuous curve per
subdomain. The legend is rebuilt every frame so it always reflects which
model (ROM or FOM) is actually active for each subdomain at that instant,
labeled ROM1 / FOM1 / ROM2 / FOM2 as appropriate.

To re-use this script for a different (non-swapping, or differently-swapping)
run, edit SUBDOMAIN_PREFIXES below: each entry can be a single string prefix
(for a subdomain that never swaps) or a list of (prefix, last_index) tuples
for a subdomain that swaps one or more times, with last_index being the last
snapshot index (inclusive) that prefix's files are used for.

By default the figure is shown live and updated frame by frame as it plots
(matching the original MATLAB-derived Plot_Soln behavior), pausing
--pause-time seconds (default 0.1) between frames, and a new snapshot is
plotted for every available CSV file (every 10 simulation steps). Pass
--save-figs to also write each plotted frame out as a PNG to --output-dir --
since every CSV file is plotted by default, this also means every CSV file
gets its own saved PNG.

Before the live loop starts, every snapshot that will be plotted is scanned
once to find the global min/max of each field (disp/velo/acce) across both
subdomains and the whole run; each subplot's y-axis limits are then fixed
to that range (with a small margin) for the entire run. The x-axis (z) is
likewise fixed for every subplot to --x-min/--x-max (default -0.5/0.5, the
full bar extent), so the axes stay put and don't rescale from frame to
frame.

The exact analytical solution (a clamped-clamped Gaussian-pulse string
problem -- see the EXACT_* derivation comments above) is plotted as a
thin black solid curve alongside each subdomain's simulated curve by
default. ROM segments are drawn with a thick dashed line (no marker);
FOM segments with a thinner solid line and circle markers -- the extra
line weight on ROM keeps it visually distinct from the thin black exact
curve even where the two nearly overlap. The relative L2 error against
the exact solution is not shown in the legend (to keep it uncluttered),
but an RMS-over-time summary (per subdomain, per field) is printed to
the console before the plotting loop starts.
Pass --no-exact to turn the exact-solution overlay and error reporting
off entirely (e.g. for a different problem where this exact solution
doesn't apply). Useful flags:

    python plot_solns_overlap_2sd_rom_fom_multi_swap.py
    python plot_solns_overlap_2sd_rom_fom_multi_swap.py --save-figs
    python plot_solns_overlap_2sd_rom_fom_multi_swap.py --save-figs --time-incr 10
    python plot_solns_overlap_2sd_rom_fom_multi_swap.py --pause-time 0.5
    python plot_solns_overlap_2sd_rom_fom_multi_swap.py --x-min -0.6 --x-max 0.6
    python plot_solns_overlap_2sd_rom_fom_multi_swap.py --no-exact
    python plot_solns_overlap_2sd_rom_fom_multi_swap.py --data-dir /path/to/run --time-incr 5
"""

from __future__ import annotations

import argparse
import os
import numpy as np
import matplotlib.pyplot as plt


# ----------------------------------------------------------------------
# Exact solution.
#
# The IC for this problem (see clamped-rom-1.yaml / clamped-fom-2.yaml) is
#     u(z,0)  =  a * exp(-z^2 / (2 s^2))
#     u_t(z,0) = -a*c*z/s^2 * exp(-z^2 / (2 s^2))  =  c * u'(z,0)
# which is exactly the initial condition for a PURE LEFT-MOVING Gaussian
# pulse in free space: u(z,t) = f(z+ct) with f(z) = a*exp(-z^2/(2 s^2))
# (verified numerically against the early-time, pre-reflection simulation
# data to ~0.01-0.1% relative error).
#
# The bar is clamped (u=0) at both z=-L and z=+L (L=0.5), so the pulse
# reflects off the wall it travels toward (z=-L, since it's left-moving)
# within the simulated time window [0, 1e-3] (it never reaches +L). The
# reflected solution is obtained via the standard method of images for a
# clamped-clamped string: extend the IC data to be ODD and PERIODIC with
# period 4L about both walls, then the free-space d'Alembert solution
# built from those extensions automatically satisfies u(-L,t)=u(L,t)=0
# for all t. Closed-form derivatives give velocity/acceleration exactly
# (no finite differencing needed). This was validated against a
# fine-mesh finite-difference reference solution of the actual PDE+BCs
# (relative error <0.02% away from points where the true solution is
# itself ~0, e.g. exactly at the moment of reflection, t=T/2) and against
# the actual ROM/FOM simulation CSVs across the whole run (residual
# 1-6%, consistent with ROM/FOM/time-discretization error, not a flaw in
# this exact solution).
# ----------------------------------------------------------------------
EXACT_A = 0.001       # IC amplitude
EXACT_S = 0.02        # IC width
EXACT_C = np.sqrt(1.0e9 / 1.0e3)  # wave speed = sqrt(E/rho)
EXACT_L = 0.5         # clamped-end position (walls at z=-L and z=+L)


def _exact_f(x):
    return EXACT_A * np.exp(-x**2 / (2.0 * EXACT_S**2))


def _exact_fprime(x):
    return EXACT_A * (-x / EXACT_S**2) * np.exp(-x**2 / (2.0 * EXACT_S**2))


def _exact_fpp(x):
    return (EXACT_A * (x**2 - EXACT_S**2) / EXACT_S**4
            * np.exp(-x**2 / (2.0 * EXACT_S**2)))


def _odd_periodic(x, direct_func, reflected_func):
    """Evaluate a function that is odd and 4L-periodic about both z=-L
    and z=+L, given its formula on the 'direct' piece (the fundamental
    domain [-L, L]) and its formula on the 'reflected' piece ((L, 3L)),
    which by construction must already be -direct_func(2L - xm)-shaped
    (or the appropriate derivative thereof) for the result to actually be
    smooth/continuous -- see the EXACT_* derivation note above."""
    period = 4.0 * EXACT_L
    xm = np.mod(x + EXACT_L, period) - EXACT_L  # fundamental range [-L, 3L)
    return np.where(xm <= EXACT_L, direct_func(xm), reflected_func(xm))


def _Phi(x):
    """Odd-periodic extension of the position IC f(z)."""
    return _odd_periodic(x, _exact_f, lambda xm: -_exact_f(2 * EXACT_L - xm))


def _Phi_prime(x):
    return _odd_periodic(x, _exact_fprime, lambda xm: _exact_fprime(2 * EXACT_L - xm))


def _Phi_pp(x):
    return _odd_periodic(x, _exact_fpp, lambda xm: -_exact_fpp(2 * EXACT_L - xm))


def _H(x):
    """Continuous antiderivative of Psi (the odd-periodic extension of the
    velocity IC's profile), used to evaluate the d'Alembert integral term
    in closed form."""
    return _odd_periodic(x, lambda xm: EXACT_C * _exact_f(xm),
                          lambda xm: EXACT_C * _exact_f(2 * EXACT_L - xm))


def _Psi(x):
    return _odd_periodic(x, lambda xm: EXACT_C * _exact_fprime(xm),
                          lambda xm: -EXACT_C * _exact_fprime(2 * EXACT_L - xm))


def _Psi_prime(x):
    return _odd_periodic(x, lambda xm: EXACT_C * _exact_fpp(xm),
                          lambda xm: EXACT_C * _exact_fpp(2 * EXACT_L - xm))


def exact_disp(z, t):
    """Exact z-displacement at position(s) z and time t."""
    z = np.asarray(z, dtype=float)
    A = z - EXACT_C * t
    B = z + EXACT_C * t
    return 0.5 * (_Phi(A) + _Phi(B)) + (_H(B) - _H(A)) / (2.0 * EXACT_C)


def exact_velo(z, t):
    """Exact z-velocity (closed-form time derivative of exact_disp)."""
    z = np.asarray(z, dtype=float)
    A = z - EXACT_C * t
    B = z + EXACT_C * t
    return (0.5 * EXACT_C * (_Phi_prime(B) - _Phi_prime(A))
            + 0.5 * (_Psi(B) + _Psi(A)))


def exact_acce(z, t):
    """Exact z-acceleration (closed-form second time derivative of
    exact_disp)."""
    z = np.asarray(z, dtype=float)
    A = z - EXACT_C * t
    B = z + EXACT_C * t
    return (0.5 * EXACT_C**2 * (_Phi_pp(B) + _Phi_pp(A))
            + 0.5 * EXACT_C * (_Psi_prime(B) - _Psi_prime(A)))


EXACT_FUNCS = {"disp": exact_disp, "velo": exact_velo, "acce": exact_acce}


# ----------------------------------------------------------------------
# Run configuration
# ----------------------------------------------------------------------
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--data-dir", default=".",
                     help="Directory containing the CSV files (default: '.')")
parser.add_argument("--save-figs", action="store_true",
                     help="Also save each frame as a PNG to --output-dir "
                          "(default: off; frames are only shown live)")
parser.add_argument("--output-dir", default="frames",
                     help="Directory to write snapshot images to, if "
                          "--save-figs is given (default: 'frames')")
parser.add_argument("--time-incr", type=int, default=1,
                     help="Plot every Nth available CSV snapshot, in addition "
                          "to every swap-boundary snapshot. Since CSVs are "
                          "already written every 10 simulation steps, the "
                          "default of 1 plots/saves every available CSV "
                          "file (i.e. every 10 steps).")
parser.add_argument("--pause-time", type=float, default=0.1,
                     help="Seconds to pause between frames while displaying "
                          "live (default: 0.1)")
parser.add_argument("--x-min", type=float, default=-0.5,
                     help="Fixed x-axis (z) lower limit for every subplot "
                          "(default: -0.5, the full bar extent)")
parser.add_argument("--x-max", type=float, default=0.5,
                     help="Fixed x-axis (z) upper limit for every subplot "
                          "(default: 0.5, the full bar extent)")
parser.add_argument("--no-exact", action="store_true",
                     help="Don't plot the exact solution or compute "
                          "relative error against it (default: shown)")
args = parser.parse_args()

data_dir = args.data_dir
base_prefix = "clamped"
nsd = 2                  # number of subdomains
istart = 0
iend = 1000               # last snapshot index in the run
istep = 10                # CSV output stride (matches CSV output interval)

dt = 1.0e-6               # simulation time step (s)

save_figs = args.save_figs
output_dir = args.output_dir
fig_format = "png"
pause_time = args.pause_time

colors = ["b", "r"]       # one color per subdomain (cycles if nsd > 2)
time_incr = args.time_incr  # plot/save every Nth available CSV snapshot
                             # (default 1 -> every CSV file, i.e. every istep
                             # simulation steps)
show_exact = not args.no_exact

# ----------------------------------------------------------------------
# Per-subdomain CSV filename-prefix segments.
#
# Most subdomains use a single fixed prefix for the whole run, e.g.
#     SUBDOMAIN_PREFIXES[k] = f'{base_prefix}-rom-{k}'
#
# A subdomain that swaps prefixes partway through the run (Norma's
# FOM<->ROM swapping capability) instead gets a list of
# (prefix, last_index) segments, in order, where last_index is the last
# snapshot index (inclusive) that segment's prefix is used for.
#
# This run's actual on-disk segments (confirmed against the CSV files
# themselves and clamped-swap.log):
#   Subdomain 1: ROM (0-100) -> FOM (110-880) -> ROM (890-1000)
#   Subdomain 2: FOM (0-350) -> ROM (360-970) -> FOM (980-1000)
#
# NOTE: clamped-swap.yaml currently specifies t_swap=0.0008/0.00085 for
# the second pair of swaps (which would move these boundaries to
# 800/810 and 850/860), but the CSV files on disk are still from a run
# made BEFORE that yaml edit, so the boundaries below (880/890, 970/980)
# match what's actually here. If you re-run the simulation with the
# current clamped-swap.yaml, re-derive these from the new files (see
# validate_segments() below, which will tell you immediately if this
# ever goes stale again).
#
# To plot a plain (non-swapping) case, just set every entry to a single
# string prefix, e.g. SUBDOMAIN_PREFIXES = {1: 'clamped-1', 2: 'clamped-2'}
# ----------------------------------------------------------------------
SUBDOMAIN_PREFIXES = {
    1: [(f"{base_prefix}-rom-1", 100),
        (f"{base_prefix}-fom-1-phase2", 850),   # was 880
        (f"{base_prefix}-rom-1-phase3", iend)],
    2: [(f"{base_prefix}-fom-2", 350),
        (f"{base_prefix}-rom-2-phase2", 500),   # was 970
        (f"{base_prefix}-fom-2-phase3", iend)],
}


def normalize_prefix_segments(spec):
    """Turn a SUBDOMAIN_PREFIXES entry (either a plain string prefix, or a
    list of (prefix, last_index) segments) into a normalized list of
    (prefix, last_index) segments, sorted by last_index."""
    if isinstance(spec, str):
        return [(spec, iend)]
    return sorted(spec, key=lambda seg: seg[1])


def prefix_for_index(segments, i):
    """Given a sorted list of (prefix, last_index) segments, return the
    prefix that should be used for snapshot index i."""
    for prefix, last_index in segments:
        if i <= last_index:
            return prefix
    # fall back to the last segment if i exceeds every last_index
    return segments[-1][0]


def is_rom_prefix(prefix):
    """True if this CSV filename prefix corresponds to a ROM segment,
    False for FOM (or anything else)."""
    return "rom" in prefix.lower()


def model_type_label(prefix, k):
    """Classify a CSV filename prefix as ROM or FOM and return the legend
    label for subdomain k, e.g. 'ROM1' or 'FOM2'."""
    if is_rom_prefix(prefix):
        return f"ROM{k}"
    if "fom" in prefix.lower():
        return f"FOM{k}"
    return f"SD{k}"


subdomain_segments = {k: normalize_prefix_segments(SUBDOMAIN_PREFIXES[k])
                       for k in range(1, nsd + 1)}


def validate_segments():
    """Check that every configured (prefix, last_index) segment actually
    has at least one matching CSV file on disk, at the first index of
    that segment. Catches stale SUBDOMAIN_PREFIXES entries (e.g. after
    a re-run shifts which '-phaseN' suffix gets applied) up front with a
    clear error, instead of silently dropping a subdomain's curve partway
    through the plot."""
    import glob

    problems = []
    for k in range(1, nsd + 1):
        segments = subdomain_segments[k]
        for seg_idx, (prefix, last_index) in enumerate(segments):
            start_index = istart if seg_idx == 0 else segments[seg_idx - 1][1] + istep
            probe = os.path.join(data_dir, f"{prefix}-disp-{start_index:04d}.csv")
            if os.path.isfile(probe):
                continue
            # Look for what prefixes actually exist with this start index,
            # to suggest the likely correct name.
            pattern = os.path.join(data_dir, f"*-disp-{start_index:04d}.csv")
            candidates = sorted(
                os.path.basename(p)[:-len(f"-disp-{start_index:04d}.csv")]
                for p in glob.glob(pattern)
            )
            problems.append(
                f"  subdomain {k}, segment {seg_idx} (prefix '{prefix}', "
                f"expected to start at index {start_index}): {probe} not found.\n"
                f"    Prefixes that DO have a snapshot at index {start_index}: "
                f"{candidates if candidates else '(none found)'}"
            )
    if problems:
        raise FileNotFoundError(
            "SUBDOMAIN_PREFIXES does not match the files on disk:\n"
            + "\n".join(problems)
            + "\n\nCheck clamped-swap.log for the actual '-phaseN' renames "
              "(grep -i swap clamped-swap.log) and update SUBDOMAIN_PREFIXES "
              "to match."
        )


validate_segments()


# ----------------------------------------------------------------------
# CSV readers
# ----------------------------------------------------------------------
def read_field_csv(prefix, field, index):
    """Read a <prefix>-<field>-XXXX.csv file (e.g. disp, velo, acce) and
    return it as an (n_nodes, 3) array of [x, y, z] components."""
    fname = os.path.join(data_dir, f"{prefix}-{field}-{index:04d}.csv")
    return np.loadtxt(fname, delimiter=",")


def read_time_csv(prefix, index):
    fname = os.path.join(data_dir, f"{prefix}-time-{index:04d}.csv")
    return float(np.loadtxt(fname))


def read_curr_csv(prefix, index):
    """Current (deformed) nodal coordinates -- used here only to recover
    the *reference* (z-only) nodal positions for plotting, taken from the
    first available snapshot of each segment (index 0 of that segment)."""
    return read_field_csv(prefix, "curr", index)


# ----------------------------------------------------------------------
# Reference z-positions for each subdomain.
#
# The z-positions of the nodes don't change between snapshots within a
# subdomain (only their displacement does), so grab them once per
# subdomain from that subdomain's very first available "curr" snapshot,
# then subtract off the z-displacement at that same snapshot to recover
# the *undeformed* z-coordinate, which is what we plot against.
# ----------------------------------------------------------------------
def reference_z_positions(k):
    segments = subdomain_segments[k]
    first_prefix, _ = segments[0]
    curr0 = read_curr_csv(first_prefix, istart)
    disp0 = read_field_csv(first_prefix, "disp", istart)
    z_ref = curr0[:, 2] - disp0[:, 2]
    return z_ref


print("Loading reference nodal z-positions for each subdomain...")
Z = {k: reference_z_positions(k) for k in range(1, nsd + 1)}
for k in range(1, nsd + 1):
    print(f"  subdomain {k}: {Z[k].size} nodes, "
          f"z in [{Z[k].min():.4f}, {Z[k].max():.4f}]")


# ----------------------------------------------------------------------
# Main plotting loop
# ----------------------------------------------------------------------
if save_figs:
    os.makedirs(output_dir, exist_ok=True)

indices = list(range(istart, iend + 1, istep))
strided = set(indices[::time_incr] if time_incr > 1 else indices)
strided.add(indices[-1])

# Always include the snapshot immediately before and after every swap
# boundary for every subdomain, so the legend transition is never skipped
# over even when --time-incr is large.
for k in range(1, nsd + 1):
    segments = subdomain_segments[k]
    for _, last_index in segments[:-1]:
        if last_index in indices:
            strided.add(last_index)
        next_index = last_index + istep
        if next_index in indices:
            strided.add(next_index)

plot_indices = sorted(strided)

print(f"Plotting {len(plot_indices)} snapshots "
      f"(every {time_incr * istep} steps)...")

field_names = ["disp", "velo", "acce"]
field_titles = ["z-displacement", "z-velocity", "z-acceleration"]
field_ylabels = ["disp_z", "velo_z", "acce_z"]


# ----------------------------------------------------------------------
# Fixed axis limits + data cache.
#
# Scan every snapshot that will be plotted once up front: read each
# field's CSV file (caching it for reuse in the plotting loop below, so
# nothing gets read from disk twice) and track the global min/max of
# each field across both subdomains and the whole run. Each subplot's
# y-axis is then locked to that range (with a small margin) so the axes
# don't rescale frame to frame. The x-axis (z-position) is fixed too,
# from the reference node positions.
# ----------------------------------------------------------------------
print("Scanning all snapshots to compute fixed axis limits...")
field_min = {field: np.inf for field in field_names}
field_max = {field: -np.inf for field in field_names}
data_cache = {}    # (i, field, k) -> simulated z-component values array
exact_cache = {}   # (i, field, k) -> exact-solution z-component values array
rel_error = {field: {k: {} for k in range(1, nsd + 1)} for field in field_names}
exact_norm_cache = {}  # (i, field, k) -> ||exact_values|| (L2 norm over nodes)
diff_norm_cache = {}   # (i, field, k) -> ||values - exact_values|| (L2 norm)

for i in plot_indices:
    t = i * dt
    for field in field_names:
        for k in range(1, nsd + 1):
            prefix = prefix_for_index(subdomain_segments[k], i)
            fname = os.path.join(data_dir, f"{prefix}-{field}-{i:04d}.csv")
            if not os.path.isfile(fname):
                raise FileNotFoundError(
                    f"Expected CSV file not found: {fname}\n"
                    f"This usually means SUBDOMAIN_PREFIXES[{k}] doesn't match "
                    f"the actual on-disk prefix/phase-suffix for snapshot {i}. "
                    f"Check clamped-swap.log for the real '-phaseN' renames "
                    f"and update SUBDOMAIN_PREFIXES accordingly."
                )
            values = read_field_csv(prefix, field, i)[:, 2]
            data_cache[(i, field, k)] = values
            field_min[field] = min(field_min[field], values.min())
            field_max[field] = max(field_max[field], values.max())

            if show_exact:
                exact_values = EXACT_FUNCS[field](Z[k], t)
                exact_cache[(i, field, k)] = exact_values
                field_min[field] = min(field_min[field], exact_values.min())
                field_max[field] = max(field_max[field], exact_values.max())
                exact_norm_cache[(i, field, k)] = np.linalg.norm(exact_values)
                diff_norm_cache[(i, field, k)] = np.linalg.norm(values - exact_values)

y_limits = {}
for field in field_names:
    lo, hi = field_min[field], field_max[field]
    span = hi - lo
    margin = 0.05 * span if span > 0 else max(abs(hi), 1.0) * 0.05
    y_limits[field] = (lo - margin, hi + margin)
    print(f"  {field}: [{lo:.4g}, {hi:.4g}] -> axis limits "
          f"[{y_limits[field][0]:.4g}, {y_limits[field][1]:.4g}]")

if show_exact:
    # Relative error is only meaningful once the exact solution's L2 norm
    # at that snapshot is a non-negligible fraction of the field's overall
    # amplitude scale -- right at a wall reflection the true solution
    # passes through ~0 everywhere, and dividing by that near-zero norm
    # would produce a meaningless, enormous "error" even though the
    # simulated and exact curves visually agree. Snapshots below the
    # threshold are excluded (legend shows no error, RMS summary skips
    # them) rather than shown as a misleading number.
    REL_ERROR_FLOOR_FRACTION = 0.01
    for field in field_names:
        scale = max(abs(field_min[field]), abs(field_max[field]))
        floor = REL_ERROR_FLOOR_FRACTION * scale
        for k in range(1, nsd + 1):
            for i in plot_indices:
                denom = exact_norm_cache[(i, field, k)]
                n_nodes = data_cache[(i, field, k)].size
                if denom > floor * np.sqrt(n_nodes):
                    rel_error[field][k][i] = diff_norm_cache[(i, field, k)] / denom

    print("\nRelative error vs. exact solution (L2 norm over nodes, RMS "
          "over plotted snapshots; snapshots where the exact solution is "
          "negligibly small everywhere -- e.g. exactly at a wall "
          "reflection -- are excluded as not meaningful):")
    for field in field_names:
        for k in range(1, nsd + 1):
            errs = np.array(list(rel_error[field][k].values()))
            n_excluded = len(plot_indices) - errs.size
            rms = np.sqrt(np.mean(errs**2)) if errs.size else np.nan
            note = f" ({n_excluded} snapshot(s) excluded)" if n_excluded else ""
            print(f"  subdomain {k}, {field}: RMS relative error = "
                  f"{rms:.4%}{note}")

z_all = np.concatenate([Z[k] for k in range(1, nsd + 1)])
print(f"  z (data range): [{z_all.min():.4f}, {z_all.max():.4f}]")
x_limits = (args.x_min, args.x_max)  # full bar extent, fixed regardless of per-subdomain data range

plt.ion()
fig, axes = plt.subplots(3, 1, figsize=(8, 10), sharex=False)
fig.show()

for i in plot_indices:
    time = i * dt

    for ax, field, title, ylabel in zip(axes, field_names, field_titles,
                                         field_ylabels):
        ax.cla()

        legend_handles = []
        legend_labels = []

        for k in range(1, nsd + 1):
            prefix = prefix_for_index(subdomain_segments[k], i)
            z = Z[k]
            values = data_cache[(i, field, k)]

            color = colors[(k - 1) % len(colors)]
            label = model_type_label(prefix, k)
            if is_rom_prefix(prefix):
                linestyle, marker, linewidth = "--", None, 3.0
            else:
                linestyle, marker, linewidth = "-", "o", 1.5
            (line,) = ax.plot(z, values, linestyle=linestyle, marker=marker,
                               color=color, markersize=3, linewidth=linewidth,
                               label=label)
            legend_handles.append(line)
            legend_labels.append(label)

            if show_exact:
                exact_values = exact_cache[(i, field, k)]
                (exact_line,) = ax.plot(
                    z, exact_values, linestyle="-", color="k", linewidth=1.0,
                    label="exact" if k == 1 else None
                )
                if k == 1:
                    legend_handles.append(exact_line)
                    legend_labels.append("exact")

        ax.set_title(f"{title}  (snapshot {i}, t = {time:.6e} s)")
        ax.set_xlabel("z")
        ax.set_ylabel(ylabel)
        ax.set_xlim(x_limits)
        ax.set_ylim(y_limits[field])
        ax.legend(legend_handles, legend_labels, loc="upper right")
        ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.canvas.draw()
    plt.pause(max(pause_time, 1e-3))

    if save_figs:
        fig.savefig(os.path.join(output_dir, f"snapshot_{i:04d}.{fig_format}"))

plt.ioff()
print("Done.")
if save_figs:
    print(f"Frames written to: {os.path.abspath(output_dir)}")
plt.show()
