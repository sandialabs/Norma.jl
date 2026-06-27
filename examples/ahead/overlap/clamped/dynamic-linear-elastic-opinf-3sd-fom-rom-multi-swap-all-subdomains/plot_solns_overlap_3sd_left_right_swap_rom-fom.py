"""
Plot the computed solution (z-disp, z-velo, z-acce) for each of the
subdomains together with the exact reference solution, snapshot by
snapshot. Generalized to an arbitrary number of subdomains, with
support for Norma's FOM/ROM swapping capability: a subdomain's CSV
files can switch from one filename prefix to another partway through
the simulation (e.g. "clamped-fom-2-*" up to some snapshot index, then
"clamped-rom-2-*" from then on), while still representing one
continuous subdomain.

This version matches clamped-swap.yaml's round trip on all three
subdomains: subdomain 1 and 3 each go ROM -> FOM (t=0.0003) -> ROM
(t=0.0006), and subdomain 2 goes FOM -> ROM (t=0.00025) -> FOM
(t=0.00055) -- three CSV-filename segments per subdomain. See the
SUBDOMAIN_PREFIXES comment below for the exact snapshot-index
boundaries and the "-phaseN" filename suffixes Norma's
uniquify_swap_output! assigns when a later swap phase reuses an
earlier phase's output mesh file.

At the end, the relative error (over all snapshots) between each
subdomain's computed solution and the exact analytical solution is
printed for all 3 fields (disp/velo/acce) and all subdomains.

Each subdomain's CSV files normally follow the naming convention:
    <prefix>-<k>-refe.csv
    <prefix>-<k>-disp-XXXX.csv
    <prefix>-<k>-velo-XXXX.csv
    <prefix>-<k>-acce-XXXX.csv
    <prefix>-<k>-time-XXXX.csv
where k is the subdomain label and XXXX is a zero-padded snapshot
index.

For a subdomain that swaps prefixes mid-run (FOM -> ROM or ROM -> FOM,
any number of times), specify a list of (prefix, last_index) segments
instead of a single prefix -- see SUBDOMAIN_PREFIXES below. The
reference geometry (*-refe.csv) is read once per subdomain, using the
*first* segment's prefix, since every later phase reuses the original
subdomain's mesh and does not write its own refe.csv.

Run this script from the directory containing the CSV files.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

# ------------------------------------------------------------------
# Settings
# ------------------------------------------------------------------
base_prefix = 'clamped'  # base file name prefix, e.g. 'clamped' -> clamped-1-disp-0000.csv
nsd         = 3          # number of subdomains
istart      = 0           # first snapshot index
istep       = 10          # snapshot index increment
iend        = 1000        # last snapshot index
pad_width   = 4           # zero-padding width of the snapshot index, e.g. 0000

save_figs  = True
pause_time = 0.5
scale = 0.001 * 2

colors = ['b', 'r', 'g', 'm', 'k', 'y']  # cycle if nsd > 6

# ------------------------------------------------------------------
# Per-subdomain file-name prefixes.
#
# Most subdomains use a single fixed prefix for the whole run, e.g.
#     SUBDOMAIN_PREFIXES[k] = f'{base_prefix}-{k}'
#
# A subdomain that swaps prefixes partway through (Norma's swapping
# capability, FOM->ROM or ROM->FOM, any number of times) instead gets
# a list of (prefix, last_index) segments, in order, where last_index
# is the last snapshot index (inclusive) that segment's prefix is used
# for.
#
# clamped-swap.yaml for this run declares a full round trip on every
# subdomain (t_swap = 0.00025, 0.0003, 0.00055, 0.0006), giving THREE
# CSV-filename segments per subdomain:
#
#   - Subdomain 2: FOM -> ROM -> FOM
#         clamped-fom-2          (index 0000 .. 0250)
#         clamped-rom-2          (index 0260 .. 0550)
#         clamped-fom-2-phase2   (index 0560 .. 1000)
#     The first swap's replacement (clamped-rom-2.yaml) writes to a
#     DIFFERENT output mesh file (clamped-3sd-rom-2.e) than the
#     original phase (clamped-3sd-fom-2.e), so no rename happens there.
#     The second swap's replacement (clamped-fom-2.yaml) reuses the
#     very first phase's output file, so THAT one gets renamed to
#     "-phase2" by uniquify_swap_output! (src/swap.jl).
#
#   - Subdomains 1 and 3: ROM -> FOM -> ROM
#         clamped-rom-1          (index 0000 .. 0300)
#         clamped-fom-1-phase2   (index 0310 .. 0600)
#         clamped-rom-1-phase3   (index 0610 .. 1000)
#     Here BOTH swaps reuse the same output mesh file (clamped-3sd-1.e
#     / clamped-3sd-3.e — clamped-rom-1.yaml and clamped-fom-1.yaml
#     always write there), so uniquify_swap_output! renames on every
#     swap: "-phase2" the first time, then "-phase3" the second time
#     (it keeps incrementing past whichever "-phaseN" suffixes already
#     exist on disk).
#
# A second swap landing back on a subsim's ORIGINAL name (e.g. the
# `subsim: clamped-fom-1` plan that fires at t=0.0006) only resolves
# correctly with the fix to apply_swap! in src/swap.jl that registers
# that original, un-renamed name as a surviving alias for the slot —
# without it, that second swap is silently skipped forever and
# subdomains 1 and 3 only ever show two segments, not three. See the
# apply_swap! comments in src/swap.jl for details.
#
# To plot a plain (non-swapping) case, just set every entry to a
# single string prefix, e.g. SUBDOMAIN_PREFIXES = {1: 'clamped-1', ...}
# ------------------------------------------------------------------
SUBDOMAIN_PREFIXES = {
    1: [(f'{base_prefix}-rom-1', 300), (f'{base_prefix}-fom-1-phase2', 600),
        (f'{base_prefix}-rom-1-phase3', iend)],
    2: [(f'{base_prefix}-fom-2', 250), (f'{base_prefix}-rom-2', 550),
        (f'{base_prefix}-fom-2-phase2', iend)],
    3: [(f'{base_prefix}-rom-3', 300), (f'{base_prefix}-fom-3-phase2', 600),
        (f'{base_prefix}-rom-3-phase3', iend)],
}


def normalize_prefix_segments(spec):
    """Turn a SUBDOMAIN_PREFIXES entry (either a plain string prefix,
    or a list of (prefix, last_index) segments) into a normalized list
    of (prefix, last_index) segments, sorted by last_index."""
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


def model_type_label(prefix, k):
    """Derive a short legend label like 'ROM1' or 'FOM3' from a prefix
    such as 'clamped-rom-1' or 'clamped-fom-3-phase2'. Looks for the
    literal substring 'rom' or 'fom' in the prefix (case-insensitive);
    if neither is present, falls back to just the subdomain number so
    the legend still renders something sensible."""
    lower = prefix.lower()
    if 'rom' in lower:
        return f'ROM{k}'
    if 'fom' in lower:
        return f'FOM{k}'
    return str(k)


subdomain_segments = {k: normalize_prefix_segments(SUBDOMAIN_PREFIXES[k])
                       for k in range(1, nsd + 1)}


def uniquetol(arr, tol):
    """Sort and collapse values within `tol` of each other, mirroring
    MATLAB's uniquetol on a sorted array."""
    arr_sorted = np.sort(arr)
    out = [arr_sorted[0]]
    for v in arr_sorted[1:]:
        if abs(v - out[-1]) > tol:
            out.append(v)
    return np.array(out)


def exact_solution(zz, t1):
    """Exact analytical 1-D wave solution (disp, velo, acce) for a
    Gaussian pulse on a bar clamped at both ends, evaluated at
    z-coordinates zz and time t1."""
    c = np.sqrt(1e9 / 1e3)
    a = scale / 2
    b = 0.0
    s = 0.02
    T = 1e-3

    d_ex = 0.5 * a * (np.exp(-(zz - c * t1 - b) ** 2 / 2 / s ** 2) + np.exp(-(zz + c * t1 - b) ** 2 / 2 / s ** 2)) \
        - 0.5 * a * (np.exp(-(zz - c * (T - t1) - b) ** 2 / 2 / s ** 2) + np.exp(-(zz + c * (T - t1) - b) ** 2 / 2 / s ** 2))

    v_ex = c / 2 * a / s ** 2 * ((zz - c * t1 - b) * np.exp(-(zz - c * t1 - b) ** 2 / 2 / s ** 2)
                                  - (zz + c * t1 - b) * np.exp(-(zz + c * t1 - b) ** 2 / 2 / s ** 2)) \
        + c / 2 * a / s ** 2 * ((zz - c * (T - t1) - b) * np.exp(-(zz - c * (T - t1) - b) ** 2 / 2 / s ** 2)
                                 - (zz + c * (T - t1) - b) * np.exp(-(zz + c * (T - t1) - b) ** 2 / 2 / s ** 2))

    a_ex = 0.5 * a * (
        -c ** 2 / s ** 2 * np.exp(-0.5 * (zz - c * t1 - b) ** 2 / s ** 2)
        + 1 / s ** 4 * c ** 2 * ((zz - c * t1 - b) ** 2) * np.exp(-0.5 * (zz - c * t1 - b) ** 2 / s ** 2)
        - c ** 2 / s ** 2 * np.exp(-0.5 * (zz + c * t1 - b) ** 2 / s ** 2)
        + 1 / s ** 4 * c ** 2 * ((zz + c * t1 - b) ** 2) * np.exp(-0.5 * (zz + c * t1 - b) ** 2 / s ** 2)
    ) - 0.5 * a * (
        -c ** 2 / s ** 2 * np.exp(-0.5 * (zz - c * (T - t1) - b) ** 2 / s ** 2)
        + 1 / s ** 4 * c ** 2 * ((zz - c * (T - t1) - b) ** 2) * np.exp(-0.5 * (zz - c * (T - t1) - b) ** 2 / s ** 2)
        - c ** 2 / s ** 2 * np.exp(-0.5 * (zz + c * (T - t1) - b) ** 2 / s ** 2)
        + 1 / s ** 4 * c ** 2 * ((zz + c * (T - t1) - b) ** 2) * np.exp(-0.5 * (zz + c * (T - t1) - b) ** 2 / s ** 2)
    )

    return d_ex, v_ex, a_ex


# ------------------------------------------------------------------
# Read reference (undeformed) coordinates for each subdomain and find
# the line of nodes along z at the (min x, min y) corner -- this picks
# out the 1-D centerline of the bar within each subdomain. Uses the
# *first* segment's prefix, since a ROM phase reuses its FOM
# subdomain's mesh and does not write its own refe.csv.
# ------------------------------------------------------------------
x = {}
y = {}
z = {}
ind = {}
zline = {}
allz = []

for k in range(1, nsd + 1):
    refe_prefix = subdomain_segments[k][0][0]
    coords = np.loadtxt(f'{refe_prefix}-refe.csv', delimiter=',')
    xk, yk, zk = coords[:, 0], coords[:, 1], coords[:, 2]
    x[k], y[k], z[k] = xk, yk, zk
    idx = np.where((xk == xk.min()) & (yk == yk.min()))[0]
    ind[k] = idx
    zline[k] = zk[idx]
    allz.append(zk)

allz = np.concatenate(allz)
zglobal = np.unique(np.sort(allz))

# per-subdomain raw centerline time histories (kept across snapshots
# in case you want them after the loop, but not required for plotting)
dispz = {k: [] for k in range(1, nsd + 1)}
veloz = {k: [] for k in range(1, nsd + 1)}
accez = {k: [] for k in range(1, nsd + 1)}

# accumulators for the per-subdomain, per-field relative error
# (sum of squared norms, over all snapshots), computed against the
# exact solution evaluated directly at each subdomain's own z-coords
num_disp = {k: 0.0 for k in range(1, nsd + 1)}
den_disp = {k: 0.0 for k in range(1, nsd + 1)}
num_velo = {k: 0.0 for k in range(1, nsd + 1)}
den_velo = {k: 0.0 for k in range(1, nsd + 1)}
num_acce = {k: 0.0 for k in range(1, nsd + 1)}
den_acce = {k: 0.0 for k in range(1, nsd + 1)}

# ------------------------------------------------------------------
# Set up the figure (3 stacked subplots: disp, velo, acce)
# ------------------------------------------------------------------
plt.ion()
fig, axes = plt.subplots(3, 1, figsize=(8, 10))

# legend handles are rebuilt every snapshot (see inside the loop below),
# since each subdomain's label ("ROM1", "FOM1", ...) can change when it
# swaps models partway through the run. Only the "exact" handle is fixed.
exact_handle = mlines.Line2D([], [], color='c', linestyle='--', label='exact')

ctr = 0
for i in range(istart, iend + 1, istep):

    idx_str = str(i).zfill(pad_width)

    # ---- read this snapshot's centerline data for every subdomain ----
    t_snap = None
    dzk = {}
    vzk = {}
    azk = {}
    label_k = {}
    for k in range(1, nsd + 1):
        file_prefix = prefix_for_index(subdomain_segments[k], i)
        label_k[k] = model_type_label(file_prefix, k)

        dk = np.loadtxt(f'{file_prefix}-disp-{idx_str}.csv', delimiter=',')
        vk = np.loadtxt(f'{file_prefix}-velo-{idx_str}.csv', delimiter=',')
        ak = np.loadtxt(f'{file_prefix}-acce-{idx_str}.csv', delimiter=',')
        tk = np.loadtxt(f'{file_prefix}-time-{idx_str}.csv', delimiter=',')

        dispz[k].append(dk[:, 2])
        veloz[k].append(vk[:, 2])
        accez[k].append(ak[:, 2])

        dzk[k] = dispz[k][ctr][ind[k]]
        vzk[k] = veloz[k][ctr][ind[k]]
        azk[k] = accez[k][ctr][ind[k]]

        t_snap = float(tk)  # time should be identical across subdomains

    t1 = t_snap

    # legend handles for this snapshot: labels reflect each subdomain's
    # CURRENT model type ("ROM1", "FOM1", ...), so the legend visibly
    # updates the moment a subdomain swaps models.
    legend_handles = [
        mlines.Line2D([], [], color=colors[(k - 1) % len(colors)], linestyle='-',
                      label=label_k[k])
        for k in range(1, nsd + 1)
    ]
    legend_handles.append(exact_handle)

    # ---- exact / reference solution, evaluated on the union of
    #      subdomain z-coords (for plotting) ----
    zz = uniquetol(np.concatenate([zline[k] for k in range(1, nsd + 1)]), 1e-5)
    d_ex, v_ex, a_ex = exact_solution(zz, t1)

    # ---- exact / reference solution, evaluated at each subdomain's
    #      own z-coords (for the per-subdomain error accumulation) ----
    for k in range(1, nsd + 1):
        d_ex_k, v_ex_k, a_ex_k = exact_solution(zline[k], t1)

        num_disp[k] += np.linalg.norm(dzk[k] - d_ex_k) ** 2
        den_disp[k] += np.linalg.norm(d_ex_k) ** 2
        num_velo[k] += np.linalg.norm(vzk[k] - v_ex_k) ** 2
        den_velo[k] += np.linalg.norm(v_ex_k) ** 2
        num_acce[k] += np.linalg.norm(azk[k] - a_ex_k) ** 2
        den_acce[k] += np.linalg.norm(a_ex_k) ** 2

    # ---- plotting: each subdomain's raw solution overlaid with the
    #      exact/reference curve ----
    ax = axes[0]
    ax.clear()
    for k in range(1, nsd + 1):
        order = np.argsort(zline[k])
        ax.plot(zline[k][order], dzk[k][order], colors[(k - 1) % len(colors)] + '-')
    ax.plot(zz, d_ex, '--c')
    ax.set_xlabel('z')
    ax.set_ylabel('z-disp')
    ax.set_title(f'displacement snapshot {i + 1} at time = {t1}')
    ax.set_xlim(zglobal.min(), zglobal.max())
    ax.set_ylim(-scale, scale)
    ax.legend(handles=legend_handles, loc='upper right', fontsize='small')

    ax = axes[1]
    ax.clear()
    for k in range(1, nsd + 1):
        order = np.argsort(zline[k])
        ax.plot(zline[k][order], vzk[k][order], colors[(k - 1) % len(colors)] + '-')
    ax.plot(zz, v_ex, '--c')
    ax.set_xlabel('z')
    ax.set_ylabel('z-velo')
    ax.set_title(f'velocity snapshot {i + 1} at time = {t1}')
    ax.set_xlim(zglobal.min(), zglobal.max())
    ax.set_ylim(-3e4 * scale, 3e4 * scale)
    ax.legend(handles=legend_handles, loc='upper right', fontsize='small')

    ax = axes[2]
    ax.clear()
    for k in range(1, nsd + 1):
        order = np.argsort(zline[k])
        ax.plot(zline[k][order], azk[k][order], colors[(k - 1) % len(colors)] + '-')
    ax.plot(zz, a_ex, '--c')
    ax.set_xlabel('z')
    ax.set_ylabel('z-acce')
    ax.set_title(f'acceleration snapshot {i + 1} at time = {t1}')
    ax.set_xlim(zglobal.min(), zglobal.max())
    ax.set_ylim(-2.5e9 * scale, 2.5e9 * scale)
    ax.legend(handles=legend_handles, loc='upper right', fontsize='small')

    fig.tight_layout()
    plt.pause(pause_time)

    if save_figs:
        fig.savefig(f'soln_{ctr + 1:04d}.png', dpi=150)

    ctr += 1

plt.ioff()

# ------------------------------------------------------------------
# Relative errors (over all snapshots) between each subdomain's
# computed solution and the exact analytical solution, for all 3
# fields and all subdomains.
#
# NOTE: printed BEFORE the final plt.show() call below, since
# plt.show() blocks (waits for the plot window to be closed) on most
# backends/setups, which would otherwise delay or prevent this from
# ever being printed.
# ------------------------------------------------------------------
print()
print('Relative errors w.r.t. exact analytical solution:')
print('-' * 55)
for k in range(1, nsd + 1):
    disp_relerr = np.sqrt(num_disp[k] / den_disp[k])
    velo_relerr = np.sqrt(num_velo[k] / den_velo[k])
    acce_relerr = np.sqrt(num_acce[k] / den_acce[k])
    print(f'Subdomain {k}:')
    print(f'  z-disp avg rel error = {disp_relerr:.6e}')
    print(f'  z-velo avg rel error = {velo_relerr:.6e}')
    print(f'  z-acce avg rel error = {acce_relerr:.6e}')

plt.show()
