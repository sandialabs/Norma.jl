"""
Plot the computed solution (z-disp, z-velo, z-acce) for each of the
subdomains together with the exact reference solution, snapshot by
snapshot. Generalized to an arbitrary number of subdomains (here: 3
-> "clamped-1", "clamped-2", "clamped-3"). At the end, the relative
error (over all snapshots) between each subdomain's computed solution
and the exact analytical solution is printed for all 3 fields
(disp/velo/acce) and all 3 subdomains.

Each subdomain's CSV files follow the naming convention:
    <prefix>-<k>-refe.csv
    <prefix>-<k>-disp-XXXX.csv
    <prefix>-<k>-velo-XXXX.csv
    <prefix>-<k>-acce-XXXX.csv
    <prefix>-<k>-time-XXXX.csv
where k = 1..nsd and XXXX is a zero-padded snapshot index.

Run this script from the directory containing the CSV files.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

# ------------------------------------------------------------------
# Settings
# ------------------------------------------------------------------
prefix    = 'clamped'    # file name prefix, e.g. 'clamped' -> clamped-1-disp-0000.csv
nsd       = 3            # number of subdomains
istart    = 0             # first snapshot index
istep     = 1             # snapshot index increment
iend      = 100           # last snapshot index
pad_width = 4             # zero-padding width of the snapshot index, e.g. 0000

save_figs  = False
pause_time = 0.5
scale = 0.001 * 2

colors = ['b', 'r', 'g', 'm', 'k', 'y']  # cycle if nsd > 6


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
# out the 1-D centerline of the bar within each subdomain.
# ------------------------------------------------------------------
x = {}
y = {}
z = {}
ind = {}
zline = {}
allz = []

for k in range(1, nsd + 1):
    coords = np.loadtxt(f'{prefix}-{k}-refe.csv', delimiter=',')
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

# legend handles: one line per subdomain, plus the exact solution
legend_handles = [
    mlines.Line2D([], [], color=colors[(k - 1) % len(colors)], linestyle='-',
                  label=f'subdomain {k}')
    for k in range(1, nsd + 1)
]
legend_handles.append(
    mlines.Line2D([], [], color='c', linestyle='--', label='exact')
)

ctr = 0
for i in range(istart, iend + 1, istep):

    idx_str = str(i).zfill(pad_width)

    # ---- read this snapshot's centerline data for every subdomain ----
    t_snap = None
    dzk = {}
    vzk = {}
    azk = {}
    for k in range(1, nsd + 1):
        dk = np.loadtxt(f'{prefix}-{k}-disp-{idx_str}.csv', delimiter=',')
        vk = np.loadtxt(f'{prefix}-{k}-velo-{idx_str}.csv', delimiter=',')
        ak = np.loadtxt(f'{prefix}-{k}-acce-{idx_str}.csv', delimiter=',')
        tk = np.loadtxt(f'{prefix}-{k}-time-{idx_str}.csv', delimiter=',')

        dispz[k].append(dk[:, 2])
        veloz[k].append(vk[:, 2])
        accez[k].append(ak[:, 2])

        dzk[k] = dispz[k][ctr][ind[k]]
        vzk[k] = veloz[k][ctr][ind[k]]
        azk[k] = accez[k][ctr][ind[k]]

        t_snap = float(tk)  # time should be identical across subdomains

    t1 = t_snap

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
plt.show()

# ------------------------------------------------------------------
# Relative errors (over all snapshots) between each subdomain's
# computed solution and the exact analytical solution, for all 3
# fields and all 3 subdomains.
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
