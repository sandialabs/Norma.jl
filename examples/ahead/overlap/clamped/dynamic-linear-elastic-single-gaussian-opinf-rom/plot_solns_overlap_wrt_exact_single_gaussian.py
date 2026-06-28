import numpy as np
import matplotlib.pyplot as plt


def uniquetol(values, tol=1.0e-5):
    """MATLAB-like uniquetol for a 1D array."""
    values = np.sort(np.asarray(values).ravel())
    if values.size == 0:
        return values

    unique_vals = [values[0]]
    for val in values[1:]:
        if abs(val - unique_vals[-1]) > tol:
            unique_vals.append(val)
    return np.asarray(unique_vals)


def snapshot_filename(domain, field, i):
    """Match MATLAB filenames such as clamped-1-disp-0010.csv."""
    return f"clamped-{domain}-{field}-{i:04d}.csv"


# -----------------------------------------------------------------------------
# Read coordinates
# -----------------------------------------------------------------------------
coords1 = np.loadtxt("clamped-1-refe.csv", delimiter=",")
x1 = coords1[:, 0]
y1 = coords1[:, 1]
z1 = coords1[:, 2]
ind1 = np.where((x1 == np.min(x1)) & (y1 == np.min(y1)))[0]

coords2 = np.loadtxt("clamped-2-refe.csv", delimiter=",")
x2 = coords2[:, 0]
y2 = coords2[:, 1]
z2 = coords2[:, 2]
ind2 = np.where((x2 == np.min(x2)) & (y2 == np.min(y2)))[0]

z = np.sort(np.unique(np.concatenate((z1, z2))))

# Storage arrays, kept as lists and stacked at the end.
dispz1 = []
veloz1 = []
accez1 = []
dispz2 = []
veloz2 = []
accez2 = []

disp_computed = []
velo_computed = []
acce_computed = []
disp_exact = []
velo_exact = []
acce_exact = []

# -----------------------------------------------------------------------------
# Problem parameters
# -----------------------------------------------------------------------------
c = np.sqrt(1.0e9 / 1.0e3)
a = 0.001
b = 0.0
s = 0.02
xL = -0.5

save_figs = True
plt.ion()
fig, axs = plt.subplots(3, 1, figsize=(8, 10))

for ctr, i in enumerate(range(0, 1001, 10), start=1):
    # -------------------------------------------------------------------------
    # Read snapshot data
    # -------------------------------------------------------------------------
    d1 = np.loadtxt(snapshot_filename(1, "disp", i), delimiter=",")
    v1 = np.loadtxt(snapshot_filename(1, "velo", i), delimiter=",")
    acc1 = np.loadtxt(snapshot_filename(1, "acce", i), delimiter=",")
    t1 = float(np.loadtxt(snapshot_filename(1, "time", i), delimiter=","))

    d2 = np.loadtxt(snapshot_filename(2, "disp", i), delimiter=",")
    v2 = np.loadtxt(snapshot_filename(2, "velo", i), delimiter=",")
    acc2 = np.loadtxt(snapshot_filename(2, "acce", i), delimiter="," )
    t2 = float(np.loadtxt(snapshot_filename(2, "time", i), delimiter=","))

    dispz1.append(d1[:, 2])
    veloz1.append(v1[:, 2])
    accez1.append(acc1[:, 2])
    dispz2.append(d2[:, 2])
    veloz2.append(v2[:, 2])
    accez2.append(acc2[:, 2])

    z1ind1 = z1[ind1]
    z2ind2 = z2[ind2]
    zz = uniquetol(np.concatenate((z1ind1, z2ind2)), tol=1.0e-5)

    dz1ind1 = dispz1[-1][ind1]
    dz2ind2 = dispz2[-1][ind2]
    vz1ind1 = veloz1[-1][ind1]
    vz2ind2 = veloz2[-1][ind2]
    az1ind1 = accez1[-1][ind1]
    az2ind2 = accez2[-1][ind2]

    dispz_merged = np.zeros_like(zz)
    veloz_merged = np.zeros_like(zz)
    accez_merged = np.zeros_like(zz)

    for k, zk in enumerate(zz):
        ii1 = np.where(np.isclose(z1ind1, zk, atol=1.0e-5, rtol=0.0))[0]
        ii2 = np.where(np.isclose(z2ind2, zk, atol=1.0e-5, rtol=0.0))[0]

        if ii1.size > 0:
            dispz_merged[k] += np.sum(dz1ind1[ii1])
            veloz_merged[k] += np.sum(vz1ind1[ii1])
            accez_merged[k] += np.sum(az1ind1[ii1])

        if ii2.size > 0:
            dispz_merged[k] += np.sum(dz2ind2[ii2])
            veloz_merged[k] += np.sum(vz2ind2[ii2])
            accez_merged[k] += np.sum(az2ind2[ii2])

        if ii1.size + ii2.size > 1:
            dispz_merged[k] /= 2.0
            veloz_merged[k] /= 2.0
            accez_merged[k] /= 2.0

    # -------------------------------------------------------------------------
    # Corrected exact solution: incident left-moving pulse plus odd reflection
    # across the homogeneous Dirichlet boundary x = xL = -0.5.
    # -------------------------------------------------------------------------
    xi_inc = zz + c * t1 - b
    xi_ref = 2.0 * xL - zz + c * t1 - b

    E_inc = np.exp(-0.5 * xi_inc**2 / s**2)
    E_ref = np.exp(-0.5 * xi_ref**2 / s**2)

    d_ex = a * E_inc - a * E_ref
    v_ex = -a * c / s**2 * xi_inc * E_inc + a * c / s**2 * xi_ref * E_ref
    a_ex = (
        a * c**2 * (xi_inc**2 / s**4 - 1.0 / s**2) * E_inc
        - a * c**2 * (xi_ref**2 / s**4 - 1.0 / s**2) * E_ref
    )

    disp_computed.append(dispz_merged)
    velo_computed.append(veloz_merged)
    acce_computed.append(accez_merged)
    disp_exact.append(d_ex)
    velo_exact.append(v_ex)
    acce_exact.append(a_ex)

    # -------------------------------------------------------------------------
    # Plot snapshots
    # -------------------------------------------------------------------------
    for ax in axs:
        ax.clear()

    axs[0].plot(z1ind1, dispz1[-1][ind1], "-b", label="Subdomain 1")
    axs[0].plot(z2ind2, dispz2[-1][ind2], "-r", label="Subdomain 2")
    axs[0].plot(zz, d_ex, "--c", label="Exact")
    axs[0].set_xlabel("z")
    axs[0].set_ylabel("z-disp")
    axs[0].set_title(f"displacement snapshot {i + 1} at time = {t1}")
    axs[0].axis([np.min(z), np.max(z), -0.001, 0.001])
    axs[0].legend(loc="best")

    axs[1].plot(z1ind1, veloz1[-1][ind1], "-b", label="Subdomain 1")
    axs[1].plot(z2ind2, veloz2[-1][ind2], "-r", label="Subdomain 2")
    axs[1].plot(zz, v_ex, "--c", label="Exact")
    axs[1].set_xlabel("z")
    axs[1].set_ylabel("z-velo")
    axs[1].set_title(f"velocity snapshot {i + 1} at time = {t1}")
    axs[1].axis([np.min(z), np.max(z), -3.0e4 * 0.001, 3.0e4 * 0.001])
    axs[1].legend(loc="best")

    axs[2].plot(z1ind1, accez1[-1][ind1], "-b", label="Subdomain 1")
    axs[2].plot(z2ind2, accez2[-1][ind2], "-r", label="Subdomain 2")
    axs[2].plot(zz, a_ex, "--c", label="Exact")
    axs[2].set_xlabel("z")
    axs[2].set_ylabel("z-acce")
    axs[2].set_title(f"acceleration snapshot {i + 1} at time = {t1}")
    axs[2].axis([np.min(z), np.max(z), -2.5e9 * 0.001, 2.5e9 * 0.001])
    axs[2].legend(loc="best")

    fig.tight_layout()
    plt.pause(0.5)

    if save_figs:
        fig.savefig(f"soln_{ctr:04d}.png", dpi=300)

# -----------------------------------------------------------------------------
# Relative errors over snapshots
# -----------------------------------------------------------------------------
numerator_disp = 0.0
denominator_disp = 0.0
numerator_velo = 0.0
denominator_velo = 0.0
numerator_acce = 0.0
denominator_acce = 0.0

for computed, exact in zip(disp_computed, disp_exact):
    numerator_disp += np.linalg.norm(computed - exact) ** 2
    denominator_disp += np.linalg.norm(exact) ** 2

for computed, exact in zip(velo_computed, velo_exact):
    numerator_velo += np.linalg.norm(computed - exact) ** 2
    denominator_velo += np.linalg.norm(exact) ** 2

for computed, exact in zip(acce_computed, acce_exact):
    numerator_acce += np.linalg.norm(computed - exact) ** 2
    denominator_acce += np.linalg.norm(exact) ** 2

dispz_relerr = np.sqrt(numerator_disp / denominator_disp)
veloz_relerr = np.sqrt(numerator_velo / denominator_velo)
accez_relerr = np.sqrt(numerator_acce / denominator_acce)

print(f"z-disp avg rel error = {dispz_relerr:e}")
print(f"z-velo avg rel error = {veloz_relerr:e}")
print(f"z-acce avg rel error = {accez_relerr:e}")

plt.ioff()
plt.show()
