#!/usr/bin/env python3
"""Meshes of test/schwarz-nonoverlap-matching-dn-symmetry.jl.

Two stacked boxes coupled by a Schwarz DN nonoverlap interface at z = H1:
  lower (Neumann side):   z in [0, H1],  NL cells per half-width in x and y
  upper (Dirichlet side): z in [H1, H],  NU cells per half-width in x and y
NL != NU, so the interface is nonconforming.

quarter/: x, y in [0, W]   -- the symmetry planes x = 0 and y = 0 cut the
                              interface
full/:    x, y in [-W, W]  -- exact mirror image of the quarter meshes, so the
                              full solution restricted to x, y >= 0 is the
                              reference for the quarter solution.

Sets: node sets nsall, nsx-, nsy- (quarter), nsz- (lower); side sets ssint,
ssx-, ssy- (quarter), ssz- (lower). Written as Exodus II through netCDF4.

Usage (from this directory): python3 make_meshes.py .
"""
import sys
from pathlib import Path

import netCDF4
import numpy as np

W, H1, H = 0.01, 0.01, 0.02
NL, NU = 3, 4  # cells per half-width (lower, upper)
NZL, NZU = 3, 4  # cells in z (lower, upper)

# Exodus HEX8 side -> local nodes (0-based)
HEX_SIDES = {1: [0, 1, 5, 4], 2: [1, 2, 6, 5], 3: [2, 3, 7, 6], 4: [0, 4, 7, 3], 5: [0, 3, 2, 1], 6: [4, 5, 6, 7]}


def box(xs, ys, zs):
    nx, ny, nz = len(xs), len(ys), len(zs)
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
    coords = np.c_[X.ravel(), Y.ravel(), Z.ravel()]
    nid = lambda i, j, k: (i * ny + j) * nz + k
    conn, cell = [], {}
    for i in range(nx - 1):
        for j in range(ny - 1):
            for k in range(nz - 1):
                cell[(i, j, k)] = len(conn)
                conn.append([nid(i, j, k), nid(i + 1, j, k), nid(i + 1, j + 1, k), nid(i, j + 1, k),
                             nid(i, j, k + 1), nid(i + 1, j, k + 1), nid(i + 1, j + 1, k + 1), nid(i, j + 1, k + 1)])
    return coords, np.array(conn), cell, (nx - 1, ny - 1, nz - 1)


def faces(cell, dims, which):
    nx, ny, nz = dims
    out = []
    for (i, j, k), e in cell.items():
        if which == "x-" and i == 0: out.append((e, 4))
        if which == "y-" and j == 0: out.append((e, 1))
        if which == "z-" and k == 0: out.append((e, 5))
        if which == "z+" and k == nz - 1: out.append((e, 6))
    return out


def write_exodus(path, title, coords, conn, block, node_sets, side_sets):
    path.unlink(missing_ok=True)
    ds = netCDF4.Dataset(path, "w", format="NETCDF3_64BIT_OFFSET")
    ds.api_version = np.float32(8.25)
    ds.version = np.float32(8.25)
    ds.floating_point_word_size = np.int32(8)
    ds.file_size = np.int32(1)
    ds.maximum_name_length = np.int32(32)
    ds.int64_status = np.int32(0)
    ds.title = title
    n, ne = len(coords), len(conn)
    for name, size in [("len_name", 256), ("time_step", None), ("num_dim", 3), ("num_nodes", n), ("num_elem", ne),
                       ("num_el_blk", 1), ("num_node_sets", len(node_sets)), ("num_side_sets", len(side_sets)),
                       ("num_el_in_blk1", ne), ("num_nod_per_el1", 8), ("four", 4), ("len_string", 33)]:
        ds.createDimension(name, size)

    def names_var(var, dim, names):
        v = ds.createVariable(var, "S1", (dim, "len_name"))
        arr = np.zeros((len(names), 256), dtype="S1")
        for r, nm in enumerate(names):
            arr[r, : len(nm)] = list(nm)
        v[:] = arr

    ds.createVariable("time_whole", "f8", ("time_step",))
    for pre, dim, count in [("eb", "num_el_blk", 1), ("ns", "num_node_sets", len(node_sets)),
                            ("ss", "num_side_sets", len(side_sets))]:
        ds.createVariable(f"{pre}_status", "i4", (dim,))[:] = np.ones(count, dtype="i4")
        p = ds.createVariable(f"{pre}_prop1", "i4", (dim,))
        p.setncattr("name", "ID")
        p[:] = np.arange(1, count + 1, dtype="i4")
    for c, a in zip("xyz", coords.T):
        ds.createVariable(f"coord{c}", "f8", ("num_nodes",))[:] = a
    names_var("eb_names", "num_el_blk", [block])
    names_var("ns_names", "num_node_sets", list(node_sets))
    names_var("ss_names", "num_side_sets", list(side_sets))
    names_var("coor_names", "num_dim", ["x", "y", "z"])
    cv = ds.createVariable("connect1", "i4", ("num_el_in_blk1", "num_nod_per_el1"))
    cv.elem_type = "HEX8"
    cv[:] = conn + 1
    for s, (nm, fl) in enumerate(side_sets.items(), start=1):
        ds.createDimension(f"num_side_ss{s}", len(fl))
        ds.createDimension(f"num_df_ss{s}", 4 * len(fl))
        ds.createVariable(f"elem_ss{s}", "i4", (f"num_side_ss{s}",))[:] = [e + 1 for e, _ in fl]
        ds.createVariable(f"side_ss{s}", "i4", (f"num_side_ss{s}",))[:] = [f for _, f in fl]
        ds.createVariable(f"dist_fact_ss{s}", "f8", (f"num_df_ss{s}",))[:] = np.ones(4 * len(fl))
    ds.createVariable("elem_map", "i4", ("num_elem",))[:] = np.arange(1, ne + 1)
    ds.createVariable("elem_num_map", "i4", ("num_elem",))[:] = np.arange(1, ne + 1)
    ds.createVariable("node_num_map", "i4", ("num_nodes",))[:] = np.arange(1, n + 1)
    for s, (nm, ids) in enumerate(node_sets.items(), start=1):
        ds.createDimension(f"num_nod_ns{s}", len(ids))
        ds.createVariable(f"node_ns{s}", "i4", (f"num_nod_ns{s}",))[:] = np.asarray(ids) + 1
        ds.createVariable(f"dist_fact_ns{s}", "f8", (f"num_nod_ns{s}",))[:] = np.ones(len(ids))
    ds.close()


def make(outdir, kind, part):
    ncell = NL if part == "lower" else NU
    nz = NZL if part == "lower" else NZU
    z0, z1 = (0.0, H1) if part == "lower" else (H1, H)
    if kind == "quarter":
        xs = np.linspace(0.0, W, ncell + 1)
    else:
        xs = np.linspace(-W, W, 2 * ncell + 1)
    xs[np.abs(xs) < 1e-15] = 0.0  # exact symmetry plane
    zs = np.linspace(z0, z1, nz + 1)
    coords, conn, cell, dims = box(xs, xs, zs)
    x, y, z = coords.T
    tol = 1e-12
    node_sets = {"nsall": np.arange(len(coords))}
    side_sets = {}
    if kind == "quarter":
        node_sets["nsx-"] = np.where(np.abs(x) < tol)[0]
        node_sets["nsy-"] = np.where(np.abs(y) < tol)[0]
        side_sets["ssx-"] = faces(cell, dims, "x-")
        side_sets["ssy-"] = faces(cell, dims, "y-")
    if part == "lower":
        node_sets["nsz-"] = np.where(np.abs(z) < tol)[0]
        side_sets["ssz-"] = faces(cell, dims, "z-")
        side_sets["ssint"] = faces(cell, dims, "z+")
    else:
        side_sets["ssint"] = faces(cell, dims, "z-")
    block = "steel" if part == "lower" else "titanium"
    d = Path(outdir) / kind
    d.mkdir(parents=True, exist_ok=True)
    write_exodus(d / f"{part}.g", f"held-dof verification {kind} {part}", coords, conn, block, node_sets, side_sets)


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else "meshes"
    for kind in ("quarter", "full"):
        for part in ("lower", "upper"):
            make(out, kind, part)
    print("wrote meshes to", out)
