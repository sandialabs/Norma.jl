# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file LICENSE.md in the
# top-level Norma.jl directory.

# The curved 10-node tetrahedral mesh of a cylinder (radius 0.5, height 2) of
# the composite tet LCM comparison: gmsh 4.15 puts the edge nodes of the
# boundary on the cylinder. Writes coords.txt and conn.txt (Exodus node
# order, 1-based) for write_cylinder.jl:
#
#   python3 make_cylinder.py 0.3 && julia --project=../../../.. write_cylinder.jl cylinder-tet10.g
import gmsh, sys, numpy as np
R, H, h = 0.5, 2.0, float(sys.argv[1]) if len(sys.argv) > 1 else 0.3
gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0)
gmsh.model.add("cyl")
gmsh.model.occ.addCylinder(0, 0, 0, 0, 0, H, R)
gmsh.model.occ.synchronize()
gmsh.option.setNumber("Mesh.MeshSizeMax", h); gmsh.option.setNumber("Mesh.MeshSizeMin", h)
gmsh.option.setNumber("Mesh.ElementOrder", 2)
gmsh.option.setNumber("Mesh.HighOrderOptimize", 0)
gmsh.model.mesh.generate(3)
tags, x, _ = gmsh.model.mesh.getNodes()
x = x.reshape(-1, 3)
idx = {t: i for i, t in enumerate(tags)}
et, etags, enodes = gmsh.model.mesh.getElements(3)
k = list(et).index(11)  # 10-node tet
conn = np.array([idx[t] for t in enodes[k]]).reshape(-1, 10)
conn = conn[:, [0, 1, 2, 3, 4, 5, 6, 7, 9, 8]]  # gmsh -> exodus edge order
used = np.unique(conn); remap = -np.ones(len(x), int); remap[used] = np.arange(len(used))
x = x[used]; conn = remap[conn]
np.savetxt("coords.txt", x, fmt="%.17g"); np.savetxt("conn.txt", conn + 1, fmt="%d")
print(len(x), "nodes", len(conn), "elements")
gmsh.finalize()
