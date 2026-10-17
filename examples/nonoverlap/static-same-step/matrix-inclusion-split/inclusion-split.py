# Cubit 17.08 Python script for the meshes of this example: a disk of radius
# 0.3 (the inclusion) in a 3 x 2 matrix, both one layer of tetrahedra of
# thickness 0.1, meshed independently so that the interface is nonconforming.
# Both bodies are cut by the planes x = 0.1 and y = -0.08 before meshing, so
# that the interface is available as one side set, arc, and as four side sets,
# arc_ll, arc_lr, arc_ul, and arc_ur, one per quadrant of the cut. The inputs
# couple the four side sets as four interfaces between the same two
# subdomains (issue #233). Run with the system python3:
#   python3 inclusion-split.py
import sys
sys.path.append('/usr/local/Cubit-17.08/bin')
import cubit
cubit.init(['cubit', '-batch', '-nographics', '-noecho', '-nojournal'])
R, H, XC, YC = 0.3, 0.1, 0.1, -0.08
def cmd(c):
    n = cubit.get_error_count(); cubit.cmd(c)
    if cubit.get_error_count() > n: raise RuntimeError(c)
quad = {'ll': (f'x_coord < {XC}', f'y_coord < {YC}'), 'lr': (f'x_coord > {XC}', f'y_coord < {YC}'),
        'ul': (f'x_coord < {XC}', f'y_coord > {YC}'), 'ur': (f'x_coord > {XC}', f'y_coord > {YC}')}
def arcs():
    out = {}
    for q, (cx, cy) in quad.items():
        ids = [s for s in cubit.parse_cubit_list('surface', f'all with {cx} and {cy}')
               if cubit.get_surface_type(s) != 'plane surface'
               and max(abs(c) for c in cubit.get_center_point('surface', s)[:2]) < 0.31]
        assert len(ids) == 1, (q, ids); out[q] = ids[0]
    return out
def split_and_name(first_ss):
    a = arcs()
    cmd(f'sideset {first_ss} surface {" ".join(map(str, a.values()))}'); cmd(f'sideset {first_ss} name "arc"')
    for i, (q, s) in enumerate(a.items()):
        cmd(f'sideset {first_ss+1+i} surface {s}'); cmd(f'sideset {first_ss+1+i} name "arc_{q}"')
cmd('reset')
cmd(f'create cylinder radius {R} z {H}')
cmd(f'webcut volume all with plane xplane offset {XC}'); cmd(f'webcut volume all with plane yplane offset {YC}')
cmd('imprint all'); cmd('merge all')
cmd('volume all size 0.07'); cmd('volume all scheme tetmesh'); cmd('mesh volume all')
cmd('block 1 volume all'); cmd('block 1 element type tetra4'); cmd('block 1 name "disk"')
split_and_name(1)
cmd('set large exodus file off'); cmd('export mesh "disk.g" overwrite')
cmd('reset')
cmd(f'create brick x 3 y 2 z {H}'); cmd(f'move volume 1 z {H/2}')
cmd(f'create cylinder radius {R} z {H}'); cmd(f'move volume 2 z {H/2}')
cmd('subtract volume 2 from volume 1')
cmd(f'webcut volume all with plane xplane offset {XC}'); cmd(f'webcut volume all with plane yplane offset {YC}')
cmd('imprint all'); cmd('merge all')
cmd('volume all size 0.14'); cmd('volume all scheme tetmesh'); cmd('mesh volume all')
cmd('block 1 volume all'); cmd('block 1 element type tetra4'); cmd('block 1 name "matrix"')
split_and_name(1)
outer = [s for s in cubit.parse_cubit_list('surface', 'all') if cubit.get_surface_type(s) == 'plane surface'
         and (abs(abs(cubit.get_center_point('surface', s)[0]) - 1.5) < 1e-6 or abs(abs(cubit.get_center_point('surface', s)[1]) - 1.0) < 1e-6)]
cmd(f'nodeset 1 surface {" ".join(map(str, outer))}'); cmd('nodeset 1 name "outer"')
cmd('set large exodus file off'); cmd('export mesh "matrix.g" overwrite')
print("MESHOK")
