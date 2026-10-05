# Staged swaps-only adaptivity of the awful cube with a size field that
# follows the element count.
#
#   julia --project=<Norma> -t <threads> size-adaptive.jl [stages]
#
# stages is the number of runs of awful-cube-swaps.yaml (default 4), each one
# outer iteration of smoothing and swaps.  After each stage the script counts
# the elements, writes the mesh as awful-cube-stage-<k>.g, and sets the size
# field of the next stage to the edge length a of the regular tetrahedron
# that fills the cube with that many elements: the cube is 2 x 2 x 2 (volume
# 8) and a regular tetrahedron of edge a has volume a^3 / (6 sqrt 2), so
# n a^3 / (6 sqrt 2) = 8 gives a = cbrt(48 sqrt 2 / n).
using Norma
const stages = isempty(ARGS) ? 4 : parse(Int, ARGS[1])
const cube_volume = 8.0
params = Norma.load_input("awful-cube-swaps.yaml")
mesh = "awful-cube.g"
elem_size = 0.0
for stage in 1:stages
    stage_params = deepcopy(params)
    stage_params["input mesh file"] = mesh
    stage_params["output mesh file"] = "awful-cube-stage-$stage.e"
    if stage > 1
        stage_params["model"]["size field"] = "$elem_size"
    end
    sim = Norma.run(stage_params)
    topology = Norma.build_topology(sim.model)
    global mesh = "awful-cube-stage-$stage.g"
    num_elements = Norma.num_alive_elements(topology)
    global elem_size = cbrt(6 * sqrt(2) * cube_volume / num_elements)
    println("stage $stage: $num_elements elements, size field for the next stage $elem_size")
    Norma.write_topology(topology, mesh)
end
