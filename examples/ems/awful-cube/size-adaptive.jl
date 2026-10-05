# julia --project=~/Repos/Norma.jl -t [threads] size-aware.jl [stages]
# stages will be the number of outer loops

# The following setup aims to honor the thought experiment performed on the mesh
# when first understanding the size field smoothing reference. 
# The process will extract the number of elements after topology operations are 
# performed, and redefine the prescribed size under "size field" using the 
# equilateral tetrahedron volume formula and the volume of the cube. s
using Norma
const stages = isempty(ARGS) ? 4 : parse(Int, ARGS[1])
params = Norma.load_input("awful-cube-adaptive.yaml")
mesh = "awful-cube.g"
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
    global elem_size = cbrt(48*sqrt(2)/num_elements)
    Norma.write_topology(topology, mesh)
end