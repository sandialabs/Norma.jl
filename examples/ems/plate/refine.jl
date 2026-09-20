# The four refinement stages of the sinusoid size field: the near-field size
# h_f is halved from 0.1 to 0.00625, each stage running the adaptivity loop
# of plate-sinusoid.yaml on the smoothed mesh of the stage before and writing
# its result as plate-stage-<k>.g.  Run from this directory:
#   julia --project=../../.. -t 8 refine.jl [stages]
using Norma

const h_c = 0.1
const stages = isempty(ARGS) ? 4 : parse(Int, ARGS[1])
params = Norma.load_input("plate-sinusoid.yaml")
mesh = "plate.g"
for stage in 1:stages
    h_f = h_c / 2^stage
    ratio = h_c / h_f
    stage_params = deepcopy(params)
    stage_params["input mesh file"] = mesh
    stage_params["output mesh file"] = "plate-stage-$stage.e"
    stage_params["model"]["size field"] = "$h_f * $ratio^(2 * abs(x - (0.1 * sin(4 * pi * y) + 0.5)))"
    Norma.norma_logf(0, :info, "Refinement stage %d: near-field size %.5f", stage, h_f)
    sim = Norma.run(stage_params)
    topology = Norma.build_topology(sim.model)
    global mesh = "plate-stage-$stage.g"
    Norma.write_topology(topology, mesh)
    Norma.norma_logf(
        0, :info, "Stage %d mesh %s: %d nodes, %d elements", stage, mesh,
        Norma.num_alive_nodes(topology), Norma.num_alive_elements(topology),
    )
end
