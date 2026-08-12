# In-memory restart: run one subdomain as a ROM for a window of steps, then
# switch it back -- no restart file, no rebuild.
#
# Set SUBDOMAIN and ROM_WINDOW below; everything else is the standard
# 2-subdomain clamped bar (sd1 604 nodes on z in [-0.50,-0.20], sd2 1604 nodes
# on [-0.30,0.50]). A Gaussian pulse starts at z=0 and travels left at
# c = sqrt(E/rho) = 1000 m/s, reaching sd1 at step 64 of 128.
#
#   julia --project=/path/to/Norma.jl restart_example1.jl
#   python3 plot_restart.py            # -> restart.gif

using Norma
using Printf
using DelimitedFiles

cd(@__DIR__)

# ---- the two inputs ----
# Only sd1 has a fitted operator here, so SUBDOMAIN = 2 needs one built first.
const SUBDOMAIN  = 1          # which subdomain switches: 1 = sd1, 2 = sd2
const ROM_WINDOW = (40, 80)   # it runs ROM over these steps, FOM outside them

const T0     = 0.0
const TFINAL = 4.0e-4
const NSTEP  = 128
const DT     = (TFINAL - T0) / NSTEP

# Shared mesh pool, named clamped-<N>sd-<i>.g
const MESH    = ("../meshes/clamped-2sd-1.g", "../meshes/clamped-2sd-2.g")
const BLOCK   = ("coarse", "fine")
const ZCLAMP  = ("nsz-", "nsz+")   # z-clamped end of each subdomain
const SIDESET = ("ssz+", "ssz-")   # Schwarz overlap face of each subdomain
const ROMFILE = "linear-opinf-operator-$SUBDOMAIN-M10.npz"

const IC = """
initial conditions:
  displacement:
    - node set: nsall
      component: z
      function: "a=0.001; s=0.02; a*exp(-z*z/s/s/2)"
  velocity:
    - node set: nsall
      component: z
      function: "a=0.001; s=0.02; E=1.0e+09; rho=1000.0; c=sqrt(E/rho); -a*c*z/s/s*exp(-z*z/s/s/2)"
"""

"The variants differ only in the `model:` block, so they discretize the same
problem and their states are mutually transferable."
function subdomain_yaml(i::Int, kind::Symbol)
    model = kind === :rom ? "  type: linear opinf rom\n  model-file: $ROMFILE" :
                            "  type: solid mechanics"
    return """
    type: single
    input mesh file: $(MESH[i])
    output mesh file: sd$i$(kind).e
    model:
    $model
      material:
        blocks:
          $(BLOCK[i]): hyperelastic
        hyperelastic:
          model: linear elastic
          elastic modulus: 1.0e+09
          Poisson's ratio: 0.0
          density: 1000.0
    time integrator:
      type: Newmark
      β: 0.25
      γ: 0.5
      time step: $DT
    $(IC)boundary conditions:
      Dirichlet:
        - {node set: nsx-, component: x, function: "0.0"}
        - {node set: nsx+, component: x, function: "0.0"}
        - {node set: nsy-, component: y, function: "0.0"}
        - {node set: nsy+, component: y, function: "0.0"}
        - {node set: $(ZCLAMP[i]), component: z, function: "0.0"}
      Schwarz overlap:
        - side set: $(SIDESET[i])
          source: sd$(3 - i)fom
          source block: $(BLOCK[3 - i])
    solver:
      type: Hessian minimizer
      step: full Newton
      minimum iterations: 1
      maximum iterations: 16
      relative tolerance: 1.0e-10
      absolute tolerance: 1.0e-06
    """
end

const TOP_YAML = """
type: multi
domains: ["sd1fom.yaml", "sd2fom.yaml"]
Exodus output interval: 1.0e+06
CSV output interval: 0
initial time: $T0
final time: $TFINAL
time step: $DT
minimum iterations: 1
maximum iterations: 16
relative tolerance: 1.0e-12
absolute tolerance: 1.0e-08
"""

isfile(ROMFILE) || error("no ROM operator $ROMFILE for sd$SUBDOMAIN; " *
                         "this example ships one for sd1 only")

foreach(rm, filter(f -> endswith(f, ".e") || endswith(f, ".log"), readdir()))
write("sd1fom.yaml", subdomain_yaml(1, :fom))
write("sd2fom.yaml", subdomain_yaml(2, :fom))
write("sd$(SUBDOMAIN)rom.yaml", subdomain_yaml(SUBDOMAIN, :rom))
write("top.yaml", TOP_YAML)

# ---- march ----
const FOM_NAME = "sd$(SUBDOMAIN)fom"
const ROM_NAME = "sd$(SUBDOMAIN)rom"
const A, B = ROM_WINDOW

log = open("restart_example1.log", "w")
r = Norma.InMemoryRestart("top.yaml")

# Build the ROM variant once; every later switch is a state transfer plus a
# pointer repoint -- no rebuild, no restart file.
Norma.register_variant!(r, SUBDOMAIN, ROM_NAME, "sd$(SUBDOMAIN)rom.yaml")

# The subdomains overlap, so average onto a shared z grid for output.
zs = [s.model.reference[3, :] for s in r.sim.subsims]
zgrid = sort(unique(round.(vcat(zs...), digits=9)))
zslot = [[searchsortedfirst(zgrid, round(z, digits=9)) for z in zs[i]] for i in eachindex(zs)]

"Full-order displacement of every subdomain, averaged onto `zgrid`."
function profile(r)
    acc = zeros(length(zgrid))
    cnt = zeros(Int, length(zgrid))
    for (i, s) in enumerate(r.sim.subsims)
        u = Norma.full_order_displacement(s)   # reconstructs a ROM's shadow FOM
        @inbounds for j in axes(u, 2)
            acc[zslot[i][j]] += u[3, j]
            cnt[zslot[i][j]] += 1
        end
    end
    return acc ./ max.(cnt, 1)
end

frames = Vector{Vector{Float64}}()
modes = String[]

redirect_stdout(log) do
    Norma.march!(r; on_step = function (r, step)
        # Decided live, one step at a time. A fixed window is just the simplest
        # thing to put in the hook; an error indicator or policy would work too.
        Norma.switch!(r, SUBDOMAIN, A <= step <= B ? ROM_NAME : FOM_NAME)
        push!(frames, profile(r))
        push!(modes, Norma.active_variant(r, SUBDOMAIN) == ROM_NAME ? "ROM" : "FOM")
    end)
end

writedlm("out_zgrid.csv", zgrid, ',')
writedlm("out_field.csv", reduce(hcat, frames)', ',')
writedlm("out_modes.csv",
         vcat(["step" "mode"], hcat(1:length(modes), modes)), ',')
Norma.close_restart!(r)
close(log)

@printf("sd%d ran ROM over steps %d-%d of %d\n", SUBDOMAIN, A, B, NSTEP)
for h in r.history
    @printf("  switch at step %3d: %s -> %s\n", h.step, h.from, h.to)
end
@printf("wrote out_field.csv (%d steps x %d z-points), out_zgrid.csv, out_modes.csv\n",
        length(frames), length(zgrid))
@printf("now run:  python3 plot_restart.py\n")
