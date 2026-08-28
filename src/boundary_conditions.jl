# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.


using Exodus

@variables t x y z
D = Differential(t)

# Swappable Dirichlet-Neumann Schwarz couplings.  Contact and non-overlap (DN)
# Schwarz are the couplings whose two sides carry complementary Dirichlet and
# Neumann roles (fields is_dirichlet / swap_bcs / dirichlet_projector /
# neumann_projector) and can exchange those roles; impedance (Robin) and overlap
# (pure Dirichlet) couplings cannot.  Contact sits outside the CouplingSchwarz
# subtree, so this capability cuts across the single-inheritance hierarchy and is
# expressed as a trait predicate rather than a shared supertype: it is the one
# place the {contact, non-overlap} pairing is defined, so the DN-swap, projector,
# and relaxation call sites test the trait instead of re-enumerating the types.
is_swappable_dn_schwarz(::SolidMechanicsBoundaryCondition) = false
is_swappable_dn_schwarz(::SolidMechanicsContactSchwarzBoundaryCondition) = true
is_swappable_dn_schwarz(::SolidMechanicsNonOverlapSchwarzBoundaryCondition) = true

function SolidMechanicsDirichletBoundaryCondition(input_mesh::ExodusDatabase, bc_params::Parameters)
    node_set_name = bc_params["node set"]
    expression = bc_params["function"]
    offset = component_offset_from_string(bc_params["component"])
    node_set_id = node_set_id_from_name(node_set_name, input_mesh)
    node_set_node_indices = Exodus.read_node_set_nodes(input_mesh, node_set_id)

    # Build symbolic expressions
    disp_num = eval(Meta.parse(expression))
    velo_num = expand_derivatives(D(disp_num))
    acce_num = expand_derivatives(D(velo_num))

    # Compile them into functions
    disp_fun = eval(build_function(disp_num, [t, x, y, z]; expression=Val(false)))
    velo_fun = eval(build_function(velo_num, [t, x, y, z]; expression=Val(false)))
    acce_fun = eval(build_function(acce_num, [t, x, y, z]; expression=Val(false)))

    return SolidMechanicsDirichletBoundaryCondition(
        node_set_name, offset, node_set_id, node_set_node_indices, disp_fun, velo_fun, acce_fun
    )
end

function SolidMechanicsSideSetDirichletBoundaryCondition(input_mesh::ExodusDatabase, bc_params::Parameters)
    side_set_name = bc_params["side set"]
    expression = bc_params["function"]
    offset = component_offset_from_string(bc_params["component"])
    side_set_id = side_set_id_from_name(side_set_name, input_mesh)
    # read_side_set_node_list repeats a node once per side it touches; the unique
    # node list is all we need to constrain the surface DOFs.
    _, side_set_node_indices = Exodus.read_side_set_node_list(input_mesh, side_set_id)
    node_indices = unique(Int64.(side_set_node_indices))

    # Build symbolic expressions
    disp_num = eval(Meta.parse(expression))
    velo_num = expand_derivatives(D(disp_num))
    acce_num = expand_derivatives(D(velo_num))

    # Compile them into functions
    disp_fun = eval(build_function(disp_num, [t, x, y, z]; expression=Val(false)))
    velo_fun = eval(build_function(velo_num, [t, x, y, z]; expression=Val(false)))
    acce_fun = eval(build_function(acce_num, [t, x, y, z]; expression=Val(false)))

    return SolidMechanicsSideSetDirichletBoundaryCondition(
        side_set_name, offset, side_set_id, node_indices, disp_fun, velo_fun, acce_fun
    )
end

function SolidMechanicsSurfaceBoundaryCondition(input_mesh::ExodusDatabase, bc_params::Parameters)
    side_set_name = bc_params["side set"]
    expression = string(bc_params["function"])
    enforcement_str = lowercase(string(get(bc_params, "enforcement", "exact")))
    if enforcement_str == "exact"
        enforcement = :exact
    elseif enforcement_str == "penalty"
        enforcement = :penalty
    else
        norma_abort(
            "Surface boundary condition on side set \"$side_set_name\": \"enforcement\" must be " *
            "\"exact\" or \"penalty\" (got \"$enforcement_str\").",
        )
    end
    penalty = Float64(get(bc_params, "penalty", 1.0e6))
    if enforcement == :penalty && penalty ≤ 0.0
        norma_abort("A penalty Surface boundary condition requires a positive \"penalty\" (got $penalty).")
    end
    side_set_id = side_set_id_from_name(side_set_name, input_mesh)
    # read_side_set_node_list repeats a node once per side it touches; the unique
    # node list is the set of DOFs the surface constraint acts on.
    _, side_set_node_indices = Exodus.read_side_set_node_list(input_mesh, side_set_id)
    node_indices = unique(Int64.(side_set_node_indices))

    # Compile g and its exact, automatic gradient ∇g from the symbolic surface
    # expression (constants are substituted into the string, exactly as for the
    # other analytic expressions).  Num() coerces a constant expression so the
    # degeneracy guard below still runs.
    g_num = Num(eval(Meta.parse(expression)))
    grad_num = Symbolics.gradient(g_num, [x, y, z])
    if all(iszero, grad_num)
        norma_abort(
            "Surface boundary condition on side set \"$side_set_name\" has a constant " *
            "level-set function \"$expression\" (∇g ≡ 0); it defines no surface to slide on.",
        )
    end
    level_set_fun = eval(build_function(g_num, [t, x, y, z]; expression=Val(false)))
    # build_function on a vector expression returns (out-of-place, in-place); the
    # out-of-place form returns ∇g as a 3-vector.
    level_set_grad = eval(build_function(grad_num, [t, x, y, z]; expression=Val(false))[1])

    return SolidMechanicsSurfaceBoundaryCondition(
        side_set_name, side_set_id, node_indices, level_set_fun, level_set_grad, enforcement, penalty
    )
end

function SolidMechanicsNeumannBoundaryCondition(input_mesh::ExodusDatabase, bc_params::Parameters)
    side_set_name = bc_params["side set"]
    expression = bc_params["function"]
    offset = component_offset_from_string(bc_params["component"])
    side_set_id = side_set_id_from_name(side_set_name, input_mesh)
    num_nodes_per_side, side_set_node_indices = Exodus.read_side_set_node_list(input_mesh, side_set_id)
    side_set_node_indices = Int64.(side_set_node_indices)

    # Build symbolic expressions
    traction_num = eval(Meta.parse(expression))

    # Compile them into functions
    traction_fun = eval(build_function(traction_num, [t, x, y, z]; expression=Val(false)))

    return SolidMechanicsNeumannBoundaryCondition(
        side_set_name, offset, side_set_id, num_nodes_per_side, side_set_node_indices, traction_fun
    )
end

function SolidMechanicsNeumannPressureBoundaryCondition(input_mesh::ExodusDatabase, bc_params::Parameters)
    side_set_name = bc_params["side set"]
    expression = bc_params["function"]
    side_set_id = side_set_id_from_name(side_set_name, input_mesh)
    num_nodes_per_side, side_set_node_indices = Exodus.read_side_set_node_list(input_mesh, side_set_id)
    side_set_node_indices = Int64.(side_set_node_indices)

    # Build symbolic expressions
    pressure_num = eval(Meta.parse(expression))

    # Compile them into functions
    pressure_fun = eval(build_function(pressure_num, [t, x, y, z]; expression=Val(false)))

    return SolidMechanicsNeumannPressureBoundaryCondition(
        side_set_name, side_set_id, num_nodes_per_side, side_set_node_indices, pressure_fun
    )
end

function SolidMechanicsRobinBoundaryCondition(input_mesh::ExodusDatabase, bc_params::Parameters)
    side_set_name = bc_params["side set"]
    expression = bc_params["function"]
    offset = component_offset_from_string(bc_params["component"])
    robin_parameter = bc_params["robin parameter"]
    tol = 1.0e-16
    if (abs(robin_parameter) < tol)
        norma_abort(
                "The robin parameter is close to zero.  Robin BC is equivalent " *
                "to Neumann BC. Use Neumann BC in input file.")
    end
    side_set_id = side_set_id_from_name(side_set_name, input_mesh)
    num_nodes_per_side, side_set_node_indices = Exodus.read_side_set_node_list(input_mesh, side_set_id)
    side_set_node_indices = Int64.(side_set_node_indices)

    # Build symbolic expressions
    rhs_num = eval(Meta.parse(expression))

    # Compile them into functions
    traction_fun = eval(build_function(rhs_num, [t, x, y, z]; expression=Val(false)))

    #We want to set traction + robin_parameter * disp = traction_fun
    return SolidMechanicsRobinBoundaryCondition(
        side_set_name, offset, side_set_id, num_nodes_per_side, side_set_node_indices,
        traction_fun, robin_parameter
    )
end


function SolidMechanicsOverlapSchwarzBoundaryCondition(
    coupled_block_name::String,
    tol::Float64,
    side_set_name::String,
    side_set_id::Int64,
    side_set_node_indices::Vector{Int64},
    num_nodes_sides::Vector{Int64},
    coupled_subsim::Simulation,
    subsim::Simulation,
    use_weak::Bool,
    compute_overlap_l2_error::String="",
)
    mesh = get_fom_model(subsim).mesh
    local_from_global_map = get_side_set_local_from_global_map(mesh, side_set_id)
    global_from_local_map = get_side_set_global_from_local_map(mesh, side_set_id)
    coupled_mesh = get_fom_model(coupled_subsim).mesh
    coupled_block_id = block_id_from_name(coupled_block_name, coupled_mesh)
    element_type_string = Exodus.read_block_parameters(coupled_mesh, coupled_block_id)[1]
    element_type = element_type_from_string(element_type_string)
    coupled_nodes_indices = Vector{Vector{Int64}}(undef, 0)
    interpolation_function_values = Vector{Vector{Float64}}(undef, 0)
    unique_node_indices = unique(side_set_node_indices)
    for node_index in unique_node_indices
        point = subsim.model.reference[:, node_index]
        node_indices, ξ, found = find_point_in_mesh(point, coupled_subsim.model, coupled_block_id, tol)
        if found == false
            norma_abortf(
                "Could not find subdomain %s point (%.4e, %.4e, %.4e) in subdomain %s",
                subsim.name,
                point[1],
                point[2],
                point[3],
                coupled_subsim.name,
            )
        end
        N = interpolate(element_type, ξ)[1]
        push!(coupled_nodes_indices, node_indices)
        push!(interpolation_function_values, N)
    end
    overlap_node_indices = Vector{Int64}(undef, 0)
    overlap_coupled_nodes_indices = Vector{Vector{Int64}}(undef, 0)
    overlap_interpolation_function_values = Vector{Vector{Float64}}(undef, 0)
    if !isempty(compute_overlap_l2_error)
        overlap_node_indices,
        overlap_coupled_nodes_indices,
        overlap_interpolation_function_values = build_overlap_l2_error_map(
            coupled_block_id,
            element_type,
            tol,
            coupled_subsim,
            subsim,
        )
    end
    dirichlet_projector = Matrix{Float64}(undef, 0, 0)
    return SolidMechanicsOverlapSchwarzBoundaryCondition(
        side_set_name,
        side_set_id,
        side_set_node_indices,
        num_nodes_sides,
        local_from_global_map,
        global_from_local_map,
        coupled_nodes_indices,
        interpolation_function_values,
        compute_overlap_l2_error,
        overlap_node_indices,
        overlap_coupled_nodes_indices,
        overlap_interpolation_function_values,
        NaN,
        coupled_block_name,
        tol,
        dirichlet_projector,
        use_weak,
        subsim.parent,
        subsim.handle,
        coupled_subsim.handle,
    )
end

function SolidMechanicsImpedanceOverlapSchwarzBoundaryCondition(
    coupled_block_name::String,
    tol::Float64,
    mesh::ExodusDatabase,
    side_set_name::String,
    side_set_id::Int64,
    side_set_node_indices::Vector{Int64},
    num_nodes_sides::Vector{Int64},
    coupled_subsim::Simulation,
    subsim::Simulation,
    impedance::Float64,
    impedance_shear::Float64,
    robin_parameter::Float64,
    impedance_scale::Vector{Float64},
    partner_traction_mode::String,
    transfer_mode::String,
    transfer_subdivisions::Int64,
    content_absorption::Bool,
    representable_dashpot::Bool,
)
    # Pointwise interpolation infrastructure (same as regular overlap)
    coupled_mesh = get_fom_model(coupled_subsim).mesh
    coupled_block_id = block_id_from_name(coupled_block_name, coupled_mesh)
    element_type_string = Exodus.read_block_parameters(coupled_mesh, coupled_block_id)[1]
    element_type = element_type_from_string(element_type_string)
    coupled_nodes_indices = Vector{Vector{Int64}}(undef, 0)
    interpolation_function_values = Vector{Vector{Float64}}(undef, 0)
    unique_node_indices = unique(side_set_node_indices)
    for node_index in unique_node_indices
        point = subsim.model.reference[:, node_index]
        node_indices, ξ, found = find_point_in_mesh(point, coupled_subsim.model, coupled_block_id, tol)
        if found == false
            norma_abortf(
                "Could not find subdomain %s point (%.4e, %.4e, %.4e) in subdomain %s",
                subsim.name, point[1], point[2], point[3], coupled_subsim.name,
            )
        end
        N = interpolate(element_type, ξ)[1]
        push!(coupled_nodes_indices, node_indices)
        push!(interpolation_function_values, N)
    end
    # Surface projector infrastructure (for weak force application)
    local_from_global_map = get_side_set_local_from_global_map(mesh, side_set_id)
    global_from_local_map = get_side_set_global_from_local_map(mesh, side_set_id)
    square_projector = Matrix{Float64}(undef, 0, 0)
    return SolidMechanicsImpedanceOverlapSchwarzBoundaryCondition(
        side_set_name,
        side_set_id,
        side_set_node_indices,
        num_nodes_sides,
        coupled_nodes_indices,
        interpolation_function_values,
        local_from_global_map,
        global_from_local_map,
        square_projector,
        impedance,
        impedance_shear,
        robin_parameter,
        impedance_scale,
        partner_traction_mode,
        nothing,
        transfer_mode,
        transfer_subdivisions,
        Matrix{Float64}(undef, 0, 0),
        content_absorption,
        Matrix{Float64}(undef, 0, 0),
        representable_dashpot,
        Matrix{Float64}(undef, 0, 0),
        coupled_block_name,
        tol,
        subsim.parent,
        subsim.handle,
        coupled_subsim.handle,
    )
end

function SolidMechanicsContactSchwarzBoundaryCondition(
    subsim::SingleDomainSimulation,
    coupled_subsim::SingleDomainSimulation,
    input_mesh::ExodusDatabase,
    bc_params::Parameters,
)
    side_set_name = bc_params["side set"]
    side_set_id = side_set_id_from_name(side_set_name, input_mesh)
    num_nodes_sides, side_set_node_indices = Exodus.read_side_set_node_list(input_mesh, side_set_id)
    side_set_node_indices = Int64.(side_set_node_indices)
    coupled_side_set_name = bc_params["source side set"]
    is_dirichlet = true
    dirichlet_projector = Matrix{Float64}(undef, 0, 0)
    neumann_projector = Matrix{Float64}(undef, 0, 0)
    local_from_global_map = get_side_set_local_from_global_map(input_mesh, side_set_id)
    global_from_local_map = get_side_set_global_from_local_map(input_mesh, side_set_id)
    coupled_bc_index = 0
    rotation_matrix = I(3)
    active_contact = false
    swap_bcs = get(bc_params, "swap BC types", false)
    friction_type_string = bc_params["friction type"]
    if friction_type_string == "frictionless"
        friction_type = 0
    elseif friction_type_string == "tied"
        friction_type = 1
    else
        norma_abort("Unknown or not implemented friction type : $friction_type_string")
    end
    return SolidMechanicsContactSchwarzBoundaryCondition(
        side_set_name,
        side_set_id,
        side_set_node_indices,
        num_nodes_sides,
        local_from_global_map,
        global_from_local_map,
        coupled_side_set_name,
        coupled_bc_index,
        dirichlet_projector,
        neumann_projector,
        is_dirichlet,
        swap_bcs,
        rotation_matrix,
        active_contact,
        friction_type,
        subsim.parent,
        subsim.handle,
        coupled_subsim.handle,
    )
end

function SolidMechanicsImpedanceNonOverlapSchwarzBoundaryCondition(
    mesh::ExodusDatabase,
    side_set_name::String,
    coupled_side_set_name::String,
    side_set_id::Int64,
    side_set_node_indices::Vector{Int64},
    num_nodes_sides::Vector{Int64},
    coupled_subsim::Simulation,
    subsim::Simulation,
    impedance::Float64,
    robin_parameter::Float64,
    adjoint_pairing::Bool,
)
    dirichlet_projector = Matrix{Float64}(undef, 0, 0)
    neumann_projector = Matrix{Float64}(undef, 0, 0)
    square_projector = Matrix{Float64}(undef, 0, 0)
    local_from_global_map = get_side_set_local_from_global_map(mesh, side_set_id)
    global_from_local_map = get_side_set_global_from_local_map(mesh, side_set_id)
    coupled_bc_index = 0
    return SolidMechanicsImpedanceNonOverlapSchwarzBoundaryCondition(
        side_set_name,
        side_set_id,
        side_set_node_indices,
        num_nodes_sides,
        local_from_global_map,
        global_from_local_map,
        coupled_side_set_name,
        coupled_bc_index,
        dirichlet_projector,
        neumann_projector,
        square_projector,
        impedance,
        robin_parameter,
        adjoint_pairing,
        subsim.parent,
        subsim.handle,
        coupled_subsim.handle,
    )
end

function SolidMechanicsNonOverlapSchwarzBoundaryCondition(
    mesh::ExodusDatabase,
    side_set_name::String,
    coupled_side_set_name::String,
    side_set_id::Int64,
    side_set_node_indices::Vector{Int64},
    num_nodes_sides::Vector{Int64},
    subsim::Simulation,
    coupled_subsim::Simulation,
    is_dirichlet::Bool,
    swap_bcs::Bool,
)
    dirichlet_projector = Matrix{Float64}(undef, 0, 0)
    neumann_projector = Matrix{Float64}(undef, 0, 0)
    square_projector = Matrix{Float64}(undef, 0, 0)
    local_from_global_map = get_side_set_local_from_global_map(mesh, side_set_id)
    global_from_local_map = get_side_set_global_from_local_map(mesh, side_set_id)
    coupled_bc_index = 0
    return SolidMechanicsNonOverlapSchwarzBoundaryCondition(
        side_set_name,
        side_set_id,
        side_set_node_indices,
        num_nodes_sides,
        local_from_global_map,
        global_from_local_map,
        coupled_side_set_name,
        coupled_bc_index,
        dirichlet_projector,
        neumann_projector,
        square_projector,
        is_dirichlet,
        swap_bcs,
        subsim.parent,
        subsim.handle,
        coupled_subsim.handle,
    )
end

function SMCouplingSchwarzBC(
    subsim::SingleDomainSimulation,
    coupled_subsim::SingleDomainSimulation,
    input_mesh::ExodusDatabase,
    bc_type::String,
    bc_params::Parameters,
)
    side_set_name = bc_params["side set"]
    coupled_block_name = get(bc_params, "source block", "")
    tol = get(bc_params, "search tolerance", 1.0e-06)
    # Only the nonoverlap variants below use a source side set; overlap variants
    # couple via a source block instead, so this key may be legitimately absent.
    coupled_side_set_name = get(bc_params, "source side set", "")
    side_set_id = side_set_id_from_name(side_set_name, input_mesh)
    num_nodes_sides, side_set_node_indices = Exodus.read_side_set_node_list(input_mesh, side_set_id)
    num_nodes_sides = Int64.(num_nodes_sides)
    side_set_node_indices = Int64.(side_set_node_indices)
    if bc_type == "Schwarz overlap"
        use_weak = get(bc_params, "weak", false)
        compute_overlap_l2_error = get(bc_params, "compute overlap L2 relative error", "")
        if !isempty(compute_overlap_l2_error) && compute_overlap_l2_error ∉ ("disp", "velo", "acce")
            norma_abort("Invalid value for 'compute overlap L2 relative error': \"$(compute_overlap_l2_error)\". Valid options are \"disp\", \"velo\", and \"acce\".")
        end
        SolidMechanicsOverlapSchwarzBoundaryCondition(
            coupled_block_name, tol, side_set_name, side_set_id, side_set_node_indices,
            num_nodes_sides, coupled_subsim, subsim, use_weak, compute_overlap_l2_error
        )
    elseif bc_type == "Schwarz DN nonoverlap"
        default_bc_type = get(bc_params, "default BC type", "Dirichlet")
        if default_bc_type == "Dirichlet"
            is_dirichlet = true
        elseif default_bc_type == "Neumann"
            is_dirichlet = false
        else
            norma_abort("Invalid string for 'default BC type'!  Valid options are 'Dirichlet' and 'Neumann'")
        end
        swap_bcs = get(bc_params, "swap BC types", false)
        SolidMechanicsNonOverlapSchwarzBoundaryCondition(
            input_mesh,
            side_set_name,
            coupled_side_set_name,
            side_set_id,
            side_set_node_indices,
            num_nodes_sides,
            subsim,
            coupled_subsim,
            is_dirichlet,
            swap_bcs,
        )
    elseif bc_type == "Schwarz RR nonoverlap" ||
           bc_type == "Schwarz impedance nonoverlap" ||
           bc_type == "Schwarz impedance overlap"
        # Characteristic impedances Z_p = √(ρ(λ + 2μ)) = ρ c_p and
        # Z_s = √(ρμ) = ρ c_s, computed from material properties.
        function _ps_impedances(sim)
            mat_params = sim.params["model"]["material"]
            mat_props = mat_params[first(values(mat_params["blocks"]))]
            E = Float64(mat_props["elastic modulus"])
            ν = Float64(mat_props["Poisson's ratio"])
            ρ = Float64(mat_props["density"])
            λ_lame = E * ν / ((1 + ν) * (1 - 2ν))
            μ = E / (2 * (1 + ν))
            return sqrt(ρ * (λ_lame + 2μ)), sqrt(ρ * μ)
        end
        # Nonoverlap variant: scalar P-impedance from this subdomain's own
        # material (unchanged legacy behavior). Overlap variant: P/S-split
        # tensor impedance from the NEIGHBOR's material, per the
        # optimized-Schwarz cross-scaling principle (each side's optimal
        # transmission operator approximates the neighbor's DtN map).
        robin_parameter = Float64(get(bc_params, "robin parameter", 0.0))
        raw_scale = get(bc_params, "impedance scale", 1.0)
        if raw_scale isa AbstractVector
            impedance_scale = Float64.(raw_scale)
        else
            impedance_scale = [Float64(raw_scale)]
        end
        if bc_type == "Schwarz RR nonoverlap"
            # Classical Robin-Robin transmission condition t + α W u = g: no
            # dashpot, and by default the legacy per-side transfer with
            # per-side Robin parameters, as in the classical Robin-Robin
            # literature. The condition has no absorbing mechanism and can
            # pump energy at the interface in elastodynamics (issue #176:
            # growing interface mode, lateral kink of the bar axis under pure
            # torsion), so it is intended for quasi-statics and for comparison
            # studies; `Schwarz impedance nonoverlap` is the recommended
            # condition for dynamics.
            if haskey(bc_params, "impedance scale")
                norma_abort(
                    "`impedance scale` is not a `Schwarz RR nonoverlap` key: this " *
                    "condition has no dashpot. For the impedance condition " *
                    "t + Z u̇ + α W u = g use `Schwarz impedance nonoverlap`.",
                )
            end
            if robin_parameter <= 0.0
                norma_abort(
                    "`Schwarz RR nonoverlap` requires a positive `robin parameter` " *
                    "(the Robin spring is this condition's only coupling term).",
                )
            end
            if get(get(subsim.params, "time integrator", Parameters()), "type", "") != "quasi static"
                norma_log(
                    0,
                    :warning,
                    "`Schwarz RR nonoverlap` applies the classical Robin condition " *
                    "t + α W u = g, which is not absorbing and can pump energy at " *
                    "the interface in elastodynamics (issue #176). For dynamics " *
                    "prefer `Schwarz impedance nonoverlap`.",
                )
            end
            # `adjoint pairing: true` remains available: it makes the Robin
            # spring a conservative interface spring built from the shared
            # cross-mass transfer, and requires ONE α per interface.
            adjoint_pairing = Bool(get(bc_params, "adjoint pairing", false))
            SolidMechanicsImpedanceNonOverlapSchwarzBoundaryCondition(
                input_mesh,
                side_set_name,
                coupled_side_set_name,
                side_set_id,
                side_set_node_indices,
                num_nodes_sides,
                coupled_subsim,
                subsim,
                0.0,
                robin_parameter,
                adjoint_pairing,
            )
        elseif bc_type == "Schwarz impedance nonoverlap"
            # Scalar dashpot scaling for the nonoverlap variant. A zero scale
            # would degenerate to the classical pure displacement-Robin
            # coupling, which has its own keyword; keep each keyword meaning
            # its condition.
            if length(impedance_scale) != 1
                norma_abort(
                    "`impedance scale` for `$bc_type` must be a scalar " *
                    "(the P/S-split schedule is an overlap-variant feature).",
                )
            end
            if impedance_scale[1] <= 0.0
                norma_abort(
                    "`impedance scale` must be positive for `$bc_type` " *
                    "(got $(impedance_scale[1])). For the classical pure Robin " *
                    "condition t + α W u = g (no dashpot) use `Schwarz RR nonoverlap`.",
                )
            end
            impedance, _ = _ps_impedances(subsim)
            impedance *= impedance_scale[1]
            # Adjoint pairing is the default: both sides derive their transfer
            # operators from one shared cross-mass matrix, share one impedance
            # and Robin parameter, and exchange the dynamically consistent
            # d'Alembert reaction (M·a + f_int − f_body). The cantilever
            # benchmark measured the legacy per-side transfer losing 9.6%/ms
            # even on conforming meshes (static reactions miss the interface
            # inertia) and injecting up to +363% on nonconforming ones, where
            # the paired coupling is sign-definite at every mesh ratio and
            # integrator combination (see docs/notes/schwarz-coupling).
            # The legacy behavior remains available as an explicit opt-out.
            adjoint_pairing = Bool(get(bc_params, "adjoint pairing", true))
            SolidMechanicsImpedanceNonOverlapSchwarzBoundaryCondition(
                input_mesh,
                side_set_name,
                coupled_side_set_name,
                side_set_id,
                side_set_node_indices,
                num_nodes_sides,
                coupled_subsim,
                subsim,
                impedance,
                robin_parameter,
                adjoint_pairing,
            )
        else
            impedance_p_coupled, impedance_s_coupled = _ps_impedances(coupled_subsim)
            coupled_block_name = bc_params["source block"]
            tol = Float64(get(bc_params, "search tolerance", 1.0e-06))
            partner_traction_mode = get(bc_params, "partner traction", "auto")
            if partner_traction_mode ∉ ("auto", "consistent traction", "recovered stress")
                norma_abort(
                    "Invalid `partner traction: $partner_traction_mode`. Valid values are " *
                    "\"auto\", \"consistent traction\", and \"recovered stress\".",
                )
            end
            # Variational (L2-projection) transfer is the default: it is the
            # contraction the nonmatching-grid Schwarz theory requires
            # (Gander-Halpern-Nataf 2003, Thm 7.4), and the cantilever
            # parametric study measured pointwise interpolation making the
            # impedance dashpot's interface work sign-indefinite on
            # nonconforming meshes (energy growth), where variational transfer
            # restores controlled dissipation (see
            # docs/notes/schwarz-coupling). On node-aligned interfaces
            # the two coincide. Pointwise remains an explicit legacy opt-in.
            transfer_mode = get(bc_params, "transfer", "variational")
            if transfer_mode ∉ ("pointwise", "variational")
                norma_abort(
                    "Invalid `transfer: $transfer_mode`. Valid values are " *
                    "\"pointwise\" and \"variational\".",
                )
            end
            transfer_subdivisions = Int64(get(bc_params, "transfer quadrature subdivisions", 1))
            if transfer_subdivisions < 1
                norma_abort("`transfer quadrature subdivisions` must be a positive integer.")
            end
            content_absorption = Bool(get(bc_params, "content aware absorption", false))
            representable_dashpot = Bool(get(bc_params, "representable dashpot", false))
            SolidMechanicsImpedanceOverlapSchwarzBoundaryCondition(
                coupled_block_name,
                tol,
                input_mesh,
                side_set_name,
                side_set_id,
                side_set_node_indices,
                num_nodes_sides,
                coupled_subsim,
                subsim,
                impedance_p_coupled,
                impedance_s_coupled,
                robin_parameter,
                impedance_scale,
                partner_traction_mode,
                transfer_mode,
                transfer_subdivisions,
                content_absorption,
                representable_dashpot,
            )
        end
    else
        norma_abort("Unknown boundary condition type : $bc_type")
    end
end

function apply_bc(model::SolidMechanics, bc::SolidMechanicsDirichletBoundaryCondition)
    for node_index in bc.node_set_node_indices
        txzy = (
            model.time, model.reference[1, node_index], model.reference[2, node_index], model.reference[3, node_index]
        )

        disp_val = bc.disp_fun(txzy)
        velo_val = bc.velo_fun(txzy)
        acce_val = bc.acce_fun(txzy)

        dof_index = 3 * (node_index - 1) + bc.offset
        model.displacement[bc.offset, node_index] = disp_val
        model.velocity[bc.offset, node_index] = velo_val
        model.acceleration[bc.offset, node_index] = acce_val
        model.free_dofs[dof_index] = false
    end
end

function apply_bc(model::SolidMechanics, bc::SolidMechanicsSideSetDirichletBoundaryCondition)
    for node_index in bc.node_indices
        txzy = (
            model.time, model.reference[1, node_index], model.reference[2, node_index], model.reference[3, node_index]
        )

        disp_val = bc.disp_fun(txzy)
        velo_val = bc.velo_fun(txzy)
        acce_val = bc.acce_fun(txzy)

        dof_index = 3 * (node_index - 1) + bc.offset
        model.displacement[bc.offset, node_index] = disp_val
        model.velocity[bc.offset, node_index] = velo_val
        model.acceleration[bc.offset, node_index] = acce_val
        model.free_dofs[dof_index] = false
    end
end

function apply_bc(model::SolidMechanics, bc::SolidMechanicsNeumannRobinBoundaryCondition)
    ss_node_index = 1
    for side in bc.num_nodes_per_side
        side_nodes = bc.side_set_node_indices[ss_node_index:(ss_node_index + side - 1)]
        side_coordinates = model.reference[:, side_nodes]
        nodal_force_component = get_side_set_nodal_forces(side_coordinates, bc.traction_fun, model.time)
        ss_node_index += side
        side_node_index = 1
        for node_index in side_nodes
            dof_index = 3 * (node_index - 1) + bc.offset
            model.boundary_force[dof_index] += nodal_force_component[side_node_index]
            side_node_index += 1
        end
    end
end

function apply_bc(model::SolidMechanics, bc::SolidMechanicsNeumannPressureBoundaryCondition)
    ss_node_index = 1
    for side in bc.num_nodes_per_side
        side_nodes = bc.side_set_node_indices[ss_node_index:(ss_node_index + side - 1)]
        side_coordinates = model.reference[:, side_nodes]
        nodal_force_component = get_side_set_nodal_pressure(side_coordinates, bc.pressure_fun, model.time)
        ss_node_index += side
        side_node_index = 1
        for node_index in side_nodes
            dof_indices = [3 * node_index - 2, 3 * node_index - 1, 3 * node_index]
            model.boundary_force[dof_indices] += nodal_force_component[:, side_node_index]
            side_node_index += 1
        end
    end
end

# The surface constraint neither prescribes DOFs nor sets a boundary force in the
# apply_bcs pass: it enters the residual as a penalty force built in
# apply_surface_penalty_internal_force! during evaluate().  The nodes remain free
# to slide, so nothing is done here.
function apply_bc(model::SolidMechanics, bc::SolidMechanicsSurfaceBoundaryCondition)
    return nothing
end

function apply_robin_bcs_internal_force!(model::SolidMechanics)
    for bc in model.boundary_conditions
        bc isa SolidMechanicsRobinBoundaryCondition || continue
        ss_node_index = 1
        for side in bc.num_nodes_per_side
            side_nodes = bc.side_set_node_indices[ss_node_index:(ss_node_index + side - 1)]
            side_coordinates = model.reference[:, side_nodes]
            K_side = bc.robin_parameter * get_side_set_nodal_stiffness(side_coordinates, model.time)
            ss_node_index += side
            for (i, node_i) in enumerate(side_nodes)
                dof_i = 3 * (node_i - 1) + bc.offset
                for (j, node_j) in enumerate(side_nodes)
                    u_j = model.displacement[bc.offset, node_j]
                    model.internal_force[dof_i] += K_side[i, j] * u_j
                end
            end
        end
    end
end

function build_robin_stiffness(model::SolidMechanics)
    num_nodes = size(model.reference, 2)
    num_dofs = 3 * num_nodes
    K_robin = spzeros(num_dofs, num_dofs)
    for bc in model.boundary_conditions
        bc isa SolidMechanicsRobinBoundaryCondition || continue
        ss_node_index = 1
        for side in bc.num_nodes_per_side
            side_nodes = bc.side_set_node_indices[ss_node_index:(ss_node_index + side - 1)]
            side_coordinates = model.reference[:, side_nodes]
            K_side = bc.robin_parameter * get_side_set_nodal_stiffness(side_coordinates, model.time)
            ss_node_index += side
            for (i, node_i) in enumerate(side_nodes)
                dof_i = 3 * (node_i - 1) + bc.offset
                for (j, node_j) in enumerate(side_nodes)
                    dof_j = 3 * (node_j - 1) + bc.offset
                    K_robin[dof_i, dof_j] += K_side[i, j]
                end
            end
        end
    end
    return K_robin
end

# Add the analytic surface penalty force ∂P/∂x = κ g ∇g (of the constraint
# P = ½ κ g²) to the internal force at each constrained node, with g and ∇g
# evaluated at the node's current position X = reference + displacement.  Being
# the gradient of ½ κ g², the force drives g → 0 (the node onto the surface)
# while leaving tangential motion free — a soft inclined roller.  Mirrors
# apply_robin_bcs_internal_force!: called after evaluate() has assembled the
# smoothing internal force.
function apply_surface_penalty_internal_force!(model::SolidMechanics)
    for bc in model.boundary_conditions
        (bc isa SolidMechanicsSurfaceBoundaryCondition && bc.enforcement == :penalty) || continue
        for node_index in bc.node_indices
            X = model.reference[:, node_index] + model.displacement[:, node_index]
            # A Vector argument (not a tuple) so the compiled vector-valued
            # gradient infers an Array output container rather than an NTuple.
            args = [model.time, X[1], X[2], X[3]]
            g = bc.level_set_fun(args)
            grad = bc.level_set_grad(args)
            factor = bc.penalty * g
            base = 3 * (node_index - 1)
            model.internal_force[base + 1] += factor * grad[1]
            model.internal_force[base + 2] += factor * grad[2]
            model.internal_force[base + 3] += factor * grad[3]
        end
    end
    return nothing
end

# Gauss–Newton tangent of the surface penalty, κ ∇g ∇gᵀ per constrained node
# (the g ∇²g term is dropped: it vanishes as g → 0).  Used by the Newton
# (HessianMinimizer) solver; the matrix-free EMS solvers need only the force.
function build_surface_penalty_stiffness(model::SolidMechanics)
    num_dofs = 3 * size(model.reference, 2)
    K_surface = spzeros(num_dofs, num_dofs)
    for bc in model.boundary_conditions
        (bc isa SolidMechanicsSurfaceBoundaryCondition && bc.enforcement == :penalty) || continue
        for node_index in bc.node_indices
            X = model.reference[:, node_index] + model.displacement[:, node_index]
            grad = bc.level_set_grad([model.time, X[1], X[2], X[3]])
            base = 3 * (node_index - 1)
            for a in 1:3
                for b in 1:3
                    K_surface[base + a, base + b] += bc.penalty * grad[a] * grad[b]
                end
            end
        end
    end
    return K_surface
end

# --- Exact local-frame roller (enforcement == :exact) --------------------------

# Closest-point Newton for the return-to-surface drift correction: stop once
# every constraint residual is below this, or after this many iterations.  On a
# planar (linear) surface one step is exact; a curved one needs a handful.
const SURFACE_RETURN_TOL = 1.0e-12
const SURFACE_RETURN_MAX_ITERS = 20
# Smallest singular value (∝ min |R| diagonal of the gradient QR) below which an
# edge/vertex is treated as degenerate: two (near-)parallel surfaces meet at a
# tangential, ill-posed intersection.
const SURFACE_FRAME_MIN_SV = 1.0e-9

# Gather, per constrained node, the (g, ∇g) closures of every exact-mode Surface
# BC whose side set contains it.  A node on two surfaces (a shared edge) collects
# both; three, a corner.  The joint set defines the node's normal subspace (the
# span of the gradients) and its tangent complement.  Returns an empty dict when
# no exact-mode surface BCs are present, so callers can early-return cheaply.
function collect_exact_surface_constraints(model::SolidMechanics)
    constraints = Dict{Int64,Vector{Tuple{Function,Function}}}()
    for bc in model.boundary_conditions
        (bc isa SolidMechanicsSurfaceBoundaryCondition && bc.enforcement == :exact) || continue
        for node_index in bc.node_indices
            push!(get!(constraints, node_index, Tuple{Function,Function}[]), (bc.level_set_fun, bc.level_set_grad))
        end
    end
    return constraints
end

# Orthonormal basis (3 × m) of a node's normal subspace — the span of its
# constraint gradients, evaluated at the current position — plus the constraint
# residual vector g and the gradient matrix A (m × 3, rows are gradients).
# Aborts if the gradients are (near-)linearly dependent (a degenerate edge or
# vertex, an input error).
function surface_normal_frame(cs::Vector{Tuple{Function,Function}}, x::AbstractVector{Float64}, time::Float64)
    m = length(cs)
    args = [time, x[1], x[2], x[3]]
    g = Vector{Float64}(undef, m)
    A = Matrix{Float64}(undef, m, 3)
    for (i, (gfun, dgfun)) in enumerate(cs)
        g[i] = gfun(args)
        A[i, :] = dgfun(args)
    end
    F = qr(permutedims(A))                       # QR of Aᵀ (3 × m): columns span the normals
    R = F.R
    if m > 3 || any(k -> abs(R[k, k]) < SURFACE_FRAME_MIN_SV, 1:m)
        norma_abort(
            "Exact surface constraint: $m surface gradients at a node are linearly dependent " *
            "(degenerate edge/vertex — nearly tangential surfaces). Check the surface definitions.",
        )
    end
    Q = Matrix(F.Q)[:, 1:m]
    return Q, g, A
end

# Closest-point-project each exact-constrained node back onto its surface(s):
# Newton iterations of g_i(x) = 0 taken in the normal subspace only (the
# minimum-norm step δx = -A⁺ g = -Aᵀ(A Aᵀ)⁻¹ g), so the node returns to the
# manifold without tangential drift.  Operates on model.displacement in place
# (x = reference + displacement).  Idempotent on an already-on-surface node.
function return_to_surface!(model::SolidMechanics)
    constraints = collect_exact_surface_constraints(model)
    isempty(constraints) && return nothing
    time = model.time
    for (node_index, cs) in constraints
        x = model.reference[:, node_index] + model.displacement[:, node_index]
        converged = false
        for _ in 1:SURFACE_RETURN_MAX_ITERS
            _, g, A = surface_normal_frame(cs, x, time)
            if maximum(abs, g) ≤ SURFACE_RETURN_TOL
                converged = true
                break
            end
            x .-= permutedims(A) * ((A * permutedims(A)) \ g)
        end
        if !converged
            norma_logf(
                4, :warning, "return_to_surface: node %d did not reach the surface tolerance", node_index
            )
        end
        model.displacement[:, node_index] = x - model.reference[:, node_index]
    end
    return nothing
end

# Zero the normal-subspace component of the residual at each exact-constrained
# node, so the matrix-free descent direction is purely tangential — the
# inclined-roller reaction carries the normal force.  Convergence is then
# measured on the tangential residual, the true first-order optimality condition
# for sliding on the surface.
function project_surface_gradient!(model::SolidMechanics, gradient::Vector{Float64})
    constraints = collect_exact_surface_constraints(model)
    isempty(constraints) && return nothing
    time = model.time
    for (node_index, cs) in constraints
        x = model.reference[:, node_index] + model.displacement[:, node_index]
        Q, _, _ = surface_normal_frame(cs, x, time)
        base = 3 * (node_index - 1)
        g_node = gradient[(base + 1):(base + 3)]
        g_node -= Q * (permutedims(Q) * g_node)
        gradient[(base + 1):(base + 3)] = g_node
    end
    return nothing
end

# Guard: exact enforcement is wired only through the matrix-free evaluate path
# (projected gradient + return to surface).  Refuse it on solvers whose evaluate
# path does not yet carry the projection, rather than silently ignore the
# constraint.  Penalty enforcement works on every solver and is the fallback.
function assert_no_exact_surface_bcs(model::SolidMechanics, solver_desc::String)
    for bc in model.boundary_conditions
        if bc isa SolidMechanicsSurfaceBoundaryCondition && bc.enforcement == :exact
            norma_abort(
                "Exact surface enforcement (enforcement: exact) requires a matrix-free solver " *
                "(steepest descent or lbfgs); $solver_desc is not supported. Use enforcement: penalty instead.",
            )
        end
    end
    return nothing
end


function compute_rotation_matrix(axis::SVector{3,Float64})::SMatrix{3,3,Float64}
    e1 = @SVector [1.0, 0.0, 0.0]
    angle_btwn = acos(dot(axis, e1))
    w = cross(axis, e1)
    s = clamp(norm(w), -1.0, 1.0)
    if isapprox(angle_btwn, 0.0; atol=1e-12)
        return @SMatrix [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    elseif isapprox(angle_btwn, π; atol=1e-12)
        return @SMatrix [
            -1.0 0.0 0.0
            0.0 -1.0 0.0
            0.0 0.0 1.0
        ]
    else
        θ = angle_btwn > π / 2 ? π - asin(s) : asin(s)
        m = normalize(w)
        rv = θ * m
        return rt_of_rv(rv)
    end
end


function collect_overlap_candidate_nodes(
    coupled_block_id::Int64,
    tol::Float64,
    coupled_subsim::Simulation,
    subsim::Simulation,
)
    dst_model = get_fom_model(subsim)
    dst_mesh = dst_model.mesh
    overlap_node_set = Set{Int64}()
    for block in Exodus.read_sets(dst_mesh, Block)
        block_id = block.id
        element_block_connectivity = get_block_connectivity(dst_mesh, block_id)
        num_block_elements, num_element_nodes = size(element_block_connectivity)
        for block_element_index in 1:num_block_elements
            connectivity_indices =
                ((block_element_index - 1) * num_element_nodes + 1):(block_element_index * num_element_nodes)
            node_indices = vec(element_block_connectivity[connectivity_indices])
            element_ref_pos = dst_model.reference[:, node_indices]
            centroid = vec(sum(element_ref_pos; dims=2) / size(element_ref_pos, 2))
            _, _, centroid_found = find_point_in_mesh(centroid, coupled_subsim.model, coupled_block_id, tol)
            node_found = false
            if !centroid_found
                for node_index in node_indices
                    _, _, node_found = find_point_in_mesh(dst_model.reference[:, node_index], coupled_subsim.model, coupled_block_id, tol)
                    node_found && break
                end
            end
            if centroid_found || node_found
                union!(overlap_node_set, node_indices)
            end
        end
    end
    return sort!(collect(overlap_node_set))
end

function build_overlap_l2_error_map(
    coupled_block_id::Int64,
    element_type::ElementType,
    tol::Float64,
    coupled_subsim::Simulation,
    subsim::Simulation,
)
    dst_model = get_fom_model(subsim)
    overlap_node_indices = collect_overlap_candidate_nodes(coupled_block_id, tol, coupled_subsim, subsim)
    overlap_coupled_nodes_indices = Vector{Vector{Int64}}(undef, 0)
    overlap_interpolation_function_values = Vector{Vector{Float64}}(undef, 0)
    mapped_overlap_node_indices = Vector{Int64}(undef, 0)
    for node_index in overlap_node_indices
        point = dst_model.reference[:, node_index]
        node_indices, ξ, found = find_point_in_mesh(point, coupled_subsim.model, coupled_block_id, tol)
        if !found
            continue
        end
        push!(mapped_overlap_node_indices, node_index)
        push!(overlap_coupled_nodes_indices, node_indices)
        push!(overlap_interpolation_function_values, interpolate(element_type, ξ)[1])
    end
    return mapped_overlap_node_indices, overlap_coupled_nodes_indices, overlap_interpolation_function_values
end

function compute_overlap_l2_error!(bc::SolidMechanicsOverlapSchwarzBoundaryCondition)
    if isempty(bc.compute_overlap_l2_error)
        bc.overlap_l2_error = NaN
        return bc.overlap_l2_error
    end

    src_model = get_fom_model(coupled_subsim_of(bc))
    dst_model = get_fom_model(self_subsim_of(bc))

    if bc.compute_overlap_l2_error == "disp"
        dst_field = dst_model.displacement
        src_field = src_model.displacement
    elseif bc.compute_overlap_l2_error == "velo"
        dst_field = dst_model.velocity
        src_field = src_model.velocity
    elseif bc.compute_overlap_l2_error == "acce"
        dst_field = dst_model.acceleration
        src_field = src_model.acceleration
    else
        error("Invalid value for compute_overlap_l2_error: \"$(bc.compute_overlap_l2_error)\". Valid options are \"disp\", \"velo\", and \"acce\".")
    end

    diff_sq = 0.0
    norm_sq = 0.0
    for i in eachindex(bc.overlap_node_indices)
        node_index = bc.overlap_node_indices[i]
        coupled_node_indices = bc.overlap_coupled_nodes_indices[i]
        N = bc.overlap_interpolation_function_values[i]
        dst_val = dst_field[:, node_index]
        src_val = src_field[:, coupled_node_indices] * N
        diff_sq += sum(abs2, dst_val - src_val)
        norm_sq += sum(abs2, dst_val)
    end
    bc.overlap_l2_error = norm_sq > 0.0 ? sqrt(diff_sq / norm_sq) : sqrt(diff_sq)
    return bc.overlap_l2_error
end

function compute_overlap_l2_error!(_)
    return NaN
end

function stores_overlap_l2_error(bc::SolidMechanicsOverlapSchwarzBoundaryCondition)
    return !isempty(bc.compute_overlap_l2_error)
end

function stores_overlap_l2_error(_)
    return false
end

function overlap_l2_error_field(bc::SolidMechanicsOverlapSchwarzBoundaryCondition)
    return bc.compute_overlap_l2_error
end

function get_overlap_l2_error(bc::SolidMechanicsOverlapSchwarzBoundaryCondition)
    return bc.overlap_l2_error
end

function get_overlap_l2_error(_)
    return NaN
end
