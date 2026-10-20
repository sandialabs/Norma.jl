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
# neumann_projector) and can exchange those roles; Robin-Robin and overlap
# (pure Dirichlet) couplings cannot.  Contact sits outside the CouplingSchwarz
# subtree, so this capability cuts across the single-inheritance hierarchy and is
# expressed as a trait predicate rather than a shared supertype: it is the one
# place the {contact, non-overlap} pairing is defined, so the DN-swap, projector,
# and relaxation call sites test the trait instead of re-enumerating the types.
is_swappable_dn_schwarz(::SolidMechanicsBoundaryCondition) = false
is_swappable_dn_schwarz(::SolidMechanicsContactSchwarzBoundaryCondition) = true
is_swappable_dn_schwarz(::SolidMechanicsNonOverlapSchwarzBoundaryCondition) = true
# The overlap matching coupling relaxes its constrained datum as the Dirichlet
# side of a nonoverlap matching pair does; both of its sides are Dirichlet.
is_swappable_dn_schwarz(::SolidMechanicsOverlapMatchingSchwarzBoundaryCondition) = true

# A prescription is static when its time derivative is identically zero.
function is_static_prescription(velocity_expression)
    value = Symbolics.value(velocity_expression)
    return value isa Number && iszero(value)
end

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
        node_set_name,
        offset,
        node_set_id,
        node_set_node_indices,
        disp_fun,
        velo_fun,
        acce_fun,
        !is_static_prescription(velo_num),
        Float64[],
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
        side_set_name,
        offset,
        side_set_id,
        node_indices,
        disp_fun,
        velo_fun,
        acce_fun,
        !is_static_prescription(velo_num),
        Float64[],
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
            "to Neumann BC. Use Neumann BC in input file.",
        )
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
        side_set_name, offset, side_set_id, num_nodes_per_side, side_set_node_indices, traction_fun, robin_parameter
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
    local_from_global_map, global_from_local_map = get_side_set_node_maps(mesh, side_set_id)
    coupled_mesh = get_fom_model(coupled_subsim).mesh
    coupled_block_id = block_id_from_name(coupled_block_name, coupled_mesh)
    element_type = block_element_type(get_fom_model(coupled_subsim), coupled_block_id)
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
        overlap_node_indices, overlap_coupled_nodes_indices, overlap_interpolation_function_values = build_overlap_l2_error_map(
            coupled_block_id, element_type, tol, coupled_subsim, subsim
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
    local_from_global_map, global_from_local_map = get_side_set_node_maps(input_mesh, side_set_id)
    coupled_bc_index = 0
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
        active_contact,
        friction_type,
        subsim.parent,
        subsim.handle,
        coupled_subsim.handle,
    )
end

function SolidMechanicsRobinNonOverlapSchwarzBoundaryCondition(
    mesh::ExodusDatabase,
    side_set_name::String,
    coupled_side_set_name::String,
    side_set_id::Int64,
    side_set_node_indices::Vector{Int64},
    num_nodes_sides::Vector{Int64},
    coupled_subsim::Simulation,
    subsim::Simulation,
    robin_parameter::Float64,
)
    dirichlet_projector = Matrix{Float64}(undef, 0, 0)
    neumann_projector = Matrix{Float64}(undef, 0, 0)
    square_projector = Matrix{Float64}(undef, 0, 0)
    local_from_global_map, global_from_local_map = get_side_set_node_maps(mesh, side_set_id)
    coupled_bc_index = 0
    return SolidMechanicsRobinNonOverlapSchwarzBoundaryCondition(
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
        robin_parameter,
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
    constrained::Bool=false,
    constraint::Symbol=:unset,
    direct_solve::Bool=false,
    window_end_constraint::Bool=false,
)
    dirichlet_projector = Matrix{Float64}(undef, 0, 0)
    neumann_projector = Matrix{Float64}(undef, 0, 0)
    square_projector = Matrix{Float64}(undef, 0, 0)
    local_from_global_map, global_from_local_map = get_side_set_node_maps(mesh, side_set_id)
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
        nothing,
        is_dirichlet,
        swap_bcs,
        constrained,
        constraint,
        direct_solve,
        Float64[],
        Tuple{Float64,Vector{Float64}}[],
        window_end_constraint,
        Float64[],
        Float64[],
        Float64[],
        Float64[],
        subsim.parent,
        subsim.handle,
        coupled_subsim.handle,
    )
end

# The value of `interface solve`: `iterative` (default) or `direct`.
function parse_direct_solve(bc_params::Parameters, side_set_name::String)
    value = String(get(bc_params, "interface solve", "iterative"))
    value == "iterative" && return false
    value == "direct" && return true
    return norma_abort(
        "Invalid `interface solve: $(value)` on side set \"$(side_set_name)\". " *
        "Valid values are `iterative` and `direct`.",
    )
end

# The value of `matching at`: `every substep` (default) or `window end`.
function parse_window_end_constraint(bc_params::Parameters, side_set_name::String)
    value = String(get(bc_params, "matching at", "every substep"))
    value == "every substep" && return false
    value == "window end" && return true
    return norma_abort(
        "Invalid `matching at: $(value)` on side set \"$(side_set_name)\". " *
        "Valid values are `every substep` and `window end`.",
    )
end

function SolidMechanicsOverlapMatchingSchwarzBoundaryCondition(
    mesh::ExodusDatabase,
    side_set_name::String,
    coupled_block_name::String,
    tol::Float64,
    side_set_id::Int64,
    side_set_node_indices::Vector{Int64},
    num_nodes_sides::Vector{Int64},
    subsim::Simulation,
    coupled_subsim::Simulation,
    constraint::Symbol,
    window_end_constraint::Bool,
    direct_solve::Bool=false,
)
    local_from_global_map, global_from_local_map = get_side_set_node_maps(mesh, side_set_id)
    empty = Matrix{Float64}(undef, 0, 0)
    return SolidMechanicsOverlapMatchingSchwarzBoundaryCondition(
        side_set_name,
        side_set_id,
        side_set_node_indices,
        num_nodes_sides,
        local_from_global_map,
        global_from_local_map,
        coupled_block_name,
        tol,
        0,
        empty,
        empty,
        empty,
        nothing,
        true,
        false,
        true,
        constraint,
        direct_solve,
        Float64[],
        Tuple{Float64,Vector{Float64}}[],
        window_end_constraint,
        Float64[],
        Float64[],
        Float64[],
        Float64[],
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
            norma_abort(
                "Invalid value for 'compute overlap L2 relative error': \"$(compute_overlap_l2_error)\". Valid options are \"disp\", \"velo\", and \"acce\".",
            )
        end
        SolidMechanicsOverlapSchwarzBoundaryCondition(
            coupled_block_name,
            tol,
            side_set_name,
            side_set_id,
            side_set_node_indices,
            num_nodes_sides,
            coupled_subsim,
            subsim,
            use_weak,
            compute_overlap_l2_error,
        )
    elseif bc_type == "Schwarz overlap matching"
        isempty(coupled_block_name) && norma_abort("`Schwarz overlap matching` requires `source block`.")
        matching_name = String(get(bc_params, "matching", "velocity"))
        constraint = if matching_name == "velocity"
            :velocity
        elseif matching_name == "displacement"
            :displacement
        else
            norma_abort(
                "Invalid `matching: $(matching_name)` on side set \"$(side_set_name)\". " *
                "Valid values are `velocity` and `displacement`.",
            )
        end
        SolidMechanicsOverlapMatchingSchwarzBoundaryCondition(
            input_mesh,
            side_set_name,
            coupled_block_name,
            tol,
            side_set_id,
            side_set_node_indices,
            num_nodes_sides,
            subsim,
            coupled_subsim,
            constraint,
            parse_window_end_constraint(bc_params, side_set_name),
            parse_direct_solve(bc_params, side_set_name),
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
        for retired in ("constrained", "constraint")
            haskey(bc_params, retired) && norma_abort(
                "`$(retired)` on side set \"$(side_set_name)\" was retired. The matching exchange is " *
                "selected by `matching: velocity` or `matching: displacement` on both sides of the pair.",
            )
        end
        # `matching` selects the matching exchange and its matched quantity;
        # without it the pair uses the standard exchange.
        constraint = if haskey(bc_params, "matching")
            matching_name = String(bc_params["matching"])
            if matching_name == "velocity"
                :velocity
            elseif matching_name == "displacement"
                :displacement
            else
                norma_abort(
                    "Invalid `matching: $(matching_name)` on side set \"$(side_set_name)\". " *
                    "Valid values are `velocity` and `displacement`.",
                )
            end
        else
            :unset
        end
        constrained = constraint != :unset
        interface_solve = String(get(bc_params, "interface solve", "iterative"))
        if interface_solve ∉ ("iterative", "direct")
            norma_abort(
                "Invalid `interface solve: $(interface_solve)` on side set \"$(side_set_name)\". " *
                "Valid values are `iterative` and `direct`.",
            )
        end
        direct_solve = interface_solve == "direct"
        if direct_solve && !constrained
            norma_abort("`interface solve: direct` on side set \"$(side_set_name)\" requires `matching`.")
        end
        if constrained && swap_bcs
            norma_abort(
                "`matching` on side set \"$(side_set_name)\" does not support " *
                "`swap BC types: true`: the matching exchange keeps the Dirichlet side fixed.",
            )
        end
        window_end_constraint = parse_window_end_constraint(bc_params, side_set_name)
        if window_end_constraint && !constrained
            norma_abort("`matching at: window end` on side set \"$(side_set_name)\" requires `matching`.")
        end
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
            constrained,
            constraint,
            direct_solve,
            window_end_constraint,
        )
    elseif bc_type == "Schwarz RR nonoverlap"
        # Classical Robin-Robin transmission condition t + α W u = g, with
        # per-side transfer operators and per-side Robin parameters. The
        # condition does not absorb interface waves and can inject energy at the
        # interface in elastodynamics (issue #176: a growing interface mode,
        # with a lateral kink of the bar axis under pure torsion), so it is
        # intended for quasi-statics and for comparison studies. For dynamics the
        # matching Dirichlet-Neumann exchange (`Schwarz DN nonoverlap` with
        # `matching: velocity`) is recommended.
        robin_parameter = Float64(get(bc_params, "robin parameter", 0.0))
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
                "t + α W u = g, which is not absorbing and can inject energy at " *
                "the interface in elastodynamics (issue #176). For dynamics " *
                "prefer `Schwarz DN nonoverlap` with `matching: velocity`.",
            )
        end
        SolidMechanicsRobinNonOverlapSchwarzBoundaryCondition(
            input_mesh,
            side_set_name,
            coupled_side_set_name,
            side_set_id,
            side_set_node_indices,
            num_nodes_sides,
            coupled_subsim,
            subsim,
            robin_parameter,
        )
    else
        norma_abort("Unknown boundary condition type : $bc_type")
    end
end

# Prescribe the constrained component at the given nodes. Expressions that do
# not involve time are evaluated once and cached; the compiled symbolic
# functions are not inferable and cost several allocations per call.
function apply_dirichlet_values!(model::SolidMechanics, node_indices::Vector{Int64}, bc)
    offset = bc.offset
    if bc.time_dependent == false
        cache = bc.cached_displacement
        if length(cache) != length(node_indices)
            resize!(cache, length(node_indices))
            for (k, node_index) in enumerate(node_indices)
                txzy = (
                    model.time,
                    model.reference[1, node_index],
                    model.reference[2, node_index],
                    model.reference[3, node_index],
                )
                cache[k] = bc.disp_fun(txzy)
            end
        end
        @inbounds for (k, node_index) in enumerate(node_indices)
            model.displacement[offset, node_index] = cache[k]
            model.velocity[offset, node_index] = 0.0
            model.acceleration[offset, node_index] = 0.0
            model.free_dofs[3 * (node_index - 1) + offset] = false
        end
        return nothing
    end
    for node_index in node_indices
        txzy = (
            model.time, model.reference[1, node_index], model.reference[2, node_index], model.reference[3, node_index]
        )
        disp_val = bc.disp_fun(txzy)
        velo_val = bc.velo_fun(txzy)
        acce_val = bc.acce_fun(txzy)
        model.displacement[offset, node_index] = disp_val
        model.velocity[offset, node_index] = velo_val
        model.acceleration[offset, node_index] = acce_val
        model.free_dofs[3 * (node_index - 1) + offset] = false
    end
    return nothing
end

function apply_bc(model::SolidMechanics, bc::SolidMechanicsDirichletBoundaryCondition)
    return apply_dirichlet_values!(model, bc.node_set_node_indices, bc)
end

function apply_bc(model::SolidMechanics, bc::SolidMechanicsSideSetDirichletBoundaryCondition)
    return apply_dirichlet_values!(model, bc.node_indices, bc)
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

# The spring forces of the Robin conditions, added to `force`, by default the
# internal force.
function apply_robin_bcs_internal_force!(model::SolidMechanics, force::AbstractVector{Float64}=model.internal_force)
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
                    force[dof_i] += K_side[i, j] * u_j
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
# P = ½ κ g²) to `force`, by default the internal force, at each constrained
# node, with g and ∇g evaluated at the node's current position
# X = reference + displacement.  Being the gradient of ½ κ g², the force drives
# g → 0 (the node onto the surface) while leaving tangential motion free — a
# soft inclined roller.  Mirrors apply_robin_bcs_internal_force!: called after
# evaluate() has assembled the smoothing internal force.
function apply_surface_penalty_internal_force!(
    model::SolidMechanics, force::AbstractVector{Float64}=model.internal_force
)
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
            force[base + 1] += factor * grad[1]
            force[base + 2] += factor * grad[2]
            force[base + 3] += factor * grad[3]
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
            norma_logf(4, :warning, "return_to_surface: node %d did not reach the surface tolerance", node_index)
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

# --- Periodic cell (SolidMechanicsPeriodicBoundaryCondition) -------------------

# Components of the macroscopic deformation gradient by name: `xy` is F̄₁₂,
# in the row of x and the column of y.
const PERIODIC_COMPONENTS = (
    "xx" => CartesianIndex(1, 1),
    "xy" => CartesianIndex(1, 2),
    "xz" => CartesianIndex(1, 3),
    "yx" => CartesianIndex(2, 1),
    "yy" => CartesianIndex(2, 2),
    "yz" => CartesianIndex(2, 3),
    "zx" => CartesianIndex(3, 1),
    "zy" => CartesianIndex(3, 2),
    "zz" => CartesianIndex(3, 3),
)

# Names of F̄ and P̄ as global variables of the Exodus output.
const PERIODIC_GLOBAL_NAMES = vcat(
    ["F_" * name for (name, _) in PERIODIC_COMPONENTS], ["P_" * name for (name, _) in PERIODIC_COMPONENTS]
)

# The nodes of a pair of node sets are matched within this distance, relative
# to the diameter of the mesh.
const PERIODIC_MATCH_TOLERANCE = 1.0e-6

# A component expression compiled to a function of the time alone.
function periodic_component_function(expression, label::String)
    num = Num(eval(Meta.parse(string(expression))))
    if any(variable -> !isequal(variable, Symbolics.unwrap(t)), Symbolics.get_variables(num))
        norma_abort("Periodic boundary condition: `$label` must be a function of t alone, got \"$expression\".")
    end
    fun = eval(build_function(num, [t, x, y, z]; expression=Val(false)))
    return time -> Float64(fun((time, 0.0, 0.0, 0.0)))
end

# The image in `upper` of every node of `lower` under the translation `shift`.
# The nodes of `upper` are binned in cells of the size of the tolerance, so a
# node within the tolerance of a point lies in the cell of the point or in one
# of its 26 neighbors.
function match_periodic_nodes(
    reference::Matrix{Float64},
    lower::Vector{Int64},
    upper::Vector{Int64},
    shift::SVector{3,Float64},
    tolerance::Float64,
    names::Vector{String},
)
    cell(X) = (floor(Int64, X[1] / tolerance), floor(Int64, X[2] / tolerance), floor(Int64, X[3] / tolerance))
    grid = Dict{NTuple{3,Int64},Vector{Int64}}()
    for node in upper
        push!(get!(grid, cell(view(reference, :, node)), Int64[]), node)
    end
    images = similar(lower)
    taken = falses(size(reference, 2))
    for (index, node) in enumerate(lower)
        target = SVector{3,Float64}(view(reference, :, node)) + shift
        center = cell(target)
        image = 0
        distance = tolerance
        for i in -1:1, j in -1:1, k in -1:1
            for candidate in get(grid, (center[1] + i, center[2] + j, center[3] + k), Int64[])
                candidate_distance = norm(SVector{3,Float64}(view(reference, :, candidate)) - target)
                if candidate_distance <= distance
                    image = candidate
                    distance = candidate_distance
                end
            end
        end
        if image == 0
            norma_abortf(
                "Periodic boundary condition: node %d of node set %s, at (%.6e, %.6e, %.6e), has no image in " *
                "node set %s. The mesh must conform across the periodic faces.",
                node,
                names[1],
                reference[1, node],
                reference[2, node],
                reference[3, node],
                names[2],
            )
        end
        if taken[image]
            norma_abortf(
                "Periodic boundary condition: node %d of node set %s is the image of two nodes of node set %s.",
                image,
                names[2],
                names[1],
            )
        end
        taken[image] = true
        images[index] = image
    end
    return images
end

# Nodes linked by the pairs form classes of periodic images. A breadth-first
# walk of each class gives every node its position in periods along each
# pair, relative to the first node of the class, and checks that the
# positions are consistent.
function periodic_image_classes(lowers::Vector{Vector{Int64}}, uppers::Vector{Vector{Int64}})
    num_pairs = length(lowers)
    links = Dict{Int64,Vector{NTuple{3,Int64}}}()
    for k in 1:num_pairs, (a, b) in zip(lowers[k], uppers[k])
        push!(get!(links, a, NTuple{3,Int64}[]), (b, k, 1))
        push!(get!(links, b, NTuple{3,Int64}[]), (a, k, -1))
    end
    position = Dict{Int64,Vector{Int64}}()
    classes = Vector{Vector{Int64}}()
    for start in sort!(collect(keys(links)))
        haskey(position, start) && continue
        position[start] = zeros(Int64, num_pairs)
        members = [start]
        head = 1
        while head <= length(members)
            node = members[head]
            head += 1
            for (other, k, step) in links[node]
                other_position = copy(position[node])
                other_position[k] += step
                if haskey(position, other) == false
                    position[other] = other_position
                    push!(members, other)
                elseif position[other] != other_position
                    norma_abort(
                        "Periodic boundary condition: node $other is reached through the pairs at two different " *
                        "positions of the periodic lattice. The node sets of the pairs are inconsistent.",
                    )
                end
            end
        end
        push!(classes, members)
    end
    return classes, position
end

# The lower of nodes a and b by x, then y, then z, with coordinates within
# `tolerance` equal; a if they coincide.
function lower_node(reference::Matrix{Float64}, a::Int64, b::Int64, tolerance::Float64)
    for j in 1:3
        difference = reference[j, b] - reference[j, a]
        abs(difference) > tolerance && return difference < 0.0 ? b : a
    end
    return a
end

function SolidMechanicsPeriodicBoundaryCondition(
    sim::SingleDomainSimulation, input_mesh::ExodusDatabase, bc_params::Parameters
)
    if sim.parent !== nothing
        norma_abort("Periodic boundary conditions apply to single-domain simulations only.")
    end
    integrator_params = get(sim.params, "time integrator", Parameters())
    integrator_type = get(integrator_params, "type", "")
    solver_type = get(get(sim.params, "solver", Parameters()), "type", "")
    if integrator_type != "quasi static" || solver_type != "Hessian minimizer"
        norma_abort(
            "Periodic boundary conditions require the `quasi static` time integrator and the " *
            "`Hessian minimizer` solver.",
        )
    end
    # The solve condenses the stiffness onto the independent degrees of
    # freedom, which the matrices of `matrix output` would not describe.
    if haskey(sim.params, "matrix output")
        norma_abort(
            "`matrix output` does not support periodic boundary conditions: the solve uses the stiffness " *
            "condensed onto the independent degrees of freedom, which the written matrices do not describe.",
        )
    end
    model = sim.model
    model isa SolidMechanics || norma_abort("Periodic boundary conditions apply to solid mechanics models only.")
    reference = model.reference
    num_nodes = size(reference, 2)
    diameter = norm(maximum(reference; dims=2) - minimum(reference; dims=2))
    tolerance = PERIODIC_MATCH_TOLERANCE * diameter

    pair_names = get(bc_params, "node set pairs", nothing)
    if !(pair_names isa AbstractVector) ||
        isempty(pair_names) ||
        length(pair_names) > 3 ||
        any(pair -> !(pair isa AbstractVector) || length(pair) != 2, pair_names)
        norma_abort(
            "Periodic boundary condition: `node set pairs` must list one to three pairs of node sets, " *
            "such as [[nsx-, nsx+], [nsy-, nsy+], [nsz-, nsz+]].",
        )
    end
    lowers = Vector{Int64}[]
    uppers = Vector{Int64}[]
    pair_axes = Int64[]
    lengths = Float64[]
    for pair in pair_names
        names = String.(pair)
        lower, upper = (
            Int64.(Exodus.read_node_set_nodes(input_mesh, node_set_id_from_name(n, input_mesh))) for n in names
        )
        if length(lower) != length(upper)
            norma_abortf(
                "Periodic boundary condition: node sets %s and %s have %d and %d nodes; the mesh must conform " *
                "across the periodic faces.",
                names[1],
                names[2],
                length(lower),
                length(upper),
            )
        end
        # The translation between the two node sets, from their centroids.
        shift = vec(sum(reference[:, upper]; dims=2) - sum(reference[:, lower]; dims=2)) / length(lower)
        axis = argmax(abs.(shift))
        if any(j -> j != axis && abs(shift[j]) > tolerance, 1:3) || shift[axis] <= tolerance
            norma_abortf(
                "Periodic boundary condition: node set %s must be the image of node set %s under a translation " *
                "along the positive x, y or z axis; the translation is (%.6e, %.6e, %.6e).",
                names[2],
                names[1],
                shift...,
            )
        end
        if axis in pair_axes
            norma_abort("Periodic boundary condition: two pairs of node sets are periodic along the same axis.")
        end
        period = shift[axis]
        exact_shift = SVector{3,Float64}(ntuple(j -> j == axis ? period : 0.0, 3))
        push!(lowers, lower)
        push!(uppers, match_periodic_nodes(reference, lower, upper, exact_shift, tolerance, names))
        push!(pair_axes, axis)
        push!(lengths, period)
    end
    num_pairs = length(pair_axes)

    # The origin is the node common to the first node sets of all pairs, and
    # the corners are its images. With one or two pairs the common nodes form
    # a face or an edge, and the origin is the lowest of them by x, then y,
    # then z, with coordinates within the matching tolerance equal: the common
    # nodes share their coordinates along the pair axes only up to roundoff,
    # which must not pick the origin.
    common = reduce(intersect, lowers)
    if isempty(common)
        norma_abort(
            "Periodic boundary condition: the first node sets of the pairs have no node in common, so the cell " *
            "has no corner to hold fixed.",
        )
    end
    origin = reduce((a, b) -> lower_node(reference, a, b, tolerance), common)
    corners = [uppers[k][findfirst(==(origin), lowers[k])] for k in 1:num_pairs]

    # Every class of images has one independent node, the one of the lowest
    # position along every pair; the others, except the corners, depend on it.
    classes, position = periodic_image_classes(lowers, uppers)
    dependent_nodes = Int64[]
    independent_nodes = Int64[]
    periods = Vector{Int64}[]
    for members in classes
        lowest = reduce((a, b) -> min.(a, b), (position[n] for n in members))
        index = findfirst(n -> position[n] == lowest, members)
        index === nothing && norma_abort(
            "Periodic boundary condition: the images of node $(members[1]) have no lowest image along every pair."
        )
        independent = members[index]
        for node in members
            (node == independent || node in corners) && continue
            push!(dependent_nodes, node)
            push!(independent_nodes, independent)
            push!(periods, position[node] - lowest)
        end
    end
    period_matrix = isempty(periods) ? zeros(Int64, num_pairs, 0) : reduce(hcat, periods)

    # u = A u: the identity on the independent degrees of freedom, and on those
    # of a dependent node the image plus n_k (u(C_k) - u(O)).
    is_dependent = falses(num_nodes)
    is_dependent[dependent_nodes] .= true
    rows = Int64[]
    columns = Int64[]
    values = Float64[]
    for node in 1:num_nodes, i in 1:3
        is_dependent[node] && continue
        push!(rows, 3 * (node - 1) + i)
        push!(columns, 3 * (node - 1) + i)
        push!(values, 1.0)
    end
    for (index, node) in enumerate(dependent_nodes), i in 1:3
        row = 3 * (node - 1) + i
        push!(rows, row)
        push!(columns, 3 * (independent_nodes[index] - 1) + i)
        push!(values, 1.0)
        for k in 1:num_pairs
            n = period_matrix[k, index]
            n == 0 && continue
            append!(rows, (row, row))
            append!(columns, (3 * (corners[k] - 1) + i, 3 * (origin - 1) + i))
            append!(values, (n, -n))
        end
    end
    transformation = dropzeros!(sparse(rows, columns, values, 3 * num_nodes, 3 * num_nodes))

    # Volume of the box: the periods along the periodic axes and the extent of
    # the mesh along the others.
    extent = vec(maximum(reference; dims=2) - minimum(reference; dims=2))
    volume = prod(j in pair_axes ? lengths[findfirst(==(j), pair_axes)] : extent[j] for j in 1:3)

    # Control of the components of F̄ along the periodic axes.
    labels = Dict(name => component for (name, component) in PERIODIC_COMPONENTS)
    deformation = get(bc_params, "deformation gradient", Parameters())
    stress = get(bc_params, "stress", Parameters())
    for (key, entries) in (("deformation gradient", deformation), ("stress", stress))
        entries isa AbstractDict ||
            norma_abort("Periodic boundary condition: `$key` must map components xx, xy, ..., zz to functions of t.")
        for name in keys(entries)
            haskey(labels, name) || norma_abort(
                "Periodic boundary condition: unknown component `$name` of `$key`; " *
                "the components are xx, xy, xz, yx, yy, yz, zx, zy and zz.",
            )
            if labels[name][2] ∉ pair_axes
                norma_abort(
                    "Periodic boundary condition: component `$name` of `$key` is in the column of the " *
                    "$(("x", "y", "z")[labels[name][2]]) axis, which is not periodic.",
                )
            end
            if key == "stress" && haskey(deformation, name)
                norma_abort(
                    "Periodic boundary condition: component `$name` is in both `deformation gradient` and `stress`."
                )
            end
        end
    end
    components = CartesianIndex{2}[]
    component_pairs = Int64[]
    stress_controlled = Bool[]
    component_funs = Function[]
    for (name, component) in PERIODIC_COMPONENTS
        k = findfirst(==(component[2]), pair_axes)
        k === nothing && continue
        push!(components, component)
        push!(component_pairs, k)
        if haskey(stress, name)
            push!(stress_controlled, true)
            push!(component_funs, periodic_component_function(stress[name], "stress: $name"))
        else
            push!(stress_controlled, false)
            identity_value = component[1] == component[2] ? "1.0" : "0.0"
            expression = get(deformation, name, identity_value)
            push!(component_funs, periodic_component_function(expression, "deformation gradient: $name"))
        end
    end
    # With F̄_ij and F̄_ji both free, a rigid rotation in the ij plane costs no
    # energy and the tangent is singular.
    for (index, component) in enumerate(components)
        stress_controlled[index] && component[1] < component[2] || continue
        transpose_index = findfirst(==(CartesianIndex(component[2], component[1])), components)
        if transpose_index !== nothing && stress_controlled[transpose_index]
            names = Dict(c => n for (n, c) in PERIODIC_COMPONENTS)
            norma_abort(
                "Periodic boundary condition: components $(names[component]) and " *
                "$(names[CartesianIndex(component[2], component[1])]) are both under `stress`, which leaves " *
                "the rigid rotation free. Prescribe one of them in `deformation gradient`.",
            )
        end
    end
    norma_logf(
        0,
        :setup,
        "Periodic: %d pair(s), %d dependent nodes, origin node %d, volume %.6e",
        num_pairs,
        length(dependent_nodes),
        origin,
        volume,
    )
    return SolidMechanicsPeriodicBoundaryCondition(
        "periodic",
        pair_axes,
        [String.(pair) for pair in pair_names],
        lengths,
        volume,
        origin,
        corners,
        dependent_nodes,
        independent_nodes,
        period_matrix,
        components,
        component_pairs,
        stress_controlled,
        component_funs,
        zeros(3, 3),
        zeros(3, 3),
        zeros(3),
        zeros(3),
        transformation,
        sparse(transpose(transformation)),
        false,
    )
end

periodic_boundary_condition(model::Model) = nothing

function periodic_boundary_condition(model::SolidMechanics)
    index = findfirst(bc -> bc isa SolidMechanicsPeriodicBoundaryCondition, model.boundary_conditions)
    return index === nothing ? nothing : model.boundary_conditions[index]
end

# F̄ - I from the corners: u(C_k) - u(O) = (F̄ - I) L_k e_k. The columns of the
# axes that are not periodic are zero.
function periodic_displacement_gradient(model::SolidMechanics, bc::SolidMechanicsPeriodicBoundaryCondition)
    H = zeros(3, 3)
    for (k, axis) in enumerate(bc.axes), i in 1:3
        H[i, axis] = (model.displacement[i, bc.corners[k]] - model.displacement[i, bc.origin]) / bc.lengths[k]
    end
    return H
end

# The macroscopic deformation gradient F̄, the volume average of F over the
# cell. Its columns along the periodic axes are those the corners carry, which
# for a periodic field equal that average. With one or two pairs, the other
# columns, which the condition does not control, are averaged over the
# elements.
function periodic_deformation_gradient(model::SolidMechanics)
    bc = periodic_boundary_condition(model)
    bc === nothing && norma_abort("The model has no periodic boundary condition.")
    H = periodic_displacement_gradient(model, bc)
    other_axes = setdiff(1:3, bc.axes)
    if !isempty(other_axes)
        H[:, other_axes] = integrated_displacement_gradient(model)[:, other_axes] / bc.volume
    end
    return I + H
end

# ∫ ∇u dV over the elements, by their quadrature.
function integrated_displacement_gradient(model::SolidMechanics)
    H = zeros(3, 3)
    for block in model.blocks
        # Function barrier over the static shape function tables of the block.
        H += integrated_displacement_gradient(model, block.N, block.dN, block.weights, block.connectivity)
    end
    return H
end

function integrated_displacement_gradient(
    model::SolidMechanics, N::SMatrix, dN::SArray, weights::AbstractVector, connectivity::Matrix{Int64}
)
    H = zeros(SMatrix{3,3,Float64,9})
    for element in axes(connectivity, 2)
        nodes = view(connectivity, :, element)
        X = gather_nodal(model.reference, nodes, N)
        u = gather_nodal(model.displacement, nodes, N)
        for point in eachindex(weights)
            dNdξ = dN[:, :, point]
            dXdξ = dNdξ * X'
            H += u * (dXdξ \ dNdξ)' * (det(dXdξ) * weights[point])
        end
    end
    return Matrix(H)
end

# The volume average of the first Piola-Kirchhoff stress over the cell,
# (1/V) Σ_a f_a ⊗ X_a over the internal forces of all nodes. The shape
# functions interpolate X exactly, so the sum is ∫ P dV of the finite element
# field. The forces of the Robin and penalty Surface conditions, which
# evaluate adds to the internal force, are not those of the material and are
# left out.
function periodic_stress(model::SolidMechanics)
    bc = periodic_boundary_condition(model)
    bc === nothing && norma_abort("The model has no periodic boundary condition.")
    P = zeros(3, 3)
    springs = zeros(length(model.internal_force))
    apply_robin_bcs_internal_force!(model, springs)
    apply_surface_penalty_internal_force!(model, springs)
    force = model.internal_force - springs
    for node in axes(model.reference, 2), j in 1:3, i in 1:3
        P[i, j] += force[3 * (node - 1) + i] * model.reference[j, node]
    end
    return P / bc.volume
end

# Prescribed conditions may hold the origin and the nodes of the first node
# sets that are on no second node set, but not the nodes that the periodic
# condition moves: the dependent nodes and the corners. With two or three
# pairs, a first node set also holds its edges on the second node sets of the
# other pairs, so a condition on a whole first node set is rejected too.
function check_periodic_conflicts!(model::SolidMechanics, bc::SolidMechanicsPeriodicBoundaryCondition)
    moved = falses(size(model.reference, 2))
    moved[bc.dependent_nodes] .= true
    moved[bc.corners] .= true
    for other in model.boundary_conditions
        nodes = if other isa SolidMechanicsDirichletBoundaryCondition
            other.node_set_node_indices
        elseif other isa SolidMechanicsSideSetDirichletBoundaryCondition
            other.node_indices
        else
            continue
        end
        held = filter(node -> moved[node], nodes)
        isempty(held) && continue
        norma_abort(
            "The Dirichlet boundary condition on \"$(other.name)\" prescribes " *
            "$(join(periodic_nodes_held(bc, held), ", ", " and ")), which the periodic boundary condition " *
            "determines: a node of the second node set of a pair follows its image on the first, and a corner " *
            "carries a column of F̄. A Dirichlet condition may hold the origin and the nodes of the first node " *
            "sets that are on no second node set; with two or three pairs, a first node set also holds its " *
            "edges on the second node sets of the other pairs.",
        )
    end
    check_periodic_rotations(model, bc)
    bc.checked = true
    return nothing
end

# What the nodes `held`, which the periodic condition moves, are: nodes of the
# second node set of a pair, or corners.
function periodic_nodes_held(bc::SolidMechanicsPeriodicBoundaryCondition, held::Vector{Int64})
    held_set = Set(held)
    descriptions = String[]
    for k in eachindex(bc.axes)
        on_second = any(
            index -> bc.periods[k, index] > 0 && bc.dependent_nodes[index] in held_set, eachindex(bc.dependent_nodes)
        )
        on_second && push!(descriptions, "nodes of $(bc.pair_names[k][2])")
    end
    for (k, corner) in enumerate(bc.corners)
        corner in held_set && push!(descriptions, "the corner across [$(join(bc.pair_names[k], ", "))]")
    end
    return descriptions
end

# The rigid rotation in the plane of axes a and b, u = W (X - X_O) with
# W_ba = -W_ab = 1, is periodic with F̄ - I = W and leaves the origin in place.
# A prescribed component of F̄ in that plane holds it. The constructor rejects
# both components of a plane under `stress`. With one or two pairs, a plane
# with an axis that is not periodic has one component of F̄, or none, and if
# none is prescribed another condition has to hold the rotation. Otherwise the
# rotation costs no energy, the tangent is singular, and F̄ ends with an
# arbitrary rotation.
function check_periodic_rotations(model::SolidMechanics, bc::SolidMechanicsPeriodicBoundaryCondition)
    X = model.reference .- model.reference[:, bc.origin]
    tolerance = PERIODIC_MATCH_TOLERANCE * maximum(abs, X)
    axis_names = ("x", "y", "z")
    component_names = Dict(c => n for (n, c) in PERIODIC_COMPONENTS)
    for (a, b) in ((1, 2), (1, 3), (2, 3))
        in_plane = findall(c -> c == CartesianIndex(a, b) || c == CartesianIndex(b, a), bc.components)
        length(in_plane) == 2 && continue
        any(index -> !bc.stress_controlled[index], in_plane) && continue
        rotation = zeros(size(X))
        rotation[a, :] = -X[b, :]
        rotation[b, :] = X[a, :]
        holds_periodic_mode(model, rotation, tolerance) && continue
        # An axis of the plane that is not periodic: a node away from the
        # origin along it can hold the rotation.
        away, other = a ∉ bc.axes ? (a, b) : (b, a)
        plane = axis_names[a] * "–" * axis_names[b]
        reason, alternative = if isempty(in_plane)
            "neither $(axis_names[a]) nor $(axis_names[b]) is periodic, so no component of F̄ holds it", ""
        else
            name = component_names[bc.components[only(in_plane)]]
            "$(axis_names[away]) is not periodic, and $name, the one component of F̄ in that plane, is under " *
            "`stress`",
            ", or prescribe $name in `deformation gradient`"
        end
        norma_abort(
            "Periodic boundary condition: the rigid rotation in the $plane plane is free: $reason. The tangent is " *
            "then singular. Hold the rotation with a Dirichlet condition on component $(axis_names[other]) of a " *
            "node away from the origin along $(axis_names[away])$alternative.",
        )
    end
    return nothing
end

# Whether a condition other than the periodic one holds the displacement field
# `mode` (3 × nodes): a Dirichlet or Robin condition on a component the field
# moves at one of the condition's nodes, or a Surface condition whose normal
# the field moves along.
function holds_periodic_mode(model::SolidMechanics, mode::Matrix{Float64}, tolerance::Float64)
    moves(component, nodes) = any(node -> abs(mode[component, node]) > tolerance, nodes)
    for other in model.boundary_conditions
        if other isa SolidMechanicsDirichletBoundaryCondition
            moves(other.offset, other.node_set_node_indices) && return true
        elseif other isa SolidMechanicsSideSetDirichletBoundaryCondition
            moves(other.offset, other.node_indices) && return true
        elseif other isa SolidMechanicsRobinBoundaryCondition
            moves(other.offset, other.side_set_node_indices) && return true
        elseif other isa SolidMechanicsSurfaceBoundaryCondition
            for node in other.node_indices
                x = model.reference[:, node] + model.displacement[:, node]
                normal = other.level_set_grad([model.time, x[1], x[2], x[3]])
                abs(dot(normal, mode[:, node])) > tolerance * norm(normal) && return true
            end
        end
    end
    return false
end

# The dependent nodes from the independent ones, u = A u. Called before every
# assembly, line search trials included, so that every evaluation sees a
# periodic field. It moves no other degree of freedom, so it leaves a state
# whose dependent nodes already follow as it is, such as the start of a step
# that the tangent predictor evaluates.
function update_periodic_dependents!(model::SolidMechanics)
    for bc in model.boundary_conditions
        bc isa SolidMechanicsPeriodicBoundaryCondition || continue
        displacement = vec(model.displacement)
        displacement .= bc.transformation * displacement
    end
    return nothing
end

# The corners with prescribed components of F̄, from the origin, and then the
# dependent nodes. Called once all the conditions are applied, so that the
# corners follow the origin also where a Dirichlet condition prescribes it.
# The translation of the origin is recorded for the prediction of the solve,
# with the increment of F̄.
function finalize_periodic_bcs!(model::SolidMechanics)
    for bc in model.boundary_conditions
        bc isa SolidMechanicsPeriodicBoundaryCondition || continue
        u = model.displacement
        bc.pending_translation += u[:, bc.origin] - bc.origin_displacement
        for (index, component) in enumerate(bc.components)
            bc.stress_controlled[index] && continue
            k = bc.component_pairs[index]
            i = component[1]
            u[i, bc.corners[k]] = u[i, bc.origin] + bc.displacement_gradient[component] * bc.lengths[k]
        end
    end
    update_periodic_dependents!(model)
    return nothing
end

# The translation of the origin and the increment of F̄ since the last solve
# applied to the whole cell, Δu_O + ΔF̄ (X - X_O) on the free degrees of
# freedom, or nothing if neither has changed; both are then cleared. It is the
# affine part of the solution: from it, Newton starts from the fluctuation of
# the previous step, without the increment concentrating on the elements next
# to the faces, or next to the origin and its images when a Dirichlet
# condition moves the origin. It is the prediction of the quasi-static solve
# when the predictor of the time integrator has none (see predict). The
# prescribed degrees of freedom already hold their values, and the dependent
# ones follow the independent ones at the next evaluation. The increments
# accumulate over applications without a solve, so that the solve applies all
# of them.
function periodic_affine_increment!(model::SolidMechanics)
    bc = periodic_boundary_condition(model)
    bc === nothing && return nothing
    increment = bc.pending_increment
    translation = bc.pending_translation
    bc.pending_increment = zeros(3, 3)
    bc.pending_translation = zeros(3)
    any(!iszero, increment) || any(!iszero, translation) || return nothing
    X = model.reference
    return vec(translation .+ increment * (X .- X[:, bc.origin])) .* model.free_dofs
end

# A model that takes over a converged state (transfer_state!) predicts
# nothing from its earlier applications, which saw another state: that of a
# swap's new model, at zero displacement before the state was copied in. The
# next application measures the increments from the converged state.
function clear_periodic_increments!(model::Model)
    bc = periodic_boundary_condition(model)
    bc === nothing && return nothing
    bc.pending_increment = zeros(3, 3)
    bc.pending_translation = zeros(3)
    return nothing
end

# The residual and the tangent on the independent degrees of freedom: Aᵀ r and
# Aᵀ K A. The rows of the dependent ones are zero, and they are not free.
function condense_periodic_bcs!(model::SolidMechanics, solver::HessianMinimizer)
    for bc in model.boundary_conditions
        bc isa SolidMechanicsPeriodicBoundaryCondition || continue
        solver.gradient = bc.transformation_transpose * solver.gradient
        solver.hessian = bc.transformation_transpose * (solver.hessian * bc.transformation)
    end
    return nothing
end

# Applies F̄(t) and the stress on its free components, and holds the origin,
# the prescribed components of the corners and the dependent nodes;
# finalize_periodic_bcs! sets their displacements once all the conditions are
# applied. The origin keeps its displacement, zero unless a Dirichlet
# condition prescribes it. The free components keep their previous values.
# The increment of F̄ is recorded for the prediction of the solve, which
# applies it to the whole cell (periodic_affine_increment!). The increments
# are measured from the state before the application: the Dirichlet
# conditions, which may move the origin, are applied after this one
# (bc_type_rank), and finalize_periodic_bcs! records the translation.
function apply_bc(model::SolidMechanics, bc::SolidMechanicsPeriodicBoundaryCondition)
    bc.checked || check_periodic_conflicts!(model, bc)
    previous = periodic_displacement_gradient(model, bc)
    bc.origin_displacement = model.displacement[:, bc.origin]
    H = copy(previous)
    for (index, component) in enumerate(bc.components)
        bc.stress_controlled[index] && continue
        H[component] = bc.component_funs[index](model.time) - (component[1] == component[2] ? 1.0 : 0.0)
    end
    bc.pending_increment += H - previous
    bc.displacement_gradient = H
    for i in 1:3
        model.free_dofs[3 * (bc.origin - 1) + i] = false
    end
    for (index, component) in enumerate(bc.components)
        k = bc.component_pairs[index]
        dof = 3 * (bc.corners[k] - 1) + component[1]
        if bc.stress_controlled[index]
            model.boundary_force[dof] += bc.volume / bc.lengths[k] * bc.component_funs[index](model.time)
        else
            model.free_dofs[dof] = false
        end
    end
    for node in bc.dependent_nodes, i in 1:3
        model.free_dofs[3 * (node - 1) + i] = false
    end
    return nothing
end

function collect_overlap_candidate_nodes(
    coupled_block_id::Int64, tol::Float64, coupled_subsim::Simulation, subsim::Simulation
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
                    _, _, node_found = find_point_in_mesh(
                        dst_model.reference[:, node_index], coupled_subsim.model, coupled_block_id, tol
                    )
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
    coupled_block_id::Int64, element_type::ElementType, tol::Float64, coupled_subsim::Simulation, subsim::Simulation
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
        error(
            "Invalid value for compute_overlap_l2_error: \"$(bc.compute_overlap_l2_error)\". Valid options are \"disp\", \"velo\", and \"acce\".",
        )
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
