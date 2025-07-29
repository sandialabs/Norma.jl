using Pkg

# Ensure required packages are installed
function ensure_installed(pkgs)
    for pkg in pkgs
        if !haskey(Pkg.project().dependencies, pkg)
            println("Installing missing package: $pkg")
            Pkg.add(pkg)
        end
    end
end

ensure_installed(["DelimitedFiles", "LinearAlgebra", "Printf", "PyPlot"])

# Now safely use the packages
using DelimitedFiles   # For readdlm
using LinearAlgebra    # For norm
using Printf           # For @sprintf
using PyPlot           # For Matlab-like plotting

function plot_solns_overlap()
    # -- Reading in initial reference coordinates --
    coords1 = readdlm("01-refe.csv", ',')
    N1 = length(coords1)
    n1 = N1 ÷ 3
    x1 = coords1[:, 1]
    y1 = coords1[:, 2]
    z1 = coords1[:, 3]
    ind1 = findall((x1 .== minimum(x1)) .& (y1 .== minimum(y1)))

    coords2 = readdlm("02-refe.csv", ',')
    N2 = length(coords2)
    n2 = N2 ÷ 3
    x2 = coords2[:, 1]
    y2 = coords2[:, 2]
    z2 = coords2[:, 3]
    ind2 = findall((x2 .== minimum(x2)) .& (y2 .== minimum(y2)))

    # Collect unique z-values across both files
    z = sort(union(z1, z2))

    # Initialize storage for displacements, velocities, accelerations
    dispz1 = zeros(n1, 0)
    veloz1 = zeros(n1, 0)
    accez1 = zeros(n1, 0)

    dispz2 = zeros(n2, 0)
    veloz2 = zeros(n2, 0)
    accez2 = zeros(n2, 0)

    # Will hold merged (computed) solutions and exact solutions
    disp_computed = zeros(0, 0)
    velo_computed = zeros(0, 0)
    acce_computed = zeros(0, 0)

    disp_exact = zeros(0, 0)
    velo_exact = zeros(0, 0)
    acce_exact = zeros(0, 0)

    # Plotting parameters
    save_figs = 1
    ctr = 1

    # List all files in the directory
    all_files = readdir()

    # Define regex pattern for the filenames
    pattern = r"^01-time-\d{4}\.csv$"

    # Filter files that match the pattern
    matching_files = filter(f -> occursin(pattern, f), all_files)

    # Count them
    num_files = length(matching_files)

    # -- Main loop --
    for i in 0:(num_files - 1)
        # Read data for subdomain #1
        d1 = readdlm(@sprintf("01-disp-%04d.csv", i), ',')
        v1 = readdlm(@sprintf("01-velo-%04d.csv", i), ',')
        a1 = readdlm(@sprintf("01-acce-%04d.csv", i), ',')
        t1 = readdlm(@sprintf("01-time-%04d.csv", i), ',')
        t = t1[1, 1]

        dispz1 = hcat(dispz1, d1[:, 3])
        veloz1 = hcat(veloz1, v1[:, 3])
        accez1 = hcat(accez1, a1[:, 3])

        # Read data for subdomain #2
        d2 = readdlm(@sprintf("02-disp-%04d.csv", i), ',')
        v2 = readdlm(@sprintf("02-velo-%04d.csv", i), ',')
        a2 = readdlm(@sprintf("02-acce-%04d.csv", i), ',')

        dispz2 = hcat(dispz2, d2[:, 3])
        veloz2 = hcat(veloz2, v2[:, 3])
        accez2 = hcat(accez2, a2[:, 3])

        # Grab z-values at ind1, ind2 and merge
        z1ind1 = z1[ind1]
        z2ind2 = z2[ind2]
        zz = unique(sort(vcat(z1ind1, z2ind2)))

        # The current column from dispz1, etc. for these specific indices
        dz1ind1 = dispz1[ind1, ctr]
        dz2ind2 = dispz2[ind2, ctr]
        vz1ind1 = veloz1[ind1, ctr]
        vz2ind2 = veloz2[ind2, ctr]
        az1ind1 = accez1[ind1, ctr]
        az2ind2 = accez2[ind2, ctr]

        # Merge subdomain #1 and #2 solutions at overlapping z
        dispz_merged = zeros(length(zz))
        veloz_merged = zeros(length(zz))
        accez_merged = zeros(length(zz))

        for j in eachindex(zz)
            # find matching points in each subdomain
            ii1 = findall(x -> x == zz[j], z1ind1)
            ii2 = findall(x -> x == zz[j], z2ind2)

            val_disp = 0.0
            val_velo = 0.0
            val_acce = 0.0
            count = 0

            if !isempty(ii1)
                val_disp += dz1ind1[ii1[1]]
                val_velo += vz1ind1[ii1[1]]
                val_acce += az1ind1[ii1[1]]
                count += 1
            end
            if !isempty(ii2)
                val_disp += dz2ind2[ii2[1]]
                val_velo += vz2ind2[ii2[1]]
                val_acce += az2ind2[ii2[1]]
                count += 1
            end
            if count > 1
                # Average if both subdomains had a point at this z
                val_disp /= 2
                val_velo /= 2
                val_acce /= 2
            end

            dispz_merged[j] = val_disp
            veloz_merged[j] = val_velo
            accez_merged[j] = val_acce
        end

        # Exact solution (Gaussian pulses)
        c = 1000.0
        a = 0.001
        scale = a
        b = 0.0
        s = 0.02
        T = 0.001

        # Displacement
        d_ex =
            0.5 .* a .* (exp.(-((zz .- c .* t .- b) .^ 2) ./ (2s^2)) .+ exp.(-((zz .+ c .* t .- b) .^ 2) ./ (2s^2))) .-
            0.5 .* a .*
            (exp.(-((zz .- c .* (T .- t) .- b) .^ 2) ./ (2s^2)) .+ exp.(-((zz .+ c .* (T .- t) .- b) .^ 2) ./ (2s^2)))

        # Velocity
        v_ex =
            (c * 0.5 * a / s^2) .* ((
                (zz .- c .* t .- b) .* exp.(-((zz .- c .* t .- b) .^ 2) / (2s^2)) .-
                (zz .+ c .* t .- b) .* exp.(-((zz .+ c .* t .- b) .^ 2) / (2s^2))
            )) .+
            (c * 0.5 * a / s^2) .* ((
                (zz .- c .* (T .- t) .- b) .* exp.(-((zz .- c .* (T .- t) .- b) .^ 2) / (2s^2)) .-
                (zz .+ c .* (T .- t) .- b) .* exp.(-((zz .+ c .* (T .- t) .- b) .^ 2) / (2s^2))
            ))

        # Acceleration
        a_ex =
            0.5 .* a .* (
                -c^2 / s^2 .* exp.(-0.5 .* ((zz .- c .* t .- b) .^ 2) ./ (s^2)) .+
                c^2 / s^4 .* ((zz .- c .* t .- b) .^ 2) .* exp.(-0.5 .* ((zz .- c .* t .- b) .^ 2) ./ (s^2)) .-
                c^2 / s^2 .* exp.(-0.5 .* ((zz .+ c .* t .- b) .^ 2) ./ (s^2)) .+
                c^2 / s^4 .* ((zz .+ c .* t .- b) .^ 2) .* exp.(-0.5 .* ((zz .+ c .* t .- b) .^ 2) ./ (s^2))
            ) .-
            0.5 .* a .* (
                -c^2 / s^2 .* exp.(-0.5 .* ((zz .- c .* (T .- t) .- b) .^ 2) ./ (s^2)) .+
                c^2 / s^4 .* ((zz .- c .* (T .- t) .- b) .^ 2) .*
                exp.(-0.5 .* ((zz .- c .* (T .- t) .- b) .^ 2) ./ (s^2)) .-
                c^2 / s^2 .* exp.(-0.5 .* ((zz .+ c .* (T .- t) .- b) .^ 2) ./ (s^2)) .+
                c^2 / s^4 .* ((zz .+ c .* (T .- t) .- b) .^ 2) .*
                exp.(-0.5 .* ((zz .+ c .* (T .- t) .- b) .^ 2) ./ (s^2))
            )

        # Initialize the storage matrices on the first iteration
        if i == 0
            disp_computed = zeros(length(zz), 0)
            velo_computed = zeros(length(zz), 0)
            acce_computed = zeros(length(zz), 0)
            disp_exact = zeros(length(zz), 0)
            velo_exact = zeros(length(zz), 0)
            acce_exact = zeros(length(zz), 0)
        end

        # Store the merged numeric solutions and the exact ones as new columns
        disp_computed = hcat(disp_computed, dispz_merged)
        velo_computed = hcat(velo_computed, veloz_merged)
        acce_computed = hcat(acce_computed, accez_merged)

        disp_exact = hcat(disp_exact, d_ex)
        velo_exact = hcat(velo_exact, v_ex)
        acce_exact = hcat(acce_exact, a_ex)

        # -- Plotting --
        clf()  # clear the figure for each loop iteration

        # subplot 1: displacement
        subplot(3, 1, 1)
        plot(z1[ind1], dispz1[ind1, ctr], "-b"; label="subdomain 1")
        plot(z2[ind2], dispz2[ind2, ctr], "-r"; label="subdomain 2")
        plot(zz, d_ex, "--c"; label="exact")
        xlabel("z")
        ylabel("z-disp")
        title("displacement snapshot $(i+1) at time = $t")
        axis([minimum(z), maximum(z), -scale, scale])
        legend(; loc="lower right")

        # subplot 2: velocity
        subplot(3, 1, 2)
        plot(z1[ind1], veloz1[ind1, ctr], "-b"; label="subdomain 1")
        plot(z2[ind2], veloz2[ind2, ctr], "-r"; label="subdomain 2")
        plot(zz, v_ex, "--c"; label="exact")
        xlabel("z")
        ylabel("z-velo")
        title("velocity snapshot $(i+1) at time = $t")
        axis([minimum(z), maximum(z), -3e4 * scale, 3e4 * scale])
        legend(; loc="lower right")

        # subplot 3: acceleration
        subplot(3, 1, 3)
        plot(z1[ind1], accez1[ind1, ctr], "-b"; label="subdomain 1")
        plot(z2[ind2], accez2[ind2, ctr], "-r"; label="subdomain 2")
        plot(zz, a_ex, "--c"; label="exact")
        xlabel("z")
        ylabel("z-acce")
        title("acceleration snapshot $(i+1) at time = $t")
        axis([minimum(z), maximum(z), -2.5e9 * scale, 2.5e9 * scale])
        legend(; loc="lower right")

        PyPlot.pause(0.001)  # short pause to visualize updates

        # optional figure saving
        if save_figs == 1
            if ctr < 10
                filename = @sprintf("soln_000%d.png", ctr)
            elseif ctr < 100
                filename = @sprintf("soln_00%d.png", ctr)
            elseif ctr < 1000
                filename = @sprintf("soln_0%d.png", ctr)
            else
                filename = @sprintf("soln_%d.png", ctr)
            end
            savefig(filename)
        end
        ctr += 1
    end

    # Compute errors
    sz = size(disp_exact)
    numerator_disp = 0.0
    denominator_disp = 0.0
    numerator_velo = 0.0
    denominator_velo = 0.0
    numerator_acce = 0.0
    denominator_acce = 0.0

    for i in 1:sz[2]
        numerator_disp += norm(disp_computed[:, i] .- disp_exact[:, i])^2
        denominator_disp += norm(disp_exact[:, i])^2

        numerator_velo += norm(velo_computed[:, i] .- velo_exact[:, i])^2
        denominator_velo += norm(velo_exact[:, i])^2

        numerator_acce += norm(acce_computed[:, i] .- acce_exact[:, i])^2
        denominator_acce += norm(acce_exact[:, i])^2
    end

    dispz_relerr = sqrt(numerator_disp / denominator_disp)
    veloz_relerr = sqrt(numerator_velo / denominator_velo)
    accez_relerr = sqrt(numerator_acce / denominator_acce)

    println("z-disp avg rel error = $dispz_relerr")
    println("z-velo avg rel error = $veloz_relerr")
    return println("z-acce avg rel error = $accez_relerr")
end

# Call the function
plot_solns_overlap()
