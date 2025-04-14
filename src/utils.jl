# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.
using Logging
using Printf
using Unicode

# Enable debugging for a specific module in the environment variable JULIA_DEBUG (all for all modules)
function configure_logger()
    debug_level = get(ENV, "JULIA_DEBUG", "")

    if debug_level != ""
        # Enable debug logging
        global_logger(ConsoleLogger(stderr, Logging.Debug))
        @info "Debugging enabled for: $debug_level"
    else
        # Default to Info level
        global_logger(ConsoleLogger(stderr, Logging.Info))
    end
    return nothing
end

# Constants for Linux/glibc FPE trap masks
const FE_INVALID = 1    # 0x01
const FE_DIVBYZERO = 4  # 0x04
const FE_OVERFLOW = 8   # 0x08

"""
    enable_fpe_traps()

Attempts to enable floating-point exceptions (FPE traps) on supported platforms.
Catches invalid operations, divide-by-zero, and overflow. Only supported on Linux x86_64.
Warns that parser behavior may break due to this setting.
"""
function enable_fpe_traps()
    if Sys.islinux() && Sys.ARCH == :x86_64
        @warn "Enabling FPE traps can break Julia's parser if done too early (e.g., before Meta.parse())."
        mask = FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW
        ccall((:feenableexcept, "libm.so.6"), Cuint, (Cuint,), mask)
        @info "Floating-point exceptions enabled (invalid, div-by-zero, overflow)"
    else
        @warn "FPE trap support not available on this platform: $(Sys.KERNEL) $(Sys.ARCH)"
    end
    return nothing
end

function format_time(seconds::Float64)::String
    days = floor(Int, seconds / 86400)  # 86400 seconds in a day
    seconds %= 86400

    hours = floor(Int, seconds / 3600)  # 3600 seconds in an hour
    seconds %= 3600

    minutes = floor(Int, seconds / 60)  # 60 seconds in a minute
    seconds %= 60

    # Build the formatted string based on non-zero values
    time_str = []
    if days > 0
        push!(time_str, @sprintf("%dd", days))
    end
    if hours > 0 || days > 0  # Show hours if days > 0 for completeness
        push!(time_str, @sprintf("%dh", hours))
    end
    if minutes > 0 || hours > 0 || days > 0  # Show minutes if higher units are non-zero
        push!(time_str, @sprintf("%dm", minutes))
    end
    push!(time_str, @sprintf("%.1fs", seconds))

    return join(time_str, " ")
end

function parse_args()
    if length(ARGS) != 1
        error("Usage: julia Norma.jl <input_file>")
    end
    return ARGS[1]
end
