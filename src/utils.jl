# Norma: Copyright 2025 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software. This software
# is released under the BSD license detailed in the file license.txt in the
# top-level Norma.jl directory.
using Logging
using Printf

function wrap_lines(msg::AbstractString, prefix::AbstractString; width::Int=80)
    words = split(msg)
    lines = String[]
    current = ""

    for word in words
        if length(current) + length(word) + 1 > width
            push!(lines, current)
            current = word
        else
            current = isempty(current) ? word : "$current $word"
        end
    end
    push!(lines, current)

    return join([i == 1 ? line : " "^length(prefix) * line for (i, line) in enumerate(lines)], "\n")
end

const NORMA_COLOR_OUTPUT =
    (isdefined(Base, :have_color) && Base.have_color !== nothing ? Base.have_color : stdout isa Base.TTY) &&
    get(ENV, "NORMA_NO_COLOR", "false") != "true"

const NORMA_COLORS = Dict(
    :abort => :light_red,
    :acceleration => :blue,
    :advance => :green,
    :contact => :light_yellow,
    :debug => :red,
    :domain => :light_blue,
    :done => :green,
    :equilibrium => :blue,
    :error => :red,
    :info => :cyan,
    :initial => :blue,
    :input => :cyan,
    :linesearch => :cyan,
    :norma => :magenta,
    :output => :cyan,
    :recover => :yellow,
    :schwarz => :light_blue,
    :setup => :magenta,
    :solve => :cyan,
    :step => :green,
    :stop => :blue,
    :summary => :magenta,
    :test => :light_green,
    :time => :light_cyan,
    :warning => :yellow,
)

function visible_length(s::AbstractString)
    return length(replace(s, r"\e\[[0-9;]*m" => ""))
end

function norma_log(level::Int, keyword::Symbol, msg::AbstractString)
    indent = " "^level
    keyword_str = uppercase(string(keyword))[1:min(end, 7)]
    bracketed = "[" * keyword_str * "]"
    padded = rpad(bracketed, 9)
    prefix = indent * padded * " "

    if NORMA_COLOR_OUTPUT
        color = get(NORMA_COLORS, keyword, :default)
        printstyled(prefix; color=color, bold=true)
    else
        print(prefix)
    end
    if visible_length(prefix * msg) <= 80
        println(msg)
    else
        wrapped = wrap_lines(msg, prefix; width=80 - length(prefix))
        println(wrapped)
    end
end

function norma_logf(level::Int, keyword::Symbol, fmt::AbstractString, args...)
    fstr = Printf.Format(fmt)
    msg = Printf.format(fstr, args...)
    norma_log(level, keyword, msg)
    return nothing
end

# Internal, testable logic
function _norma_abort_message(msg::AbstractString)
    norma_log(0, :abort, msg)
    return norma_log(0, :norma, "SIMULATION ABORTED")
end

function _norma_abort_messagef(fmt::AbstractString, args...)
    norma_logf(0, :abort, fmt, args...)
    return norma_log(0, :norma, "SIMULATION ABORTED")
end

const NORMA_TEST_MODE = Ref(false)  # top-level constant

struct NormaAbortException <: Exception
    msg::String
end

function norma_abort(msg::AbstractString)
    _norma_abort_message(msg)
    if NORMA_TEST_MODE[]
        throw(NormaAbortException(msg))
    else
        error("Norma aborted: $msg")
    end
end

function norma_abortf(fmt::AbstractString, args...)
    _norma_abort_messagef(fmt, args...)
    f = Printf.Format(fmt)
    msg = Printf.format(f, args...)
    if NORMA_TEST_MODE[]
        throw(NormaAbortException(msg))
    else
        error("Norma aborted: $msg")
    end
end

function log_matrix(level::Int, keyword::Symbol, label::AbstractString, M::AbstractMatrix)
    norma_log(level, keyword, "$label:")
    for i in axes(M, 1)
        row = M[i, :]
        row_str = join(map(val -> @sprintf("%.3e", val), row), " ")
        norma_log(level, keyword, row_str)
    end
end

# Enable debugging for a specific module in the environment variable JULIA_DEBUG (all for all modules)
function configure_logger()
    debug_level = get(ENV, "JULIA_DEBUG", "")

    if debug_level != ""
        # Enable debug logging
        global_logger(ConsoleLogger(stderr, Logging.Debug))
        norma_log(0, :info, "Debugging enabled for: $debug_level")
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
        norma_log(0, :warning, "Enabling FPE traps can break Julia's parser if done too early")
        norma_log(0, :warning, "(e.g., before Meta.parse()).")
        mask = FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW
        ccall((:feenableexcept, "libm.so.6"), Cuint, (Cuint,), mask)
        norma_log(0, :info, "Floating-point exceptions enabled (invalid, div-by-zero, overflow)")
    else
        norma_log(0, :warning, "FPE trap support not available on this platform: $(Sys.KERNEL) $(Sys.ARCH)")
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
        norma_abort("Usage: julia Norma.jl <input_file>")
    end
    return ARGS[1]
end

function colored_status(status::String)
    if status == "[WAIT]"
        return NORMA_COLOR_OUTPUT ? "\e[33m[WAIT]\e[39m" : "[WAIT]"  # yellow
    elseif status == "[DONE]"
        return NORMA_COLOR_OUTPUT ? "\e[32m[DONE]\e[39m" : "[DONE]"  # green
    elseif status == "[CONVERGING]"
        return NORMA_COLOR_OUTPUT ? "\e[33m[CONVERGING]\e[39m" : "[CONVERGING]"  # yellow
    elseif status == "[CONVERGED]"
        return NORMA_COLOR_OUTPUT ? "\e[32m[CONVERGED]\e[39m" : "[CONVERGED]"  # green
    else
        return status  # fallback (no color)
    end
end

function stripped_name(file::AbstractString)
    return first(splitext(basename(file)))
end
