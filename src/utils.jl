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

@inline _use_color() =
    get(ENV, "NORMA_NO_COLOR", "false") != "true" &&
    (get(ENV, "FORCE_COLOR", "") != "" ||
     (isdefined(Base, :have_color) && Base.have_color === true) ||
     stdout isa Base.TTY)

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
    :progress => :green,
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

const NORMA_WRITE_LOG_FILE = Ref(true)
const NORMA_LOG_FILE = Ref{Union{IOStream,Nothing}}(nothing)

function open_log_file(input_file::AbstractString)
    NORMA_WRITE_LOG_FILE[] || return nothing
    NORMA_LOG_FILE[] === nothing || return nothing  # outermost run() owns the file
    path = first(splitext(input_file)) * ".log"
    NORMA_LOG_FILE[] = open(path, "w")
    return nothing
end

function close_log_file()
    io = NORMA_LOG_FILE[]
    io === nothing && return nothing
    close(io)
    NORMA_LOG_FILE[] = nothing
    return nothing
end

function norma_log(level::Int, keyword::Symbol, msg::AbstractString)
    indent = " "^level
    keyword_str = uppercase(string(keyword))[1:min(end, 7)]
    bracketed = "[" * keyword_str * "]"
    padded = rpad(bracketed, 9)
    prefix = indent * padded * " "

    if _use_color()
        color = get(NORMA_COLORS, keyword, :default)
        printstyled(prefix; color=color, bold=true)
    else
        print(prefix)
    end
    if visible_length(prefix * msg) <= 120
        println(msg)
    else
        wrapped = wrap_lines(msg, prefix; width=120 - length(prefix))
        println(wrapped)
    end

    io = NORMA_LOG_FILE[]
    if io !== nothing
        plain = replace(msg, r"\e\[[0-9;]*m" => "")
        if visible_length(prefix * plain) <= 120
            println(io, prefix, plain)
        else
            println(io, prefix, wrap_lines(plain, prefix; width=120 - length(prefix)))
        end
        flush(io)
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
    push!(time_str, @sprintf("%.2fs", seconds))

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
        return _use_color() ? "\e[33m[WAIT]\e[39m" : "[WAIT]"  # yellow
    elseif status == "[DONE]"
        return _use_color() ? "\e[32m[DONE]\e[39m" : "[DONE]"  # green
    elseif status == "[CONVERGING]"
        return _use_color() ? "\e[33m[CONVERGING]\e[39m" : "[CONVERGING]"  # yellow
    elseif status == "[CONVERGED]"
        return _use_color() ? "\e[32m[CONVERGED]\e[39m" : "[CONVERGED]"  # green
    else
        return status  # fallback (no color)
    end
end

function stripped_name(file::AbstractString)
    return first(splitext(basename(file)))
end
