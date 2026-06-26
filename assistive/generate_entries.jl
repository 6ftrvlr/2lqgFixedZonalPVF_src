# generate_entries.jl
using CSV

# =========================
# 这里填你的不变代码段
# =========================

const CODE1 = raw"""
#=
Entrance script for zonally-averaged PV-flux 2LQG problems.

Two modes:
1) mode = :fixedU
   Input:  U
   Output: F1, F2

2) mode = :fixedF1
   Input:  F1
   Output: U, F2
=#

########################################
#=======FILE SYSTEM====================#
########################################
"""
#=    TARGET_DIR = "./ZPVF_F7500_K30" =#
const CODE2 = raw"""

if isdir(TARGET_DIR)
    println("Directory already exists: $TARGET_DIR")
else
    mkdir(TARGET_DIR)
    println("Created directory: $TARGET_DIR")
end

println("ᗜˬᗜ Fumofumo~ We have arrived at the target directory.")
cd(TARGET_DIR)

OUTPUT_DIR = "./data_Saving"

########################################
#=======NUMERICAL CONFIGURATION========#
########################################
nx, ny = 512, 512

stepper         = "ETDRK4"
"""
# step_per_second = 20000
const CODE2B = raw"""
# include a "step_per_sec(ond)" column in the CSV to include a non-default(5000) value for this
seconds_wished  = 60
snap_density    = 100

dt              = 1 / step_per_second
nsteps          = step_per_second * seconds_wished
nsubs           = step_per_second ÷ snap_density

"""
#=    is_continuation = false =#
const CODE3 = raw"""
data_origin = 1
"""
#=    continuation_point = 1044 =#
const CODE4 = raw"""

if is_continuation == false
    continuation_point = 0
end

########################################
#=======PHYSICAL PARAMETERS============#
########################################
Lx, Ly = 2π, 2π
"""
#=    mode = :fixedF1=#
const CODE5 = raw"""
K = 

"""

const CODE5U = raw"""
# if mode == :fixedU
U_fixed = 

"""#=3=#

const CODE5F = raw"""
# if mode == :fixedF1
F1_fixed = 

"""#=3=#

const CODE6 = raw"""

########################################
#=======N19N PROPERTIES================#
########################################
N19N_Length         = 1
N19N_Velocity       = 1
N19N_Time           = 1
N19N_StreamFunction = 1
N19N_Vorticity      = 1

########################################
#=======PACKAGES=======================#
########################################
include(joinpath(@__DIR__, "TwoLayerQGWithZonalPVFlux.jl"))
using .TwoLayerQGWithZonalPVFlux
using FourierFlows, Printf
using Random: seed!
using JLD2

########################################
#=======UTILS==========================#
########################################
function gaussf2dxz(Nx, Nz, ld)
    omegax = (0 : Nx ÷ 2) |> collect
    omegax = omegax'
    omegax = [omegax (-omegax[Nx ÷ 2:-1:2])']

    omegaz = (0 : Nz ÷ 2) |> collect
    omegaz = omegaz'
    omegaz = [omegaz (-omegaz[Nz ÷ 2:-1:2])']

    k = repeat(reshape(omegaz, 1, :), length(omegax), 1)
    l = repeat(reshape(omegax, :, 1), 1, length(omegaz))
    Karr = sqrt.(k .^ 2 .+ l .^ 2)

    sigma = sqrt(2 * pi) / ld
    Chat = exp.(-sigma^2 .* Karr.^2 ./ 2)

    f = fft(sqrt.(Chat) .* (randn(Nx, Nz) .+ 1im .* randn(Nx, Nz)))
    f = real(f)
    frms = sqrt(sum(f.^2) / Nx / Nz)
    f = f ./ frms
    return f
end

########################################
#=======PROBLEM DEFINITION=============#
########################################
dev = CPU()

if mode == :fixedU
    prob = TwoLayerQGWithZonalPVFlux.Problem(dev;
        nx, ny, Lx, Ly,
        mode=mode, U=U_fixed, K=K,
        dt=dt, stepper=stepper)
elseif mode == :fixedF1
    prob = TwoLayerQGWithZonalPVFlux.Problem(dev;
        nx, ny, Lx, Ly,
        mode=mode, F1=F1_fixed, K=K,
        dt=dt, stepper=stepper)
else
    error("Unknown mode = $mode")
end

println("ᗜˬᗜ Fumofumo~ The zonal-PV-flux QG problem has been set successfully.")

########################################
#=======SHORTCUTS======================#
########################################
sol, clock, params, vars, grid = prob.sol, prob.clock, prob.params, prob.vars, prob.grid

########################################
#=======INITIAL CONDITIONS=============#
########################################
if is_continuation == false
    seed!(1234)
    q10 = 0.05 * gaussf2dxz(nx, ny, floor(Int, 1.25*K))
    q10 .-= sum(q10) / nx / ny
    q20 = 0.05 * gaussf2dxz(nx, ny, floor(Int, 1.25*K))
    q20 .-= sum(q20) / nx / ny
else
    checkpoint_file = OUTPUT_DIR * "/raw/" * "t_" * string(continuation_point) * ".jld2"
    checkpoint_data = JLD2.load(checkpoint_file)
    q10 = checkpoint_data["VAR_q1"]
    q20 = checkpoint_data["VAR_q2"]
end

TwoLayerQGWithZonalPVFlux.set_q!(prob, q10, q20)
TwoLayerQGWithZonalPVFlux.updatevars!(prob; compute_streamfunctions=true)

println("ᗜˬᗜ Fumofumo~ Initial values have been set successfully.")

########################################
#=======DIAGNOSTICS====================#
########################################
E   = Diagnostic(TwoLayerQGWithZonalPVFlux.energies_tlqg, prob; nsteps)
FLX = Diagnostic(TwoLayerQGWithZonalPVFlux.pv_fluxes, prob; nsteps)
UBG = Diagnostic(TwoLayerQGWithZonalPVFlux.bg_velocity, prob; nsteps)

diags = [E, FLX, UBG]

########################################
#=======OUTPUT=========================#
########################################
if !isdir(OUTPUT_DIR)
    mkdir(OUTPUT_DIR)
end
if !isdir(OUTPUT_DIR * "/raw")
    mkdir(OUTPUT_DIR * "/raw")
end
if !isdir(OUTPUT_DIR * "/fig")
    mkdir(OUTPUT_DIR * "/fig")
end

JLD2.save(OUTPUT_DIR * "/t_parameters.jld2",
          "PARAMS",               prob.params,
          "MODE",                 String(mode),
          "U_FIXED",              U_fixed,
          "F1_FIXED",             F1_fixed,
          "K",                    K,
          "N19N_Length",          N19N_Length,
          "N19N_Velocity",        N19N_Velocity,
          "N19N_Time",            N19N_Time,
          "N19N_StreamFunction",  N19N_StreamFunction,
          "N19N_Vorticity",       N19N_Vorticity)

println("ᗜˬᗜ Fumofumo~ Timestepping...")
startwalltime = time()
frames = 0:round(Int, nsteps / nsubs)

########################################
#=======TIME-STEPPING==================#
########################################
for j = frames
    TwoLayerQGWithZonalPVFlux.updatevars!(prob; compute_streamfunctions=true)

    current_U = TwoLayerQGWithZonalPVFlux.bg_velocity(prob)
    current_F1, current_F2 = TwoLayerQGWithZonalPVFlux.pv_fluxes(prob)

    cfl = clock.dt * maximum([
        maximum(abs.(vars.u1)) / grid.dx,
        maximum(abs.(vars.v1)) / grid.dy,
        maximum(abs.(vars.u2)) / grid.dx,
        maximum(abs.(vars.v2)) / grid.dy
    ])

    log = @sprintf("step: %06d, t: %.6f, cfl: %.2f, U: %.6e, F1: %.6e, F2: %.6e, KE1: %.3e, KE2: %.3e, PE: %.3e, walltime: %.2f min",
                   clock.step, clock.t, cfl,
                   current_U, current_F1, current_F2,
                   E.data[E.i][1], E.data[E.i][2], E.data[E.i][3],
                   (time() - startwalltime)/60)
    println(log)

    if any(isnan.([current_U, current_F1, current_F2, E.data[E.i][1], E.data[E.i][2], E.data[E.i][3]]))
        println("A NaN error is figured out and the program hence aborts. Make your dt smaller!")
        break
    end

    JLD2.save(OUTPUT_DIR * "/raw/" * "t_" * string(j + continuation_point) * ".jld2",
              "VAR_q1",   prob.vars.q1,
              "VAR_q2",   prob.vars.q2,
              "VAR_psi1", prob.vars.ψ1,
              "VAR_psi2", prob.vars.ψ2,
              "E",        E.data[E.i],
              "F1",       current_F1,
              "F2",       current_F2,
              "U",        current_U,
              "CLOCK_STEP", clock.step,
              "CLOCK_T",    clock.t)

    println("$j: Saved data at $(E.i) at $(OUTPUT_DIR * "/raw")")

    FourierFlows.stepforward!(prob, diags, nsubs)
    TwoLayerQGWithZonalPVFlux.updatevars!(prob; compute_streamfunctions=true)
end


"""

# =========================
# 工具函数
# =========================

clean(x) = x === missing ? "" : strip(String(x))

function row_to_dict(row)
    d = Dict{String,String}()
    for name in propertynames(row)
        d[lowercase(String(name))] = clean(getproperty(row, name))
    end
    return d
end

function require_value(d::Dict{String,String}, name::String)
    key = lowercase(name)
    if !haskey(d, key) || isempty(d[key])
        error("CSV 缺少必填列或该列为空: $(name)")
    end
    return d[key]
end

function optional_value(d::Dict{String,String}, names::Vector{String}, default::String="")
    for name in names
        key = lowercase(name)
        if haskey(d, key) && !isempty(d[key])
            return d[key]
        end
    end
    return default
end

function normalize_bool_literal(s::String)
    t = lowercase(strip(s))
    if t in ["true", "t", "1", "yes", "y"]
        return "true"
    elseif t in ["false", "f", "0", "no", "n"]
        return "false"
    else
        error("is_continuation 只能写 true/false、1/0、yes/no；当前值: $(s)")
    end
end

function ensure_jl_filename(name::String)
    base = strip(name)
    base = replace(base, r"[\/\\:*?\"<>|]" => "_")
    return endswith(base, ".jl") ? base : base * ".jl"
end

file_token(x::String) = replace(x, "." => "p")

function generate_one(row, outdir::String)
    d = row_to_dict(row)

    mode = require_value(d, "mode")
    F = require_value(d, "F")
    U = require_value(d, "U")
    K = require_value(d, "K")
    step_per_sec = optional_value(d, ["step_per_sec","step_per_second"], "5000")

    is_continuation = normalize_bool_literal(
        require_value(d, "is_continuation")
    )

    continuation_point = optional_value(
        d,
        ["continuation_point", "continue_from", "cp"],
        "nothing",
    )

    if mode == "fixedF1"
        dir_name = "ZPVF_F$(file_token(F))_K$(file_token(K))"
    elseif mode == "fixedU"
        dir_name = "VS_U$(file_token(U))_K$(file_token(K))"
    end
    # 兼容 script_name / die_name / file_name / output_name
    script_name = optional_value(
        d,
        ["script_name", "die_name", "file_name", "output_name"],
        dir_name,
    )
    dir_name = script_name

    output_file = joinpath(outdir, ensure_jl_filename(script_name))

    content = string(
        CODE1, "\n",
        "TARGET_DIR = \"./$(dir_name)\"\n",
        CODE2, "\n",
        "step_per_second = $(step_per_sec)\n",
        CODE2B, "\n",
        "is_continuation = $(is_continuation)\n",
        CODE3, "\n",
        "continuation_point = $(continuation_point)\n",
        CODE4, "\n",
        "mode = :$(mode)\n",
        "K = $(K)\n",
        "U_fixed = $(U)\n",
        "F1_fixed =  -$(F)\n",
        CODE6, "\n",
    )

    mkpath(dirname(output_file))
    write(output_file, content)

    return output_file
end

const DEFAULT_CORES = 8

const SBATCH_HEADER = raw"""#!/bin/bash
#SBATCH -J julia_batch
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c {{CORES}}
#SBATCH --array={{ARRAY_SPEC}}
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err
"""

const SBATCH_PARAMETERS = raw"""
# ===== parameters =====
cores={{CORES}} # care the budget and the outcome
LOG_DIR="./logs"
mkdir -p "${LOG_DIR}"
"""

const SBATCH_CONSTRUCT_CASE = raw"""
# ===== constructing case =====
case1=${cases[$SLURM_ARRAY_TASK_ID]}

echo "Running case: ${case1}"
echo "Job ID: ${SLURM_JOB_ID}, Task ID: ${SLURM_ARRAY_TASK_ID}"
"""

const SBATCH_RUN = raw"""
# ===== run =====
srun julia --threads=${cores} "${case1}"
"""

function print_usage()
    println("""
    用法:
      julia generate_entries.jl params.csv [output_dir] [cores] [sbatch_path]

    参数:
      params.csv      参数 CSV 文件
      output_dir      入口脚本输出目录，默认 generated_scripts
      cores           每个 array task 使用的核心数，默认 $(DEFAULT_CORES)
      sbatch_path     sbatch 文件输出路径，默认 output_dir/run.sbatch

    示例:
      julia generate_entries.jl params.csv
      julia generate_entries.jl params.csv generated_scripts
      julia generate_entries.jl params.csv generated_scripts 16
      julia generate_entries.jl params.csv generated_scripts 16 submit.sbatch
    """)
end

function parse_main_args(args)
    if length(args) < 1 || length(args) > 4
        print_usage()
        exit(1)
    end

    csv_path = args[1]
    outdir = length(args) >= 2 ? args[2] : "generated_scripts"
    cores = length(args) >= 3 ? parse(Int, args[3]) : DEFAULT_CORES
    sbatch_path = length(args) >= 4 ? args[4] : joinpath(outdir, "run.sbatch")

    if cores <= 0
        error("cores 必须是正整数，目前得到: $(cores)")
    end

    return csv_path, outdir, cores, sbatch_path
end

function bash_single_quote(s::AbstractString)
    # Bash 单引号安全转义:
    # abc'def -> 'abc'\''def'
    return "'" * replace(String(s), "'" => "'\\''") * "'"
end

function render_array_spec(n::Integer)
    n <= 0 && error("没有生成任何任务，无法生成 sbatch array")
    return n == 1 ? "0" : "0-$(n - 1)"
end

function render_case_list(files::Vector{String})
    lines = String[]

    push!(lines, "# ===== case list =====")
    push!(lines, "cases=(")

    for f in files
        push!(lines, "    $(bash_single_quote(f))")
    end

    push!(lines, ")")

    return join(lines, "\n") * "\n"
end

function render_sbatch(files::Vector{String}, cores::Integer)
    array_spec = render_array_spec(length(files))

    header = replace(
        SBATCH_HEADER,
        "{{CORES}}" => string(cores),
        "{{ARRAY_SPEC}}" => array_spec,
    )

    parameters = replace(
        SBATCH_PARAMETERS,
        "{{CORES}}" => string(cores),
    )

    case_list = render_case_list(files)

    return join(
        [
            header,
            parameters,
            case_list,
            SBATCH_CONSTRUCT_CASE,
            SBATCH_RUN,
        ],
        "\n",
    )
end

function write_sbatch(sbatch_path::AbstractString, files::Vector{String}, cores::Integer)
    sbatch_text = render_sbatch(files, cores)

    parent = dirname(sbatch_path)
    if !isempty(parent)
        mkpath(parent)
    end

    write(sbatch_path, sbatch_text)

    return sbatch_path
end


function main()
    csv_path, outdir, cores, sbatch_path = parse_main_args(ARGS)

    mkpath(outdir)

    files = String[]

    for row in CSV.File(csv_path; types=String, normalizenames=true)
        push!(files, generate_one(row, outdir))
    end

    write_sbatch(sbatch_path, files, cores)

    println("生成完成，共 $(length(files)) 个入口脚本：")
    for f in files
        println("  ", f)
    end

    println()
    println("已生成 sbatch 脚本：")
    println("  ", sbatch_path)
    println("核心数 cores：", cores)
    println("array 范围：", render_array_spec(length(files)))
end


main()
