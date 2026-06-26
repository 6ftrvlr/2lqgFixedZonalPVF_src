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
TARGET_DIR = "./ZPVF_F3000_K30"

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
step_per_second = 20000
seconds_wished  = 60
snap_density    = 100

dt              = 1 / step_per_second
nsteps          = step_per_second * seconds_wished
nsubs           = step_per_second ÷ snap_density

is_continuation = true
data_origin     = 1
continuation_point = 1044

if is_continuation == false
    continuation_point = 0
end

########################################
#=======PHYSICAL PARAMETERS============#
########################################
Lx, Ly = 2π, 2π

# choose one:
mode = :fixedF1
# mode = :fixedU

K = 30.0

# if mode == :fixedU
U_fixed = 3

# if mode == :fixedF1
F1_fixed = 3000.0

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
