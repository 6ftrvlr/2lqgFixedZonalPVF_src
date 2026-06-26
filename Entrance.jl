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

TARGET_DIR = "./data/F72000_SampTest"

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

stepper            = "ETDRK4"

step_per_second    = 200000
# include a "step_per_sec(ond)" column in the CSV to include a non-default(5000) value for this
seconds_wished     = 60
snap_density       = 1              # snapshots per second
diag_print_density = 200000               # diagnostic printouts per second (NEW)

dt                 = 1 / step_per_second
nsteps             = step_per_second * seconds_wished

# snapshot interval in steps
nsubs_snap         = step_per_second ÷ snap_density

# diagnostic-print interval in steps
nsubs_diag         = step_per_second ÷ diag_print_density

# unified outer-loop interval in steps: minimal change while decoupling print/snapshot cadence
outer_nsubs        = gcd(nsubs_snap, nsubs_diag)
nouter             = nsteps ÷ outer_nsubs

is_continuation = true
data_origin = 1

continuation_point = 13000

if is_continuation == false
    continuation_point = 0
end

########################################
#=======PHYSICAL PARAMETERS============#
########################################
Lx, Ly = 2π, 2π

mode = :fixedF1
K = 30
U_fixed = 0
F1_fixed = -67000

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
          "STEP_PER_SECOND",      step_per_second,
          "SNAP_DENSITY",         snap_density,
          "DIAG_PRINT_DENSITY",   diag_print_density,
          "N19N_Length",          N19N_Length,
          "N19N_Velocity",        N19N_Velocity,
          "N19N_Time",            N19N_Time,
          "N19N_StreamFunction",  N19N_StreamFunction,
          "N19N_Vorticity",       N19N_Vorticity)

println("ᗜˬᗜ Fumofumo~ Timestepping...")
startwalltime = time()

########################################
#=======HELPERS========================#
########################################
function print_diagnostics!(prob, E, startwalltime)
    vars  = prob.vars
    grid  = prob.grid
    clock = prob.clock

    TwoLayerQGWithZonalPVFlux.updatevars!(prob; compute_streamfunctions=true)

    current_U = TwoLayerQGWithZonalPVFlux.bg_velocity(prob)
    current_F1, current_F2 = TwoLayerQGWithZonalPVFlux.pv_fluxes(prob)

    cfl = clock.dt * maximum([
        maximum(abs.(vars.u1)) / grid.dx,
        maximum(abs.(vars.v1)) / grid.dy,
        maximum(abs.(vars.u2)) / grid.dx,
        maximum(abs.(vars.v2)) / grid.dy
    ])

    log = @sprintf(
        "step: %06d, t: %.6f, cfl: %.2f, U: %.6e, F1: %.6e, F2: %.6e, KE1: %.3e, KE2: %.3e, PE: %.3e, walltime: %.2f min",
        clock.step, clock.t, cfl,
        current_U, current_F1, current_F2,
        E.data[E.i][1], E.data[E.i][2], E.data[E.i][3],
        (time() - startwalltime) / 60
    )
    println(log)

    if any(isnan.([current_U, current_F1, current_F2, E.data[E.i][1], E.data[E.i][2], E.data[E.i][3]]))
        println("A NaN error is figured out and the program hence aborts. Make your dt smaller!")
        return false
    end
    return true
end

function save_snapshot!(prob, E, output_dir, snapshot_index)
    clock = prob.clock

    TwoLayerQGWithZonalPVFlux.updatevars!(prob; compute_streamfunctions=true)
    current_U = TwoLayerQGWithZonalPVFlux.bg_velocity(prob)
    current_F1, current_F2 = TwoLayerQGWithZonalPVFlux.pv_fluxes(prob)

    JLD2.save(output_dir * "/raw/" * "t_" * string(snapshot_index) * ".jld2",
              "VAR_q1",      prob.vars.q1,
              "VAR_q2",      prob.vars.q2,
              "VAR_psi1",    prob.vars.ψ1,
              "VAR_psi2",    prob.vars.ψ2,
              "E",           E.data[E.i],
              "F1",          current_F1,
              "F2",          current_F2,
              "U",           current_U,
              "CLOCK_STEP",  clock.step,
              "CLOCK_T",     clock.t)

    println("Saved snapshot $(snapshot_index) at step $(clock.step) to $(output_dir * "/raw")")
    return nothing
end

########################################
#=======TIME-STEPPING==================#
########################################
# Keep snapshot cadence unchanged, but decouple diagnostic print cadence.
# Advance with the smallest common interval and trigger print/save only when needed.

# optional: preserve the old behavior of output at initial time
if !print_diagnostics!(prob, E, startwalltime)
    error("NaN detected before timestepping starts.")
end
save_snapshot!(prob, E, OUTPUT_DIR, continuation_point)

snap_count = 0

for outer = 1:nouter
    FourierFlows.stepforward!(prob, diags, outer_nsubs)

    # print diagnostics at user-specified cadence
    if mod(clock.step, nsubs_diag) == 0
        ok = print_diagnostics!(prob, E, startwalltime)
        if !ok
            break
        end
    end

    # save snapshots at original cadence
    if mod(clock.step, nsubs_snap) == 0
        snap_count += 1
        save_snapshot!(prob, E, OUTPUT_DIR, snap_count + continuation_point)
    end
end
