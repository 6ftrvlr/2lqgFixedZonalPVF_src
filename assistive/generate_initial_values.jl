#!/usr/bin/env julia

using FFTW
using JLD2
using Random
using Printf

############################################################
# Help / usage
############################################################

function print_usage_and_exit(code::Int=0)
    println("""
Usage:
  julia Generate_initial_values.jl F kD k_sat l_sat pert_type [pert_args...] [options]

Required positional arguments:
  F           : Float64   -- prescribed F1 / background flux parameter
  kD          : Int       -- deformation wavenumber
  k_sat       : Int       -- background saturated-state x-wavenumber
  l_sat       : Int       -- background saturated-state y-wavenumber
  pert_type   : Int       -- 1 (Fourier-mode perturbations) or 2 (Gaussian perturbation)

If pert_type == 1:
  Perturbation arguments must be given in blocks of 6 numbers:
    k_mode l_mode amp1 phase1 amp2 phase2
  meaning:
    k_mode : Int
    l_mode : Int
    amp1   : Float64   amplitude added to q1 on this mode
    phase1 : Float64   phase added to q1 on this mode
    amp2   : Float64   amplitude added to q2 on this mode
    phase2 : Float64   phase added to q2 on this mode

  Example:
    julia Generate_initial_values.jl -67000 48 30 0 1 \\
      1 0 1e-6 0.0 2e-6 1.57079632679 \\
      2 1 5e-7 0.3 5e-7 -0.2

If pert_type == 2:
  Perturbation arguments:
    amp_gauss relParam_gauss

  meaning:
    amp_gauss      : Float64   RMS prefactor of Gaussian random perturbation
    relParam_gauss : Float64   correlation-scale parameter, ld = floor(Int, relParam_gauss * kD)

  Example:
    julia Generate_initial_values.jl -67000 48 30 0 2 1e-4 1.25

Optional options:
  --nx N            grid size in x (default: 512)
  --ny N            grid size in y (default: 512)
  --seed S          RNG seed for reproducibility (default: 1234)
  --outdir DIR      output directory (default: .)
  --outfile FILE    output file name (default: t_0.jld2)
  --force           overwrite existing output file
  --help, -h        print this help

Output file content is written in a format readable by Entrance.jl continuation:
  VAR_q1, VAR_q2, VAR_psi1, VAR_psi2
plus metadata fields.
""")
    exit(code)
end

############################################################
# Utilities
############################################################

function gaussf2dxz(Nx::Int, Nz::Int, ld::Int)
    ld >= 1 || error("gaussf2dxz: ld must be >= 1, got $ld")

    omegax = (0:Nx ÷ 2) |> collect
    omegax = omegax'
    omegax = [omegax (-omegax[Nx ÷ 2:-1:2])']

    omegaz = (0:Nz ÷ 2) |> collect
    omegaz = omegaz'
    omegaz = [omegaz (-omegaz[Nz ÷ 2:-1:2])']

    k = repeat(reshape(omegaz, 1, :), length(omegax), 1)
    l = repeat(reshape(omegax, :, 1), 1, length(omegaz))
    Karr = sqrt.(k.^2 .+ l.^2)

    sigma = sqrt(2π) / ld
    Chat = exp.(-sigma^2 .* Karr.^2 ./ 2)

    f = fft(sqrt.(Chat) .* (randn(Nx, Nz) .+ 1im .* randn(Nx, Nz)))
    f = real(f)
    frms = sqrt(sum(f.^2) / Nx / Nz)
    frms > 0 || error("gaussf2dxz: zero RMS encountered unexpectedly")
    f ./= frms
    return f
end

function l_to_col(l::Int, Ny::Int)
    col = l >= 0 ? l + 1 : l + Ny + 1
    (1 <= col <= Ny) || error("Invalid y-wavenumber l = $l for Ny = $Ny")
    return col
end

# add Re{A exp(i(kx+ly))} into rfft storage
function add_real_mode!(uh::AbstractMatrix{ComplexF64}, k::Int, l::Int, A::ComplexF64, Nx::Int, Ny::Int)
    0 <= k <= Nx ÷ 2 || error("For rfft storage, x-wavenumber k must satisfy 0 <= k <= Nx/2, got k = $k with Nx = $Nx")
    col = l_to_col(l, Ny)
    uh[k + 1, col] += A * Nx * Ny / 2
    return nothing
end

function ensure_parent_dir(path::String)
    dir = dirname(path)
    if !isdir(dir)
        mkpath(dir)
    end
    return nothing
end

function parse_float(s::String, name::String)
    try
        return parse(Float64, s)
    catch
        error("Failed to parse $name as Float64: \"$s\"")
    end
end

function parse_int(s::String, name::String)
    try
        return parse(Int, s)
    catch
        error("Failed to parse $name as Int: \"$s\"")
    end
end

############################################################
# Parse command line
############################################################

if any(x -> x == "--help" || x == "-h", ARGS)
    print_usage_and_exit(0)
end

length(ARGS) >= 5 || print_usage_and_exit(1)

# separate positional and optional arguments
first_opt = findfirst(x -> startswith(x, "--"), ARGS)
positional = first_opt === nothing ? copy(ARGS) : ARGS[1:first_opt-1]
options    = first_opt === nothing ? String[]    : ARGS[first_opt:end]

length(positional) >= 5 || error("Need at least 5 positional arguments. Use --help for usage.")

F         = parse_float(positional[1], "F")
kD        = parse_int(positional[2], "kD")
k_sat     = parse_int(positional[3], "k_sat")
l_sat     = parse_int(positional[4], "l_sat")
pert_type = parse_int(positional[5], "pert_type")

# defaults
Nx      = 512
Ny      = 512
seed    = 1234
outdir  = "."
outfile = "t_0.jld2"
force   = false

# parse options
i = 1
while i <= length(options)
    opt = options[i]
    if opt == "--nx"
        i + 1 <= length(options) || error("--nx requires an integer value")
        Nx = parse_int(options[i+1], "--nx")
        i += 2
    elseif opt == "--ny"
        i + 1 <= length(options) || error("--ny requires an integer value")
        Ny = parse_int(options[i+1], "--ny")
        i += 2
    elseif opt == "--seed"
        i + 1 <= length(options) || error("--seed requires an integer value")
        seed = parse_int(options[i+1], "--seed")
        i += 2
    elseif opt == "--outdir"
        i + 1 <= length(options) || error("--outdir requires a directory path")
        outdir = options[i+1]
        i += 2
    elseif opt == "--outfile"
        i + 1 <= length(options) || error("--outfile requires a file name")
        outfile = options[i+1]
        i += 2
    elseif opt == "--force"
        force = true
        i += 1
    else
        error("Unknown option: $opt. Use --help for usage.")
    end
end

Nx > 0 || error("Nx must be positive, got $Nx")
Ny > 0 || error("Ny must be positive, got $Ny")
iseven(Nx) || error("Nx must be even for this rfft-based script, got Nx = $Nx")
iseven(Ny) || error("Ny must be even for this script's spectral indexing, got Ny = $Ny")
kD > 0 || error("kD must be positive, got $kD")
k_sat != 0 || error("k_sat must be nonzero because U = sqrt(S) * (k^2+l^2)/k")
abs(l_sat) <= Ny ÷ 2 || error("Background l_sat = $l_sat exceeds representable range for Ny = $Ny")
0 <= k_sat <= Nx ÷ 2 || error("Background k_sat = $k_sat must satisfy 0 <= k_sat <= Nx/2 = $(Nx ÷ 2)")
(pert_type == 1 || pert_type == 2) || error("pert_type must be 1 or 2, got $pert_type")

pert_args = positional[6:end]

Random.seed!(seed)

############################################################
# Build saturated background state
############################################################

κ2 = k_sat^2 + l_sat^2
kD^2 > κ2 || error("Need kD^2 > k_sat^2 + l_sat^2. Got kD^2 = $(kD^2), κ2 = $κ2")

S = (kD^2 + κ2) / (kD^2 - κ2)
S > 0 || error("Computed S must be positive, got S = $S")

U = sqrt(S) * κ2 / k_sat

A1_sq = -(2 * (1 + S^2) / (k_sat * kD^2 * S)) * (F + kD^2 * U)
isfinite(A1_sq) || error("Computed A1^2 is not finite: $A1_sq")
A1_sq > 0 || error("Computed A1^2 <= 0. Got A1^2 = $A1_sq. Check parameter regime.")

A1 = sqrt(A1_sq) + 0im
A2 = -A1 * (1 + 1im*S) / (1 - 1im*S)

Q1p = -(κ2 + 0.5*kD^2) * A1 + 0.5*kD^2 * A2
Q2p = -(κ2 + 0.5*kD^2) * A2 + 0.5*kD^2 * A1

all(isfinite.([real(A1), imag(A1), real(A2), imag(A2), real(Q1p), imag(Q1p), real(Q2p), imag(Q2p), U])) ||
    error("Non-finite background quantities encountered")

# spectral arrays
q1h = zeros(ComplexF64, Nx ÷ 2 + 1, Ny)
q2h = zeros(ComplexF64, Nx ÷ 2 + 1, Ny)
ψ1h = zeros(ComplexF64, Nx ÷ 2 + 1, Ny)
ψ2h = zeros(ComplexF64, Nx ÷ 2 + 1, Ny)

add_real_mode!(ψ1h, k_sat, l_sat, A1,  Nx, Ny)
add_real_mode!(ψ2h, k_sat, l_sat, A2,  Nx, Ny)
add_real_mode!(q1h, k_sat, l_sat, Q1p, Nx, Ny)
add_real_mode!(q2h, k_sat, l_sat, Q2p, Nx, Ny)

############################################################
# Add perturbation
############################################################

if pert_type == 1
    length(pert_args) > 0 || error("For pert_type = 1, you must provide at least one Fourier mode block")
    length(pert_args) % 6 == 0 || error("For pert_type = 1, perturbation args must come in groups of 6: k l amp1 phase1 amp2 phase2")

    nmodes = length(pert_args) ÷ 6

    for j in 0:nmodes-1
        km     = parse_int(pert_args[1 + 6j], "k_mode[$(j+1)]")
        lm     = parse_int(pert_args[2 + 6j], "l_mode[$(j+1)]")
        amp1   = parse_float(pert_args[3 + 6j], "amp1[$(j+1)]")
        phase1 = parse_float(pert_args[4 + 6j], "phase1[$(j+1)]")
        amp2   = parse_float(pert_args[5 + 6j], "amp2[$(j+1)]")
        phase2 = parse_float(pert_args[6 + 6j], "phase2[$(j+1)]")

        0 <= km <= Nx ÷ 2 || error("Mode $(j+1): k_mode = $km must satisfy 0 <= k <= Nx/2 = $(Nx ÷ 2)")
        abs(lm) <= Ny ÷ 2 || error("Mode $(j+1): l_mode = $lm exceeds representable range for Ny = $Ny")
        amp1 >= 0 || error("Mode $(j+1): amp1 must be nonnegative, got $amp1")
        amp2 >= 0 || error("Mode $(j+1): amp2 must be nonnegative, got $amp2")
        isfinite(phase1) || error("Mode $(j+1): phase1 is not finite")
        isfinite(phase2) || error("Mode $(j+1): phase2 is not finite")

        A_q1 = ComplexF64(amp1 * exp(1im * phase1))
        A_q2 = ComplexF64(amp2 * exp(1im * phase2))

        add_real_mode!(q1h, km, lm, A_q1, Nx, Ny)
        add_real_mode!(q2h, km, lm, A_q2, Nx, Ny)
    end

elseif pert_type == 2
    length(pert_args) == 2 || error("For pert_type = 2, need exactly 2 perturbation arguments: amp_gauss relParam_gauss")

    amp_gauss      = parse_float(pert_args[1], "amp_gauss")
    relParam_gauss = parse_float(pert_args[2], "relParam_gauss")

    amp_gauss >= 0 || error("amp_gauss must be nonnegative, got $amp_gauss")
    relParam_gauss > 0 || error("relParam_gauss must be positive, got $relParam_gauss")

    ld = floor(Int, relParam_gauss * kD)
    ld >= 1 || error("Computed ld = floor(relParam_gauss * kD) must be >= 1, got ld = $ld")

    q1p = amp_gauss * gaussf2dxz(Nx, Ny, ld)
    q2p = amp_gauss * gaussf2dxz(Nx, Ny, ld)

    q1p .-= sum(q1p) / (Nx * Ny)
    q2p .-= sum(q2p) / (Nx * Ny)

    q1h .+= rfft(q1p)
    q2h .+= rfft(q2p)
end

############################################################
# Transform to physical space
############################################################

q1 = irfft(q1h, Nx, (1,2))
q2 = irfft(q2h, Nx, (1,2))
ψ1 = irfft(ψ1h, Nx, (1,2))
ψ2 = irfft(ψ2h, Nx, (1,2))

size(q1) == (Nx, Ny) || error("Internal error: q1 has wrong size")
size(q2) == (Nx, Ny) || error("Internal error: q2 has wrong size")
size(ψ1) == (Nx, Ny) || error("Internal error: ψ1 has wrong size")
size(ψ2) == (Nx, Ny) || error("Internal error: ψ2 has wrong size")

all(isfinite, q1) || error("Non-finite values found in q1")
all(isfinite, q2) || error("Non-finite values found in q2")
all(isfinite, ψ1) || error("Non-finite values found in ψ1")
all(isfinite, ψ2) || error("Non-finite values found in ψ2")

############################################################
# Output
############################################################

isdir(outdir) || mkpath(outdir)
outpath = joinpath(outdir, outfile)

if isfile(outpath) && !force
    error("Output file already exists: $outpath. Use --force to overwrite.")
end

ensure_parent_dir(outpath)

E_placeholder = ([0.0, 0.0, 0.0], [0.0])

JLD2.save(outpath,
    "VAR_q1",     q1,
    "VAR_q2",     q2,
    "VAR_psi1",   ψ1,
    "VAR_psi2",   ψ2,
    "F1",         F,
    "U",          U,
    "A1",         A1,
    "A2",         A2,
    "Q1p",        Q1p,
    "Q2p",        Q2p,
    "kD",         kD,
    "k_sat",      k_sat,
    "l_sat",      l_sat,
    "PERT_TYPE",  pert_type,
    "NX",         Nx,
    "NY",         Ny,
    "SEED",       seed,
    "E",          E_placeholder
)

@printf("Saved initial condition to: %s\n", outpath)
@printf("Grid: Nx = %d, Ny = %d\n", Nx, Ny)
@printf("Background: F1 = %.12e, kD = %d, k_sat = %d, l_sat = %d, U = %.12e\n", F, kD, k_sat, l_sat, U)