module TwoLayerQGWithZonalPVFlux

using FourierFlows
using LinearAlgebra: mul!, ldiv!
using Statistics: mean

export Problem, set_q!, updatevars!,
       energies_tlqg, pv_fluxes, bg_velocity,
       Params, Vars

# ============================================================
# Model:
#
# (q1)_t + U ( (q1)_x + kD^2 (ψ1)_x ) + [ψ1,q1] - overline([ψ1,q1]) = Δ q1
# (q2)_t - U ( (q2)_x + kD^2 (ψ2)_x ) + [ψ2,q2] - overline([ψ2,q2]) = Δ q2
#
# q_i = Δ ψ_i + (1/2) kD^2 ( ψ_{3-i} - ψ_i )
#
# PV fluxes:
# F1 = overline( (ψ1)_x q1 ) - kD^2 U
# F2 = overline( (ψ2)_x q2 ) + kD^2 U
#
# Two modes:
# 1) :fixedU   -> U fixed, F1/F2 diagnostic outputs
# 2) :fixedF1  -> F1 fixed, U = ( <ψ1_x q1> - F1 ) / kD^2, F2 diagnostic output
# ============================================================

# ----------------
# Parameters
# ----------------
struct Params{T, A<:AbstractArray} <: AbstractParams
    mode::Symbol          # :fixedU or :fixedF1
    U::T                  # used if mode = :fixedU
    F1::T                 # used if mode = :fixedF1
    K2::T
    invKrsq::A
end

function Params_(mode::Symbol, U, F1, K, grid::TwoDGrid)
    Dev = typeof(grid.device)
    T   = eltype(grid)

    K2 = T(K)^2
    @devzeros Dev T (grid.nkr, grid.nl) invKrsq
    @. invKrsq = ifelse(grid.Krsq == 0, zero(T), 1 / grid.Krsq)

    return Params(mode, T(mode == :fixedU ? U : 0),
                  T(mode == :fixedF1 ? F1 : 0),
                  K2,
                  invKrsq)
end

# ----------------
# Variables
# ----------------
mutable struct Vars{Aphys, Ah, T} <: AbstractVars
    q1::Aphys
    q2::Aphys
    ψ1::Aphys
    ψ2::Aphys
    u1::Aphys
    v1::Aphys
    u2::Aphys
    v2::Aphys
    fluxx::Aphys
    fluxy::Aphys

    q1h::Ah
    q2h::Ah
    ψ1h::Ah
    ψ2h::Ah
    u1h::Ah
    v1h::Ah
    u2h::Ah
    v2h::Ah
    fluxxh::Ah
    fluxyh::Ah

    Uinst::T
    F1inst::T
    F2inst::T
end

function Vars(grid::TwoDGrid)
    Dev = typeof(grid.device)
    T   = eltype(grid)

    @devzeros Dev T (grid.nx, grid.ny) q1 q2 ψ1 ψ2 u1 v1 u2 v2 fluxx fluxy
    @devzeros Dev Complex{T} (grid.nkr, grid.nl) q1h q2h ψ1h ψ2h u1h v1h u2h v2h fluxxh fluxyh

    return Vars(q1, q2, ψ1, ψ2, u1, v1, u2, v2, fluxx, fluxy,
                q1h, q2h, ψ1h, ψ2h, u1h, v1h, u2h, v2h, fluxxh, fluxyh,
                zero(T), zero(T), zero(T))
end

# ----------------
# FFT wrappers
# ----------------
fwdtransform!(varh, var, grid::TwoDGrid) = mul!(varh, grid.rfftplan, var)
invtransform!(var, varh, grid::TwoDGrid) = ldiv!(var, grid.rfftplan, varh)

# ----------------
# Streamfunction inversion
# q1 = Δψ1 + 1/2 K2 (ψ2 - ψ1)
# q2 = Δψ2 + 1/2 K2 (ψ1 - ψ2)
#
# In Fourier:
# q1h = -(Krsq + K2/2) ψ1h + (K2/2) ψ2h
# q2h = (K2/2) ψ1h - (Krsq + K2/2) ψ2h
# ----------------
function invert_streamfunctions!(vars::Vars, params::Params, grid::TwoDGrid)
    K2 = params.K2
    @. vars.q1h = vars.q1h
    @. vars.q2h = vars.q2h

    A = grid.Krsq .+ K2/2
    B = K2/2
    det = A .* A .- B .* B

    @. vars.ψ1h = (-A * vars.q1h - B * vars.q2h) / det
    @. vars.ψ2h = (-B * vars.q1h - A * vars.q2h) / det

    # zero mode gauge
    vars.ψ1h[1,1] = 0
    vars.ψ2h[1,1] = 0
    return nothing
end

# ----------------
# Update helper
# ----------------
function updatevars!(vars::Vars, params::Params, grid::TwoDGrid, sol;
                     compute_streamfunctions::Bool=false)

    @. vars.q1h = sol[:, :, 1]
    @. vars.q2h = sol[:, :, 2]

    invert_streamfunctions!(vars, params, grid)

    invtransform!(vars.q1, vars.q1h, grid)
    invtransform!(vars.q2, vars.q2h, grid)

    kr = reshape(grid.kr, grid.nkr, 1)
    l  = reshape(grid.l,  1, grid.nl)

    @. vars.u1h = -im * l  * vars.ψ1h
    @. vars.v1h =  im * kr * vars.ψ1h
    @. vars.u2h = -im * l  * vars.ψ2h
    @. vars.v2h =  im * kr * vars.ψ2h

    invtransform!(vars.u1, vars.u1h, grid)
    invtransform!(vars.v1, vars.v1h, grid)
    invtransform!(vars.u2, vars.u2h, grid)
    invtransform!(vars.v2, vars.v2h, grid)

    if compute_streamfunctions
        invtransform!(vars.ψ1, vars.ψ1h, grid)
        invtransform!(vars.ψ2, vars.ψ2h, grid)
    end

    # instantaneous U and PV fluxes
    if params.mode == :fixedU
        U = params.U
    elseif params.mode == :fixedF1
        # F1 = < ψ1_x q1 > - K2 U
        # U = ( <ψ1_x q1> - F1 ) / K2
        # ψ1_x = v1
        U = (mean(vars.v1 .* vars.q1) - params.F1) / params.K2
    else
        error("Unknown mode $(params.mode)")
    end

    vars.Uinst = U
    vars.F1inst = mean(vars.v1 .* vars.q1) - params.K2 * U
    vars.F2inst = mean(vars.v2 .* vars.q2) + params.K2 * U

    return nothing
end

function updatevars!(prob::FourierFlows.Problem; compute_streamfunctions::Bool=false)
    return updatevars!(prob.vars, prob.params, prob.grid, prob.sol;
                       compute_streamfunctions=compute_streamfunctions)
end

# ----------------
# Nonlinear RHS
# ----------------
function calcN!(N, sol, t, clock, vars::Vars, params::Params, grid::TwoDGrid)

    dealias!(sol, grid)
    updatevars!(vars, params, grid, sol; compute_streamfunctions=false)

    kr = reshape(grid.kr, grid.nkr, 1)
    l  = reshape(grid.l,  1, grid.nl)

    @views q1h  = sol[:, :, 1]
    @views q2h  = sol[:, :, 2]
    @views N1h  = N[:, :, 1]
    @views N2h  = N[:, :, 2]

    U = vars.Uinst
    K2 = params.K2

    # ---------- layer 1 ----------
    # bracket [ψ1,q1] = ∂x(u1 q1) + ∂y(v1 q1)
    @. vars.fluxx = vars.u1 * vars.q1
    @. vars.fluxy = vars.v1 * vars.q1
    fwdtransform!(vars.fluxxh, vars.fluxx, grid)
    fwdtransform!(vars.fluxyh, vars.fluxy, grid)

    @. N1h = -(im * kr * vars.fluxxh + im * l * vars.fluxyh)

    # subtract zonal mean of bracket: kx = 0 modes removed
    @views N1h[1, :] .= 0

    # add linear background-advection term: +U(q1_x + K2 ψ1_x)
    @. N1h += -U * (im * kr * q1h + K2 * im * kr * vars.ψ1h)

    # ---------- layer 2 ----------
    @. vars.fluxx = vars.u2 * vars.q2
    @. vars.fluxy = vars.v2 * vars.q2
    fwdtransform!(vars.fluxxh, vars.fluxx, grid)
    fwdtransform!(vars.fluxyh, vars.fluxy, grid)

    @. N2h = -(im * kr * vars.fluxxh + im * l * vars.fluxyh)
    @views N2h[1, :] .= 0

    # add linear term: -U(q2_x + K2 ψ2_x)
    @. N2h += +U * (im * kr * q2h + K2 * im * kr * vars.ψ2h)

    dealias!(N, grid)
    return nothing
end

# ----------------
# Equation & Problem
# ----------------
function Equation(params::Params, grid::TwoDGrid)
    dev = grid.device
    T   = eltype(grid)

    L = zeros(dev, T, (grid.nkr, grid.nl, 2))
    @views L[:, :, 1] .= -grid.Krsq
    @views L[:, :, 2] .= -grid.Krsq

    return FourierFlows.Equation(L, calcN!, grid)
end

function Problem(dev::Device=CPU(); nx=128, ny=nx, Lx=2π, Ly=Lx,
                 mode::Symbol=:fixedU, U=0.0, F1=0.0, K=1.0,
                 dt=1e-3, stepper="ETDRK4",
                 aliased_fraction=1/3, T=Float64)

    grid     = TwoDGrid(dev; nx, Lx, ny, Ly, aliased_fraction, T)
    params   = Params_(mode, U, F1, K, grid)
    vars     = Vars(grid)
    equation = Equation(params, grid)

    return FourierFlows.Problem(equation, stepper, dt, grid, vars, params)
end

# ----------------
# Initial conditions
# ----------------
function set_q!(prob::FourierFlows.Problem, q10, q20)
    vars, grid, sol = prob.vars, prob.grid, prob.sol
    A = typeof(vars.q1)

    fwdtransform!(vars.q1h, A(q10), grid)
    fwdtransform!(vars.q2h, A(q20), grid)

    @views sol[:, :, 1] .= vars.q1h
    @views sol[:, :, 2] .= vars.q2h

    updatevars!(prob; compute_streamfunctions=false)
    return nothing
end

# ----------------
# Diagnostics
# ----------------
function energies_tlqg(prob::FourierFlows.Problem)
    updatevars!(prob; compute_streamfunctions=true)
    v = prob.vars
    grid = prob.grid
    K2 = prob.params.K2

    KE1 = 0.5 * sum(@. v.u1^2 + v.v1^2) * grid.dx * grid.dy
    KE2 = 0.5 * sum(@. v.u2^2 + v.v2^2) * grid.dx * grid.dy
    PE  = 0.25 * K2 * sum(@. (v.ψ1 - v.ψ2)^2) * grid.dx * grid.dy
    return (KE1, KE2, PE)
end

function pv_fluxes(prob::FourierFlows.Problem)
    updatevars!(prob; compute_streamfunctions=false)
    return (prob.vars.F1inst, prob.vars.F2inst)
end

function bg_velocity(prob::FourierFlows.Problem)
    updatevars!(prob; compute_streamfunctions=false)
    return prob.vars.Uinst
end

end # module
