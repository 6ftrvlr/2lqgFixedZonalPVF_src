########################################
#====== SpectralEnergySeries.jl =======#
########################################
using CairoMakie
using JLD2
using Printf
using FFTW
using Statistics
using LinearAlgebra
using Colors

########################################
#============= COLORMAPS ==============#
########################################
function build_yuyuko(nColor::Int=256)
    nColor <= 0 && error("nColor 必须是正整数")

    blue    = [0.40 0.50 0.95]
    neutral = [0.95 0.95 0.95]
    red     = [0.95 0.50 0.80]
    color_nodes = [blue; neutral; red]

    nNode = size(color_nodes, 1)

    if nColor == 1
        c = color_nodes[1, :]
        return [RGB(c[1], c[2], c[3])]
    end

    yuyuko = Vector{RGB{Float64}}(undef, nColor)

    for k in 1:nColor
        # x 从 0 到 nNode - 1
        x = (k - 1) / (nColor - 1) * (nNode - 1)

        # 当前所在色段
        i = min(floor(Int, x) + 1, nNode - 1)

        # 当前色段内部插值参数，范围 0 到 1
        t = x - (i - 1)

        c1 = color_nodes[i, :]
        c2 = color_nodes[i + 1, :]

        r = (1 - t) * c1[1] + t * c2[1]
        g = (1 - t) * c1[2] + t * c2[2]
        b = (1 - t) * c1[3] + t * c2[3]

        yuyuko[k] = RGB(r, g, b)
    end

    return yuyuko
end


function positive_yuyuko(nColor::Int=256)
    nColor <= 0 && error("nColor 必须是正整数")

    blue    = [0.40 0.50 0.95]
    neutral = [0.80 0.60 0.85]
    red     = [0.95 0.50 0.80]
    color_nodes = [blue; neutral; red]

    nNode = size(color_nodes, 1)

    if nColor == 1
        c = color_nodes[1, :]
        return [RGB(c[1], c[2], c[3])]
    end

    yuyuko = Vector{RGB{Float64}}(undef, nColor)

    for k in 1:nColor
        # x 从 0 到 nNode - 1
        x = (k - 1) / (nColor - 1) * (nNode - 1)

        # 当前所在色段
        i = min(floor(Int, x) + 1, nNode - 1)

        # 当前色段内部插值参数，范围 0 到 1
        t = x - (i - 1)

        c1 = color_nodes[i, :]
        c2 = color_nodes[i + 1, :]

        r = (1 - t) * c1[1] + t * c2[1]
        g = (1 - t) * c1[2] + t * c2[2]
        b = (1 - t) * c1[3] + t * c2[3]

        yuyuko[k] = RGB(r, g, b)
    end

    return yuyuko
end

function sample_curve_colors(n::Int)
    cmap = positive_yuyuko(512)
    idx = round.(Int, range(1, length(cmap), length=n))
    return cmap[idx]
end

########################################
#============= UTILITIES ==============#
########################################
function parse_mode1_args()
    if length(ARGS) < 10
        error("""
        模式 1 用法:
            julia SpectralEnergySeries.jl 1 file_dir start_num end_num num_stride dt N N_1 t_0 k_D

        例如:
            julia SpectralEnergySeries.jl 1 data 0 3255 1 0.005 20 8 0.0 10.0
        """)
    end
    file_dir   = ARGS[2]
    start_num  = parse(Int, ARGS[3])
    end_num    = parse(Int, ARGS[4])
    num_stride = parse(Int, ARGS[5])
    dt         = parse(Float64, ARGS[6])
    N          = parse(Int, ARGS[7])
    N_1        = parse(Int, ARGS[8])
    t_0        = parse(Float64, ARGS[9])
    kD         = parse(Float64, ARGS[10])

    isdir(file_dir)      || error("找不到文件夹: $file_dir")
    start_num <= end_num || error("start_num 必须 <= end_num")
    num_stride > 0       || error("num_stride 必须 > 0")
    dt > 0               || error("dt 必须 > 0")
    N > 0                || error("N 必须 > 0")
    N_1 > 0              || error("N_1 必须 > 0")
    kD >= 0              || error("k_D 必须 >= 0")

    return file_dir, start_num, end_num, num_stride, dt, N, N_1, t_0, kD
end

function parse_mode2_args()
    if length(ARGS) < 3
        error("""
        模式 2 用法:
            julia SpectralEnergySeries.jl 2 spectral_energy_jld2_path N_1
        """)
    end
    infile = ARGS[2]
    N_1    = parse(Int, ARGS[3])
    isfile(infile) || error("找不到文件: $infile")
    N_1 > 0        || error("N_1 必须 > 0")
    return infile, N_1
end

function ensure_real_matrix(A, key::String, num::Int)
    M = Array(A)
    ndims(M) == 2 || error("文件 t_$num.jld2 中的 $key 不是二维矩阵")
    all(isreal, M) || error("文件 t_$num.jld2 中的 $key 含有复数")
    return Float64.(real.(M))
end

function load_snapshot(infile::String, num::Int)
    isfile(infile) || error("找不到文件: $infile")
    data = JLD2.load(infile)

    q1   = ensure_real_matrix(data["VAR_q1"],   "VAR_q1",   num)
    q2   = ensure_real_matrix(data["VAR_q2"],   "VAR_q2",   num)
    psi1 = ensure_real_matrix(data["VAR_psi1"], "VAR_psi1", num)
    psi2 = ensure_real_matrix(data["VAR_psi2"], "VAR_psi2", num)

    U = haskey(data, "U") ? Float64(data["U"]) : NaN

    return (q1=q1, q2=q2, psi1=psi1, psi2=psi2, U=U)
end

# 全复谱 fft 对应波数
function fftfreq_like_matched(n::Int, L::Float64)
    n > 0 || error("n 必须 > 0")
    L > 0 || error("L 必须 > 0")

    kk = Vector{Float64}(undef, n)

    if iseven(n)
        nh = n ÷ 2
        for i in 1:(nh + 1)
            kk[i] = i - 1
        end
        for i in (nh + 2):n
            kk[i] = i - 1 - n
        end
    else
        nh = (n - 1) ÷ 2
        for i in 1:(nh + 1)
            kk[i] = i - 1
        end
        for i in (nh + 2):n
            kk[i] = i - 1 - n
        end
    end

    return (2π / L) .* kk
end

# 按你的数据约定修正：
# 半谱实际尺寸是 (Nx÷2+1, Ny)，也就是第一维压缩，第二维完整
# 因而这里返回第一维对应的 rfft 波数
function rfftfreq_like_matched(n::Int, L::Float64)
    n > 0 || error("n 必须 > 0")
    L > 0 || error("L 必须 > 0")

    nr = fld(n, 2) + 1
    kk = Vector{Float64}(undef, nr)
    for i in 1:nr
        kk[i] = i - 1
    end
    return (2π / L) .* kk
end

function spectral_grid_rfft(Nx::Int, Ny::Int; Lx=2π, Ly=2π)
    # 对应半谱尺寸：(Nx÷2+1, Ny)
    kx = rfftfreq_like_matched(Nx, Lx)   # 长度 Nx÷2+1，对应第1维
    ky = fftfreq_like_matched(Ny, Ly)    # 长度 Ny，对应第2维

    Nxh = fld(Nx, 2) + 1

    length(kx) == Nxh || error("length(kx) = $(length(kx)) != Nxh = $Nxh")
    length(ky) == Ny  || error("length(ky) = $(length(ky)) != Ny = $Ny")

    KX = repeat(reshape(kx, Nxh, 1), 1, Ny)
    KY = repeat(reshape(ky, 1, Ny), Nxh, 1)
    K2 = KX.^2 .+ KY.^2
    K  = sqrt.(K2)

    return kx, ky, KX, KY, K2, K
end

function modal_energy_from_streamfunctions(psi1::Matrix{Float64}, psi2::Matrix{Float64}, kD::Float64; Lx=2π, Ly=2π)
    size(psi1) == size(psi2) || error("psi1 与 psi2 尺寸不一致")

    Nx, Ny = size(psi1)

    Nx > 0 || error("Nx 必须 > 0")
    Ny > 0 || error("Ny 必须 > 0")
    iseven(Nx) || error("为保证 2^n 网格与 rfft 半谱安全，Nx 必须为偶数；当前 Nx = $Nx")
    iseven(Ny) || error("为保证 2^n 网格与逻辑一致，Ny 必须为偶数；当前 Ny = $Ny")

    psibt = 0.5 .* (psi1 .+ psi2)
    psibc = 0.5 .* (psi1 .- psi2)

    # 按你的数据约定，半谱需要是 (Nx÷2+1, Ny) = (257, 512)
    # 因此对第1维做 rfft
    ψbth = rfft(psibt)
    ψbch = rfft(psibc)

    expected_size = (fld(Nx, 2) + 1, Ny)
    size(ψbth) == expected_size || error("rfft(psibt, 1) 尺寸异常: $(size(ψbth))，预期 = $expected_size")
    size(ψbch) == expected_size || error("rfft(psibc, 1) 尺寸异常: $(size(ψbch))，预期 = $expected_size")

    kx, ky, KX, KY, K2, K = spectral_grid_rfft(Nx, Ny; Lx=Lx, Ly=Ly)

    size(K2) == size(ψbth) || error("K2 尺寸 $(size(K2)) 与 ψbth 尺寸 $(size(ψbth)) 不一致")
    size(K2) == size(ψbch) || error("K2 尺寸 $(size(K2)) 与 ψbch 尺寸 $(size(ψbch)) 不一致")

    normfac = (Nx * Ny)^2

    Ebt   = 0.5 .* K2 .* abs2.(ψbth) ./ normfac
    Ebc   = 0.5 .* K2 .* abs2.(ψbch) ./ normfac
    Eape  = 0.5 .* (kD^2) .* abs2.(ψbch) ./ normfac

    Ebt[1,1]  = 0.0
    Ebc[1,1]  = 0.0
    Eape[1,1] = 0.0

    return Ebt, Ebc, Eape, kx, ky, K
end

function mode_label(ix::Int, iy::Int, kx, ky)
    return @sprintf("(k=%.3f, l=%.3f)", kx[ix], ky[iy])
end

function top_modes_by_time_mean(E3d::Array{Float64,3}, N::Int)
    @show size(E3d)
    Emean = dropdims(mean(E3d, dims=1), dims=1)
    @show size(Emean)
    inds = vec(CartesianIndices(Emean))
    vals = [Emean[I] for I in inds]
    p = sortperm(vals, rev=true)
    topI = inds[p[1:min(N, length(p))]]
    return topI, Emean
end

function union_mode_sets(mode_sets...)
    d = Dict{Tuple{Int,Int},Nothing}()
    for ms in mode_sets
        for I in ms
            d[(I[1], I[2])] = nothing
        end
    end
    return collect(keys(d))
end

function build_series_from_modes(E3d::Array{Float64,3}, mode_set::Vector{Tuple{Int,Int}})
    nt = size(E3d, 1)
    Nxh = size(E3d, 2)
    Ny  = size(E3d, 3)

    out = Dict{Tuple{Int,Int}, Vector{Float64}}()
    for mode in mode_set
        ix, iy = mode
        1 <= ix <= Nxh || error("mode ix = $ix 越界，允许范围 1:$Nxh")
        1 <= iy <= Ny  || error("mode iy = $iy 越界，允许范围 1:$Ny")
        out[mode] = [E3d[it, ix, iy] for it in 1:nt]
    end
    return out
end

function compute_ylim_from_topN(series_dict::Dict{Tuple{Int,Int},Vector{Float64}},
                                mode_strength_order::Vector{Tuple{Int,Int}},
                                N_1::Int)
    sel = mode_strength_order[1:min(N_1, length(mode_strength_order))]
    vals = Float64[]
    for m in sel
        append!(vals, series_dict[m])
    end
    vals = filter(x -> isfinite(x) && x > 0, vals)

    if isempty(vals)
        return (1e-16, 1.0)
    end

    lv = log10.(vals)
    m = minimum(lv)
    M = maximum(lv)
    δ = M - m
    if δ == 0
        return (10.0^(m - 0.5), 10.0^(M + 0.5))
    else
        return (10.0^(m - 0.1*δ), 10.0^(M + 0.1*δ))
    end
end

function sort_modes_by_mean_strength(series_dict::Dict{Tuple{Int,Int},Vector{Float64}})
    modes = collect(keys(series_dict))
    vals = [mean(series_dict[m]) for m in modes]
    p = sortperm(vals, rev=true)
    return modes[p]
end

function isotropic_shell_spectrum(E2d::Matrix{Float64}, K::Matrix{Float64})
    size(E2d) == size(K) || error("E2d 尺寸 $(size(E2d)) 与 K 尺寸 $(size(K)) 不一致")

    Kmax = floor(Int, maximum(K))
    spec = zeros(Float64, Kmax + 1)
    counts = zeros(Int, Kmax + 1)

    Nxh, Ny = size(E2d)
    for iy in 1:Ny, ix in 1:Nxh
        n = round(Int, K[ix, iy])
        if 0 <= n <= Kmax
            spec[n+1] += E2d[ix, iy]
            counts[n+1] += 1
        end
    end
    return collect(0:Kmax), spec, counts
end

function isotropic_shell_spectrum(E2d::Matrix{Float64}, K::Matrix{Float64})
    @assert size(E2d) == size(K) "E2d and K must have the same size"
    dk = 0.1   # 可在这里改成你想要的固定分桶间隔
    Kmax = maximum(K)
    nbins = floor(Int, Kmax / dk) + 1
    k_shell = collect(0:nbins-1) .* dk
    E_shell = zeros(Float64, nbins)
    for I in CartesianIndices(E2d)
        kval = K[I]
        e = E2d[I]
        ibin = floor(Int, kval / dk) + 1
        if 1 <= ibin <= nbins
            E_shell[ibin] += e
        end
    end
    return k_shell, E_shell
end


function plot_series(time, series_dict, mode_order, kx, ky, title_str, outfile, N_1)
    fig = Figure(size=(1100, 700))
    ax = Axis(fig[1,1],
        title=title_str,
        xlabel="t",
        ylabel="modal energy",
        yscale=log10
    )

    colors = sample_curve_colors(length(mode_order))

    for (i, mode) in enumerate(mode_order)
        ix, iy = mode
        lines!(ax, time, series_dict[mode],
            color=colors[i],
            linewidth=2,
            label=mode_label(ix, iy, kx, ky))
    end

    ylims!(ax, compute_ylim_from_topN(series_dict, mode_order, N_1))
    xlims!(ax, minimum(time), maximum(time))
    axislegend(ax, position=:rb, labelsize=12)
    save(outfile, fig)
end

function plot_mean_spectra_2d(Ebt_mean, Ebc_mean, Eape_mean, outfile; Lk=2π, Ll=2π, cmp=positive_yuyuko(512))
    nk, nl = size(Ebt_mean)
    uniform_maxwn = 30;

    @assert size(Ebc_mean) == (nk, nl) "Ebc_mean size mismatch"
    @assert size(Eape_mean) == (nk, nl) "Eape_mean size mismatch"

    # 第1维：半谱方向，只保留非负波数
    l = collect(0:nk-1) .* (2π / Lk)

    # 第2维：完整谱方向，先构造 FFTW 顺序，再做 shift
    l_raw = vcat(0:(nl ÷ 2 - 1), -nl ÷ 2:-1)
    p = fftshift(1:nl)
    k = l_raw[p] .* (2π / Ll)

    # 只对完整谱方向做 shift
    Ebt_plot  = Ebt_mean[:, p]
    Ebc_plot  = Ebc_mean[:, p]
    Eape_plot = Eape_mean[:, p]


    fig = Figure(size=(2400, 1200))

    ticks = -12:2:2
    ticklabels = [L"10^{%$t}" for t in ticks]

    ax1 = Axis(
        fig[1, 1],
        title = "mean BT EKE 2D spectrum",
        xlabel = "k",
        ylabel = "l",
        aspect = nk / nl,
    )
    xlims!(ax1, 0, uniform_maxwn)
    ylims!(ax1, -uniform_maxwn, uniform_maxwn)

    hm1 = heatmap!(ax1, l, k, log10.(Ebt_plot .+ eps()), colorrange=(minimum(ticks), maximum(ticks)), colormap = cmp)
    Colorbar(fig[1, 2], hm1, ticks=(ticks, ticklabels))

    ax2 = Axis(
        fig[1, 3],
        title = "mean BC EKE 2D spectrum",
        xlabel = "k",
        ylabel = "l",
        aspect = nk / nl,
    )
    xlims!(ax2, 0, uniform_maxwn)
    ylims!(ax2, -uniform_maxwn, uniform_maxwn)

    hm2 = heatmap!(ax2, l, k, log10.(Ebc_plot .+ eps()), colorrange =(minimum(ticks), maximum(ticks)), colormap = cmp)
    Colorbar(fig[1, 4], hm2, ticks=(ticks, ticklabels))

    ax3 = Axis(
        fig[1, 5],
        title = "mean EAPE 2D spectrum",
        xlabel = "k",
        ylabel = "l",
        aspect = nk / nl,
    )
    xlims!(ax3, 0, uniform_maxwn)
    ylims!(ax3, -uniform_maxwn, uniform_maxwn) 

    hm3 = heatmap!(ax3, l, k, log10.(Eape_plot .+ eps()), colorrange =(minimum(ticks), maximum(ticks)), colormap = cmp)
    Colorbar(fig[1, 6], hm3, ticks=(ticks, ticklabels))

    save(outfile, fig)
end

function plot_isotropic_spectra(nshell, bt_spec, bc_spec, ape_spec, outfile)
    fig = Figure(size=(1000, 700))
    ax = Axis(fig[1,1],
        title="isotropic shell-integrated mean spectra",
        xlabel="shell index n",
        ylabel="shell-integrated energy",
        yscale=log10
    )

    lines!(ax, nshell, bt_spec, linewidth=2, label="BT EKE")
    lines!(ax, nshell, bc_spec, linewidth=2, label="BC EKE")
    lines!(ax, nshell, ape_spec, linewidth=2, label="EAPE")
    
    ylims!(ax, (10.0^-5,10.0^5))

    axislegend(ax, position=:rb)
    save(outfile, fig)
end

########################################
#========== MODE-1 ANALYSIS ===========#
########################################
function run_mode1()
    file_dir, start_num, end_num, num_stride, dt, N, N_1, t_0, kD = parse_mode1_args()

    raw_dir = joinpath(file_dir, "raw")
    fig_dir = joinpath(file_dir, "fig")
    isdir(raw_dir) || error("当前目录下找不到 $raw_dir")
    mkpath(fig_dir)

    nums = Int[]
    Ebt_list  = Array{Float64,2}[]
    Ebc_list  = Array{Float64,2}[]
    Eape_list = Array{Float64,2}[]
    kx = nothing
    ky = nothing
    K  = nothing

    println("Reading snapshots from: $(abspath(raw_dir))")

    expected_phys_size = nothing

    for num in start_num:num_stride:end_num
        infile = joinpath(raw_dir, @sprintf("t_%d.jld2", num))
        if !isfile(infile)
            @warn "跳过缺失文件" infile
            continue
        end

        snap = load_snapshot(infile, num)

        size(snap.psi1) == size(snap.psi2) || error("t_$num: psi1 与 psi2 尺寸不一致")
        Nx, Ny = size(snap.psi1)

        if expected_phys_size === nothing
            expected_phys_size = (Nx, Ny)
        else
            size(snap.psi1) == expected_phys_size || error("t_$num: 物理空间尺寸 $(size(snap.psi1)) 与前面快照 $expected_phys_size 不一致")
        end

        iseven(Nx) || error("t_$num: Nx = $Nx 不是偶数；为保证 rfft 半谱安全，输入必须是偶数长度")
        iseven(Ny) || error("t_$num: Ny = $Ny 不是偶数；为保证 2^n × 2^n 逻辑一致，输入必须是偶数长度")

        #@show num size(snap.psi1) size(snap.psi2)

        Ebt, Ebc, Eape, kx_, ky_, K_ = modal_energy_from_streamfunctions(snap.psi1, snap.psi2, kD)

        expected_half_size = (fld(Nx,2)+1, Ny)
        size(Ebt) == size(Ebc) == size(Eape) || error("t_$num: 三个谱能量矩阵尺寸不一致")
        size(Ebt) == expected_half_size || error("t_$num: 半谱尺寸异常，得到 $(size(Ebt))，预期 = $expected_half_size")

        push!(nums, num)
        push!(Ebt_list, Ebt)
        push!(Ebc_list, Ebc)
        push!(Eape_list, Eape)

        if kx === nothing
            kx = kx_
            ky = ky_
            K  = K_
        else
            length(kx_) == length(kx) || error("t_$num: kx 长度变化")
            length(ky_) == length(ky) || error("t_$num: ky 长度变化")
            size(K_) == size(K) || error("t_$num: K 尺寸变化")
        end

        println(@sprintf("Loaded spectral energies from t_%d", num))
    end

    isempty(nums) && error("指定范围内没有成功读取到任何 jld2 文件")

    nt = length(nums)
    Nxh, Ny = size(Ebt_list[1])

    Ebt_3d  = Array{Float64,3}(undef, nt, Nxh, Ny)
    Ebc_3d  = Array{Float64,3}(undef, nt, Nxh, Ny)
    Eape_3d = Array{Float64,3}(undef, nt, Nxh, Ny)

    for it in 1:nt
        size(Ebt_list[it])  == (Nxh, Ny) || error("第 $it 个 Ebt 尺寸不一致")
        size(Ebc_list[it])  == (Nxh, Ny) || error("第 $it 个 Ebc 尺寸不一致")
        size(Eape_list[it]) == (Nxh, Ny) || error("第 $it 个 Eape 尺寸不一致")

        Ebt_3d[it, :, :]  = Ebt_list[it]
        Ebc_3d[it, :, :]  = Ebc_list[it]
        Eape_3d[it, :, :] = Eape_list[it]
    end

    bt_top,   Ebt_mean  = top_modes_by_time_mean(Ebt_3d, N)
    bc_top,   Ebc_mean  = top_modes_by_time_mean(Ebc_3d, N)
    ape_top,  Eape_mean = top_modes_by_time_mean(Eape_3d, N)

    union_modes = union_mode_sets(bt_top, bc_top, ape_top)

    bt_series  = build_series_from_modes(Ebt_3d, union_modes)
    bc_series  = build_series_from_modes(Ebc_3d, union_modes)
    ape_series = build_series_from_modes(Eape_3d, union_modes)

    bt_order  = sort_modes_by_mean_strength(bt_series)
    bc_order  = sort_modes_by_mean_strength(bc_series)
    ape_order = sort_modes_by_mean_strength(ape_series)

    time = t_0 .+ ((nums .- first(nums)) .* dt)

    nshell, bt_iso  = isotropic_shell_spectrum(Ebt_mean, K)
    _,      bc_iso  = isotropic_shell_spectrum(Ebc_mean, K)
    _,      ape_iso = isotropic_shell_spectrum(Eape_mean, K)

    outjld = joinpath(fig_dir, @sprintf("spectral_energy_series_%d_%d.jld2", first(nums), last(nums)))
    JLD2.save(outjld,
        "mode", "snapshot_to_spectral_series",
        "nums", nums,
        "time", time,
        "kD", kD,
        "kx", kx,
        "ky", ky,
        "N", N,
        "N_1", N_1,
        "union_modes", union_modes,
        "bt_order", bt_order,
        "bc_order", bc_order,
        "ape_order", ape_order,
        "bt_series", bt_series,
        "bc_series", bc_series,
        "ape_series", ape_series,
        "Ebt_mean", Ebt_mean,
        "Ebc_mean", Ebc_mean,
        "Eape_mean", Eape_mean,
        "nshell", nshell,
        "bt_iso", bt_iso,
        "bc_iso", bc_iso,
        "ape_iso", ape_iso
    )
    println("Saved spectral energy data to: $outjld")

    plot_series(time, bt_series,  bt_order,  kx, ky, "Barotropic modal kinetic energy", joinpath(fig_dir, "BT_spectral_energy_series.png"), N_1)
    plot_series(time, bc_series,  bc_order,  kx, ky, "Baroclinic modal kinetic energy", joinpath(fig_dir, "BC_spectral_energy_series.png"), N_1)
    plot_series(time, ape_series, ape_order, kx, ky, "Modal EAPE",                        joinpath(fig_dir, "EAPE_spectral_energy_series.png"), N_1)

    plot_mean_spectra_2d(Ebt_mean, Ebc_mean, Eape_mean, joinpath(fig_dir, "mean_2D_spectra.png"))
    plot_isotropic_spectra(nshell, bt_iso, bc_iso, ape_iso, joinpath(fig_dir, "mean_isotropic_spectra.png"))

    println("Saved figures to: $(abspath(fig_dir))")
    println("Done.")
end

########################################
#========== MODE-2 REDRAW =============#
########################################
function run_mode2()
    infile, N_1 = parse_mode2_args()
    data = JLD2.load(infile)

    time       = Vector{Float64}(data["time"])
    kx         = Vector{Float64}(data["kx"])
    ky         = Vector{Float64}(data["ky"])
    bt_series  = data["bt_series"]
    bc_series  = data["bc_series"]
    ape_series = data["ape_series"]
    bt_order   = data["bt_order"]
    bc_order   = data["bc_order"]
    ape_order  = data["ape_order"]

    Ebt_mean   = Matrix{Float64}(data["Ebt_mean"])
    Ebc_mean   = Matrix{Float64}(data["Ebc_mean"])
    Eape_mean  = Matrix{Float64}(data["Eape_mean"])
    nshell     = Vector{Float64}(data["nshell"])
    bt_iso     = Vector{Float64}(data["bt_iso"])
    bc_iso     = Vector{Float64}(data["bc_iso"])
    ape_iso    = Vector{Float64}(data["ape_iso"])

    size(Ebt_mean) == size(Ebc_mean) == size(Eape_mean) || error("重绘数据中的二维平均谱尺寸不一致")

    fig_dir = dirname(infile)

    plot_series(time, bt_series,  bt_order,  kx, ky, "Barotropic modal kinetic energy", joinpath(fig_dir, "BT_spectral_energy_series_redraw.png"), N_1)
    plot_series(time, bc_series,  bc_order,  kx, ky, "Baroclinic modal kinetic energy", joinpath(fig_dir, "BC_spectral_energy_series_redraw.png"), N_1)
    plot_series(time, ape_series, ape_order, kx, ky, "Modal EAPE",                        joinpath(fig_dir, "EAPE_spectral_energy_series_redraw.png"), N_1)

    plot_mean_spectra_2d(Ebt_mean, Ebc_mean, Eape_mean, joinpath(fig_dir, "mean_2D_spectra_redraw.png"))
    plot_isotropic_spectra(nshell, bt_iso, bc_iso, ape_iso, joinpath(fig_dir, "mean_isotropic_spectra_redraw.png"))

    println("Redrawn figures in: $(abspath(fig_dir))")
end

########################################
#=============== MAIN =================#
########################################
function main()
    length(ARGS) >= 1 || error("请先给出模式编号 1 或 2")
    mode = parse(Int, ARGS[1])
    if mode == 1
        run_mode1()
    elseif mode == 2
        run_mode2()
    else
        error("模式编号只能是 1 或 2")
    end
end

main()