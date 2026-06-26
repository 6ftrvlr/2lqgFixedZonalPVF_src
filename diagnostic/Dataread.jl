########################################
#=========== PACKAGES =================#
########################################
using CairoMakie
using JLD2
using CSV
using DataFrames
using Printf
using Colors   # 若 reisana 来自 ColorSchemes，需要这一行
include("Aesthetics.jl")
########################################
#=========== USER CONFIG ==============#
########################################
# 运行方式（在 data_0 文件夹下）:
#   julia Dataread.jl data 0 3255 1 0.005
#
# 表示读取:
#   data/raw/t_0.jld2, data/raw/t_1.jld2, ..., data/raw/t_3255.jld2
#
# 输出图片到:
#   data/fig/
#
# 输出诊断量 CSV 到:
#   data/fig/

const CMAP_Q   = reisana
const CMAP_P   = reisana
const CMAP_PSI = reisana
const CMAP_PHI = reisana

########################################
#=========== BASIC UTILS ==============#
########################################
function parse_args()
    if length(ARGS) < 5
        error("""
        用法:
            julia Dataread.jl file_dir start_num end_num num_stride dt

        例如:
            julia Dataread.jl data 0 3255 1 0.005
        """)
    end

    file_dir   = ARGS[1]
    start_num  = parse(Int, ARGS[2])
    end_num    = parse(Int, ARGS[3])
    num_stride = parse(Int, ARGS[4])
    dt         = parse(Float64, ARGS[5])

    isdir(file_dir)      || error("找不到文件夹: $file_dir")
    start_num <= end_num || error("start_num 必须 <= end_num")
    num_stride > 0       || error("num_stride 必须 > 0")
    dt > 0               || error("dt 必须 > 0")

    return file_dir, start_num, end_num, num_stride, dt
end

function ensure_real_matrix(A, key::String, num::Int)
    M = Array(A)
    ndims(M) == 2 || error("文件 t_$num.jld2 中的 $key 不是二维矩阵")

    all(x -> x isa Number, M) || error("文件 t_$num.jld2 中的 $key 含有非数值元素")

    if !all(isreal, M)
        error("文件 t_$num.jld2 中的 $key 含有复数，当前脚本只直接画实数矩阵")
    end

    return Float64.(real.(M))
end

function extract_energies(Eobj, num::Int)
    # 预期:
    # E = ([KE1, KE2], [PE])
    length(Eobj) >= 3 || error("文件 t_$num.jld2 中的 E 格式不符合预期")

    KE1 = Eobj[1]
    KE2 = Eobj[2]
    PE  = Eobj[3]

    # length(KE_block) >= 2 || error("文件 t_$num.jld2 中的 E[1] 长度不足 2")

    # KE1 = Float64(KE_block[1])
    # KE2 = Float64(KE_block[2])

    # PE = if PE_block isa Number
    #     Float64(PE_block)
    # else
    #     Float64(first(PE_block))
    # end

    return KE1, KE2, PE
end

function symmetric_range(A)
    vmax = maximum(abs, A)
    vmax = vmax == 0 ? 1.0 : Float64(vmax)
    return (-vmax, vmax)
end

########################################
#=========== IO =======================#
########################################
function load_snapshot(infile::String, num::Int)
    isfile(infile) || error("找不到文件: $infile")

    data = JLD2.load(infile)

    Q   = ensure_real_matrix(data["VAR_q1"],   "VAR_q1",   num)
    P   = ensure_real_matrix(data["VAR_q2"],   "VAR_q2",   num)
    PSI = ensure_real_matrix(data["VAR_psi1"], "VAR_psi1", num)
    PHI = ensure_real_matrix(data["VAR_psi2"], "VAR_psi2", num)

    EHF = Float64(data["F1"])
    U = Float64(data["U"])
    KE1, KE2, PE = extract_energies(data["E"], num)

    return (
        Q   = Q,
        P   = P,
        PSI = PSI,
        PHI = PHI,
        U   = U,
        EHF = EHF,
        KE1 = KE1,
        KE2 = KE2,
        PE  = PE,
    )
end

########################################
#=========== PLOTTING =================#
########################################
function build_figure(snap, num::Int)
    
    fig = Figure(size = (1024, 768))

    Label(fig[0, :], @sprintf("t_%d", num), fontsize = 26)

    Lx, Ly = 2*pi, 2*pi
    Nx, Ny = size(snap.Q)[2], size(snap.Q)[1]
    dx, dy = Lx/Nx, Ly/Ny
    x,y    = (0:dx:dx*(Nx-1)).-Lx/2, (0:dy:dy*(Ny-1)).-Ly/2

    function symmetric_range(O)
        max_val = maximum(abs, O)
        limit = max_val == 0 ? 1.0 : max_val
        return (-limit, limit)
    end

    Q_range   = symmetric_range(snap.Q)
    PSI_range = symmetric_range(snap.PSI)
    P_range   = symmetric_range(snap.P)
    PHI_range = symmetric_range(snap.PHI)


    axis_kwargs = (xlabel = L"x/L_x",
                ylabel = L"y/L_y",
                aspect = 1,
                limits = ((-Lx/2, Lx/2), (-Ly/2, Ly/2)),
                alignmode = Mixed(left = 0, bottom = 0)
                )

    axQ = Axis(fig[1, 1][1, 1]; title = L"q_1^*", axis_kwargs...) # 
    axPSI = Axis(fig[2, 1][1, 1]; title = L"\psi_1^*", axis_kwargs...)
    axP = Axis(fig[1, 2][1, 1]; title = L"q_2^*", axis_kwargs...)
    axPHI = Axis(fig[2, 2][1, 1]; title = L"\phi_2^*", axis_kwargs...)

    cb_kwargs = (
        width = 15,           
        ticklabelspace = 35.0, 
        tellwidth = true    
    )

    hmQ = heatmap!(axQ, x, y, snap.Q; colormap = CMAP_Q, colorrange = Q_range)
    Colorbar(fig[1, 1][1, 2], hmQ; cb_kwargs...)
    hmP = heatmap!(axP, x, y, snap.P; colormap = CMAP_P, colorrange = P_range)
    Colorbar(fig[1, 2][1, 2], hmP; cb_kwargs...)
    hmPSI = heatmap!(axPSI, x, y, snap.PSI; colormap = CMAP_PSI, colorrange = PSI_range)
    Colorbar(fig[2, 1][1, 2], hmPSI; cb_kwargs...)
    hmPHI = heatmap!(axPHI, x, y, snap.PHI; colormap = CMAP_PHI, colorrange = PHI_range)
    Colorbar(fig[2, 2][1, 2], hmPHI; cb_kwargs...)

    colsize!(fig.layout, 1, Fixed(320))
    colsize!(fig.layout, 2, Fixed(320))

    resize_to_layout!(fig)

    return fig
end

########################################
#=========== MAIN =====================#
########################################
function main()
    file_dir, start_num, end_num, num_stride, dt = parse_args()

    raw_dir = joinpath(file_dir, "raw")
    fig_dir = joinpath(file_dir, "fig")

    isdir(raw_dir) || error("当前目录下找不到 $raw_dir")
    mkpath(fig_dir)

    nums_loaded = Int[]
    EHF_arr = Float64[]
    U_arr = Float64[]
    KE1_arr = Float64[]
    KE2_arr = Float64[]
    PE_arr  = Float64[]

    println("Start reading snapshots from: $(abspath(raw_dir))")
    println("Figures will be saved to:     $(abspath(fig_dir))")

    for num in start_num:num_stride:end_num
        infile = joinpath(raw_dir, @sprintf("t_%d.jld2", num))

        if !isfile(infile)
            @warn "跳过缺失文件" infile
            continue
        end

        snap = load_snapshot(infile, num)

        fig = build_figure(snap, num)
        outfile = joinpath(fig_dir, @sprintf("t_%d.png", num))
        save(outfile, fig)

        push!(nums_loaded, num)
        push!(EHF_arr, snap.EHF)
        push!(U_arr, snap.U)
        push!(KE1_arr, snap.KE1)
        push!(KE2_arr, snap.KE2)
        push!(PE_arr,  snap.PE)

        println(@sprintf("Saved figure: %s", outfile))
    end

    isempty(nums_loaded) && error("指定范围内没有成功读取到任何 jld2 文件")

    # 时间数组从 0 开始，步长为 dt
    t_arr = collect(0:length(nums_loaded)-1) .* dt

    # 如果你想让时间直接和文件编号对应，可改成:
    # t_arr = Float64.(nums_loaded) .* dt

    df = DataFrame(
        number = nums_loaded,
        t      = t_arr,
        EHF    = EHF_arr,
        U      = U_arr,
        KE1    = KE1_arr,
        KE2    = KE2_arr,
        PE     = PE_arr,
    )

    outcsv = joinpath(
        fig_dir,
        @sprintf("diagnostics_%d_%d.csv", first(nums_loaded), last(nums_loaded))
    )
    CSV.write(outcsv, df)

    println(@sprintf("Saved diagnostics CSV: %s", outcsv))
    println("Done.")
end

main()
