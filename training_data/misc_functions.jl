using Plots, CSV, DataFrames, ScatteredInterpolation, JLD2

include("dat_to_sdf.jl")

function plot_random_space(;n=2^8, save=false)
    AoA = rand(n)
    AoA = (AoA .- 0.5) * 40

    logRe = rand(n)
    logRe = logRe*4 .+ 2 
    Re = exp10.(logRe)

    n1 = rand(0:6, n)
    n2 = rand(0:6, n)
    n34 = rand(1:24, n)
    naca = 1000n1 .+ 100n2 .+ n34

    R = [AoA Re naca]

    f = scatter(AoA[1:Int(7/8*n)], logRe[1:Int(7/8*n)], naca[1:Int(7/8*n)]; 
        legend=false, size=(600, 600), markerstrokewidth=0.5)

    scatter!(AoA[Int(7/8*n):end], logRe[Int(7/8*n):end], naca[Int(7/8*n):end]; 
        markershape=:utriangle, markerstrokewidth=0.5)

    yticks!(2:6, ["10²", "10³", "10⁴", "10⁵", "10⁶"])
    zticks!([0006, 2412, 4424, 6409], ["0006", "2412", "4424", "6409"])
    xlabel!("α")
    ylabel!("Re")
    zlabel!("naca")

    if save
        savefig("Paper/figures/sample_space.png")
    end
    return f
end

function get_itp(df::DataFrame, column)
    max = maximum(df[!, column])
    min = minimum(df[!, column])

    itp = ScatteredInterpolation.interpolate(ThinPlate(), [df.x df.y]', df[!, column])

    return itp, min, max
end

function get_grid(df::DataFrame, column, n::Int)
    max = maximum(df[!, column])
    min = minimum(df[!, column])

    itp = ScatteredInterpolation.interpolate(ThinPlate(), [df.x df.y]', df[!, column])

    A = [evaluate(itp, [i,j])[1] for i in LinRange(-0.5, 1.5, n), j in LinRange(-0.5, 0.5, n)]

    return itp, A, min, max
end

function plotUVP(path::String; n=128)
    df = CSV.File(path) |> DataFrame
    rename!(df, ["Velocity: Magnitude (m/s)", "Velocity[i] (m/s)", "Velocity[j] (m/s)", 
        "Pressure (Pa)", "Vorticity: Magnitude (/s)", "X (m)", "Y (m)", "Z (m)"] .=> 
        ["Vm", "Vx", "Vy", "P", "Vorticity", "x", "y", "z"])

    subset!(df, :x => x -> -0.5 .<= x .<= 1.5, :y => y -> -0.5 .<= y .<= 0.5)

    xs = LinRange(-0.5, 1.5, n)
    ys = LinRange(-0.5, 0.5, n)

    Vx_i, Vx, Vx_min, Vx_max = get_grid(df, :Vx, n)
    Vy_i, Vy, Vy_min, Vy_max = get_grid(df, :Vy, n)
    P_i, P, P_min, P_max = get_grid(df, :P, n)

    return Vx_i,Vx, Vx_min, Vx_max, Vy_i, Vy, Vy_min, Vy_max, P_i, P, P_min, P_max
end

function interpolate_star(path::String)
    df = CSV.File(path) |> DataFrame
    rename!(df, ["Velocity: Magnitude (m/s)", "Velocity[i] (m/s)", "Velocity[j] (m/s)", 
        "Pressure (Pa)", "Vorticity: Magnitude (/s)", "X (m)", "Y (m)", "Z (m)"] .=> 
        ["Vm", "Vx", "Vy", "P", "Vorticity", "x", "y", "z"])

    subset!(df, :x => x -> -0.5 .<= x .<= 1.5, :y => y -> -0.5 .<= y .<= 0.5)

    Vx_i, Vx_min, Vx_max = get_itp(df, :Vx)
    Vy_i, Vy_min, Vy_max = get_itp(df, :Vy)
    P_i, P_min, P_max = get_itp(df, :P)
    return Vx_i, Vx_min, Vx_max, Vy_i, Vy_min, Vy_max, P_i, P_min, P_max
end


Vx_i, Vx_min, Vx_max, Vy_i, Vy_min, Vy_max, P_i, P_min, P_max = interpolate_star("StarCCM/Output/test export 2.csv")
n = 2^9;
xs = LinRange(-0.5, 1.5, n);
ys = LinRange(-0.5, 0.5, n);

@time contour(xs, ys, (x,y) -> evaluate(Vx_i, [x,y])[1]; fill=(true,), clims=(Vx_min, Vx_max))



Vx_i,Vx, Vx_min, Vx_max, Vy_i, Vy, Vy_min, Vy_max, P_i, P, P_min, P_max = plotUVP("StarCCM/Output/test export 2.csv"; n)


contour(xs, ys, Vx'; fill=(true,), clims=(Vx_min, Vx_max))
contour(xs, ys, Vy'; fill=(true,), clims=(Vy_min, Vy_max))
contour(xs, ys, P'; fill=(true,), clims=(P_min, P_max))
foil = Airfoil("Airfoils/naca/naca6622.dat")

VVx = deepcopy(Vx);
VVy = deepcopy(Vy);
PP = deepcopy(P);

for i in 1:n, j in 1:n
    if (foil.leftmost + 1/n) <= xs[i] <= (foil.rightmost - 1/n)
        vu = findvertdist(foil.upperspline, [xs[i], ys[j]])
        vl = findvertdist(foil.lowerspline, [xs[i], ys[j]])
        if vu < (0 - 1/n) && vl > (0 + 1/n)
            VVx[i,j] = -10^8
            VVy[i,j] = -10^8
            PP[i,j] = -10^8
        end
    end
end

contour(xs, ys, VVx'; fill=(true, :jet), clims=(Vx_min, Vx_max))
contour(xs, ys, VVy'; fill=(true, :jet), clims=(Vy_min, Vy_max))
contour(xs, ys, PP'; fill=(true, :jet), clims=(P_min, P_max))

savefig("Paper/figures/P_n128.png")


function plot_training_data_6(test_num)
    @load "training_data/sdf_puv_params_128x128x241_Mk1.jld2" SDF_UV PUV Naca AoA Re

    SDF = SDF_UV[:, :, 1, test_num]
    Ui = SDF_UV[:, :, 2, test_num]
    Vi = SDF_UV[:, :, 3, test_num]

    P = PUV[:, :, 1, test_num]
    U = PUV[:, :, 2, test_num]
    V = PUV[:, :, 3, test_num]

    cs = contour(SDF'; size=(1000,800), fill=(true, :jet), grid=false, xaxis=false, yaxis=false, title="SDF")
    cui = contour(Ui'; size=(1000,800), fill=(true, :jet), grid=false, xaxis=false, yaxis=false, title="Initial U Velocity")
    cvi = contour(Vi'; size=(1000,800), fill=(true, :jet), grid=false, xaxis=false, yaxis=false,title="Initial V Velocity")
    cp = contour(P'; size=(1000,800), fill=(true, :jet), grid=false, xaxis=false, yaxis=false, title="Pressure")
    cu = contour(U'; size=(1000,800), fill=(true, :jet), grid=false, xaxis=false, yaxis=false, title="U Velocity")
    cv = contour(V'; size=(1000,800), fill=(true, :jet), grid=false, xaxis=false, yaxis=false, title="V Velocity")

    plot(cs, cp, cui, cu, cvi, cv; layout=(3,2))
end

struct X_SDF
    sdf::Array{Float32, 4}
    mask::Array{Bool, 4}
    re::Vector{Float32}
    aoa::Vector{Float32}
    u_ini::Vector{Float32}
    v_ini::Vector{Float32}
end

function convert_to_trainable()
    @load "training_data/sdf_puv_params_128x128x526, zeroed.jld2" SDF_UV PUV Re AoA;

    shuff = shuffle(1:size(PUV, 4));
    SDF_UV = SDF_UV[:,:,:,shuff];
    PUV = PUV[:,:,:,shuff];
    Re = Re[shuff]
    AoA = AoA[shuff]

    sdf::Array{Float32, 4} = SDF_UV[:,:,1:1,:];

    mask::Array{Bool, 4} = sdf .> 0;

    μ = 1.85508e-5 # Pa-s
    ρ = 1.18415 # kg/m^3

    u::Vector{Float32} = Re*μ/ρ .* cos.(AoA*pi/180)
    v::Vector{Float32} = Re*μ/ρ .* sin.(AoA*pi/180)

    re = Re .|> Float32
    aoa = AoA .|> Float32

    x_sdf = (sdf, mask, re, aoa, u, v)
    y_puv = PUV;

    # x_sdf = X_SDF(sdf, mask, re, aoa, u, v)

    jldsave("training_data/airfoil_training_data.jld2"; x_sdf, y_puv)
end

function convert_to_small()
    @load "training_data/sdf_puv_params_128x128x526, zeroed.jld2" SDF_UV PUV Re AoA;

    sdf_full = SDF_UV[1:2:end,1:2:end,1:1,:];
    puv_full = PUV[1:2:end,1:2:end,1:3,:];
    sdf_small = zeros(Float32,64,64,1,0)
    puv_small = zeros(Float32,64,64,3,0)
    re_small = zeros(Float32,0)
    aoa_small = zeros(Float32,0)
    
    for i in axes(Re, 1)
        if 20000 < Re[i] < 80000 && -16 < AoA[i] < 16
            sdf_small = [sdf_small ;;;; sdf_full[:,:,:,i]]
            puv_small = [puv_small ;;;; puv_full[:,:,:,i]]
            re_small = [re_small ; Re[i]]
            aoa_small = [aoa_small ; AoA[i]]
        end
    end
    
    mask::Array{Bool, 4} = sdf_small .> 0;
    
    μ = 1.85508e-5 # Pa-s
    ρ = 1.18415 # kg/m^3
    
    u::Vector{Float32} = re_small*μ/ρ .* cos.(aoa_small*pi/180)
    v::Vector{Float32} = re_small*μ/ρ .* sin.(aoa_small*pi/180)
    
    re = re_small .|> Float32
    aoa = aoa_small .|> Float32
    
    x_sdf = (sdf_small, mask, re, aoa, u, v)
    y_puv = puv_small;
    
    # x_sdf = X_SDF(sdf, mask, re, aoa, u, v)
    
    jldsave("training_data/airfoil_training_data_small_aoa.jld2"; x_sdf, y_puv)
end
convert_to_small()

