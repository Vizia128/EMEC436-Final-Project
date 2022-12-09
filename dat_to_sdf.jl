using LinearAlgebra, Plots, Interpolations, CSV, DataFrames, JLD2

"Rotates a 2D point about a pivot point."
function rotate(pivot::Vector, point::Vector, θ::Real)
    s, c = sin(θ), cos(θ)
    point -= pivot
    point = [c -s; s c] * point
    point += pivot
    return point
end

struct Airfoil
    name::String
    points::Matrix
    upperpoints::Matrix
    lowerpoints::Matrix
    upperspline
    lowerspline
    leftmost
    rightmost
end

function save_as_csv(name, xy)
    z = zeros(size(xy, 1))
    xyz = [xy z]
    
    CSV.write("Airfoils/CSV_naca/$name.csv", DataFrame(xyz, :auto); writeheader=false)
end

"""
    Airfoil(address::String)

Finds the upper and lower splines of an airfoil from a dat file. 

Assumptions: The points range from x=1 to x=0 for the upper section and 
then range from x=0 to x=1 on the lower section.
"""
function Airfoil(address::String)
    file = open(address)
    lines = readlines(file)
    name = first(popfirst!(lines), 8)
    slines = strip.(lines, [' '])
    splitlines = split.(slines, " "; limit = 2)
    numpnts = size(splitlines)[1]

    points = zeros(size(splitlines, 1), 2)

    upperpoints = 0
    for i in 1:numpnts
        sp = splitlines[i]
        points[i, :] = parse.(Float64, splitlines[i])
        if points[i, 1] == 0.0
            upperpoints = points[1:i, :]
        end
    end
    if upperpoints == 0
        error("could not find leading edge")
    end
    lowerpoints = points[size(upperpoints, 1):end, :]

    upperspline = Interpolations.scale(Interpolations.interpolate(upperpoints, 
        (BSpline(Cubic(Natural(OnGrid()))), NoInterp())), 1:size(upperpoints, 1), 1:2)
    lowerspline = Interpolations.scale(Interpolations.interpolate(lowerpoints, 
        (BSpline(Cubic(Natural(OnGrid()))), NoInterp())), 1:size(lowerpoints, 1), 1:2)

    leftmost = 0.0
    rightmost = 1.0

    save_as_csv(name, points)

    return Airfoil(name, points, upperpoints, lowerpoints, upperspline, lowerspline, 
        leftmost, rightmost)
end

function Airfoil(name::String, upperpoints::Matrix, lowerpoints::Matrix)
    leftmost = lowerpoints[1,1]
    rightmost = upperpoints[1,1]

    upperspline = Interpolations.scale(Interpolations.interpolate(upperpoints, 
        (BSpline(Cubic(Natural(OnGrid()))), NoInterp())), 1:size(upperpoints, 1), 1:2)
    lowerspline = Interpolations.scale(Interpolations.interpolate(lowerpoints, 
        (BSpline(Cubic(Natural(OnGrid()))), NoInterp())), 1:size(lowerpoints, 1), 1:2)

    return Airfoil(name, [upperpoints; lowerpoints], upperpoints, lowerpoints, 
        upperspline, lowerspline, leftmost, rightmost)
end

function Airfoil(address::String, angle::Real)
    foil = Airfoil(address)

    for i in size(foil.upperpoints, 1)
        foil.upperpoints[i,:] = rotate([0,0], foil.upperpoints[i,:], angle)
    end
    for i in size(foil.lowerpoints, 1)
        foil.lowerpoints[i,:] = rotate([0,0], foil.lowerpoints[i,:], angle)
    end
    return Airfoil(foil.name, foil.lowerpoints, foil.upperpoints)
end

function findclosest2points(points, pos::Vector)
    pnt1index::Int = 0
    pnt2index::Int = 0
    pnt1dist = 1.0e8
    pnt2dist = 1.0e8
    for (i, point) in enumerate(eachrow(points))
        tempdist = norm(point - pos, 2)
        if tempdist <= pnt1dist
            pnt1index, pnt2index = i, pnt1index
            pnt1dist, pnt2dist = tempdist, pnt1dist
        elseif tempdist < pnt2dist
            pnt2index = i
            pnt2dist = tempdist
        end
    end
    return [pnt1index pnt2index; pnt1dist pnt2dist]
end

function findclosest2points(pointsA::Matrix, pointsB::Matrix)
    # pntindex = [A1 B1; A2 B2]
    pntindexA = [0 0]
    pntindexB = [0 0]
    pntdist = 1e8*[1,1]

    for (a, pntA) in enumerate(eachrow(pointsA)), (b, pntB) in enumerate(eachrow(pointsB))
        tempdist = norm(pntA - pntB, 2)
        if tempdist < pntdist[1]
            pntdist[1], pntdist[2] = tempdist, pntdist[1]
            pntindexA[1], pntindexA[2] = a, pntindexA[1]
            pntindexB[1], pntindexB[2] = b, pntindexB[1]
        elseif tempdist < pntdist[2]
            pntdist[2] = tempdist
            pntindexA[2] = a
            pntindexB[2] = b
        end
    end

    if pointsA[pntindexA[1]] > pointsA[pntindexA[2]]
        pntindexA[1], pntindexA[2] = pntindexA[2], pntindexA[1]
    end
    if pointsB[pntindexB[1]] > pointsB[pntindexB[2]]
        pntindexB[1], pntindexB[2] = pntindexB[2], pntindexB[1]
    end
    return pntindexA, pntindexB
end

function remove_overlap(body, flap)
    ibody, iflap = findclosest2points(body, flap)
    xbody = [body[ibody[1]], body[ibody[2]]]
    xflap = [flap[iflap[1]], flap[iflap[2]]]
    
    # Check if there is no overlap
    if xbody[2] < xflap[1]
        return body, flap
    end

    # Check direction of spline
    if ibody[1] < ibody[2]
        return body[1:ibody[1], :], flap[iflap[2]:end, :]
    end

    return body[ibody[2]:end, :], flap[1:iflap[1], :]
end

function setflap_angle(foil::Airfoil, pivot::Vector, angle::Real)
    upperbreak = findclosest2points(foil.upperpoints[:,1], [pivot[1]])
    lowerbreak = findclosest2points(foil.lowerpoints[:,1], [pivot[1]])

    if upperbreak[1,1] > upperbreak[1,2]
        upperbreak[:,1], upperbreak[:,2] = upperbreak[:,2], upperbreak[:,1]
    end
    if lowerbreak[1,1] < lowerbreak[1,2]
        lowerbreak[:,1], lowerbreak[:,2] = lowerbreak[:,2], lowerbreak[:,1]
    end

    ubb = Int(upperbreak[1,1])
    ufb = Int(upperbreak[1,2])
    lbb = Int(lowerbreak[1,2])
    lfb = Int(lowerbreak[1,1])

    upperbody = foil.upperpoints[1:ubb, :]
    upperflap = foil.upperpoints[ufb:end, :]
    lowerbody = foil.lowerpoints[1:lbb, :]
    lowerflap = foil.lowerpoints[lfb:end, :]

    pivot[2] = (foil.upperpoints[ufb, 2] + foil.lowerpoints[lfb, 2])/2

    for (i, uf) in enumerate(eachrow(upperbody))
        upperbody[i,:] = rotate(pivot, upperbody[i,:], angle)
    end
    for (i, lf) in enumerate(eachrow(lowerflap))
        lowerflap[i,:] = rotate(pivot, lowerflap[i,:], angle)
    end

    upperbody, upperflap = remove_overlap(upperbody, upperflap)
    lowerbody, lowerflap = remove_overlap(lowerbody, lowerflap)
    return Airfoil(foil.name, [upperbody; upperflap], [lowerbody; lowerflap])
end

function setflap_angle(foil::Airfoil, pivot::Real, angle)
    (;name, upperpoints, lowerpoints) = foil

    for (i, pnt) in enumerate(eachrow(upperpoints))
        if pnt[1] < pivot
            continue
        end
        upperpoints[i,2] += (upperpoints[i,1] - pivot)*sin(angle)
        upperpoints[i,1] -= (upperpoints[i,1] - pivot)*abs(sin(angle))
    end
    for (i, pnt) in enumerate(eachrow(lowerpoints))
        if pnt[1] < pivot
            continue
        end
        lowerpoints[i,2] += (lowerpoints[i,1] - pivot)*sin(angle)
        lowerpoints[i,1] -= (lowerpoints[i,1] - pivot)*abs(sin(angle))
    end
    return Airfoil(name, upperpoints, lowerpoints)
end

function findclosestpoint(spline, c2pnts::Matrix, pos::Vector, 
    maxdist::Float64, maxiter::Int)

    for iter in 1:maxiter
        midindex = (c2pnts[1,1] + c2pnts[1,2])/2
        midpoint = [spline(midindex, 1), spline(midindex, 2)]
        mid_dist = norm(midpoint - pos, 2)

        if mid_dist <= c2pnts[3]
            c2pnts[:, 1], c2pnts[:, 2] = [midindex, mid_dist], c2pnts[:, 1]
        elseif mid_dist <= c2pnts[4]
            c2pnts[:, 2] = [midindex, mid_dist]
        end

        if (c2pnts[2,1] - c2pnts[2,2]) < maxdist
            pp = [spline(c2pnts[1,1], 1), spline(c2pnts[1,1], 2)]
            dd = c2pnts[2,1]
            return pp, dd
        end 
    end
    return [spline(c2pnts[1,1], 1), spline(c2pnts[1,1], 2)], c2pnts[2,1]
end


"""
Finds the vertical distance between a point and a B-spline. 
This function assumes the point is within the horizontal domain of the spline.
"""
function findvertdist(spline, pos::Vector; maxdist::Float64 = 0.001, maxiter::Int = 100)
    t1 = 1
    t2 = size(spline, 1)
    tm = t2/2

    x1 = spline(t1, 1)
    x2 = spline(t2, 1)

    if x1 > x2
        t1, x1, t2, x2 = t2, x2, t1, x1
    end

    for iter in 1:maxiter
        tm = (t1 + t2)/2
        xm = spline(tm, 1)

        if abs(x1 - pos[1]) < maxdist
            return pos[2] - spline(t1, 2)

        elseif abs(x2 - pos[1]) < maxdist
            return pos[2] - spline(t2, 2)

        elseif abs(xm - pos[1]) < maxdist
            return pos[2] - spline(tm, 2)

        elseif x1 < pos[1] < xm
            t2, x2 = tm, xm
        elseif xm < pos[1] < x2
            t1, x1 = tm, xm
        else
            error("out of range! (t1, tm, t2) = ($t1, $tm, $t2) , 
                (x1, xm, x2) = ($x1, $xm, $x2) , pos[1] = $(pos[1])")
        end
    end
    error("Never converged! (t1, t2) = ($t1, $t2) , 
        (x1, x2) = ($x1, $x2) , pos[1] = $(pos[1])")
end

function sdf(foil::Airfoil, pos::Vector; maxdist::Float64 = 0.001, maxiter::Int = 100)
    u2pnts = findclosest2points(foil.upperpoints, pos)
    l2pnts = findclosest2points(foil.lowerpoints, pos)
    upper = l2pnts[2,1] + l2pnts[2,2] > u2pnts[2,1] + u2pnts[2,2] ? true : false

    if upper
        cpnt, cdist = findclosestpoint(foil.upperspline, u2pnts, pos, maxdist, maxiter)
        if foil.leftmost <= pos[1] <= foil.rightmost
            vdist = findvertdist(foil.upperspline, pos; maxdist, maxiter)
            if vdist >= 0
                return cdist
            else
                return -cdist
            end
        else
            return cdist
        end
    else
        cpnt, cdist = findclosestpoint(foil.lowerspline, l2pnts, pos, maxdist, maxiter)
        if foil.leftmost <= pos[1] <= foil.rightmost
            vdist = findvertdist(foil.lowerspline, pos; maxdist, maxiter)
            if vdist <= 0
                return cdist
            else
                return -cdist
            end
        else
            return cdist
        end
    end
end

function sdf_array(foil::Airfoil; n = 16, m = 16, 
    maxdist::Float64 = 1/(n+m), maxiter::Int = 10000)

    A = zeros(n,m)
    for i in 1:n, j in 1:m
        A[i,j] = sdf(foil, [2i/n-0.5, 2j/m-1]; maxdist, maxiter)
    end
    return A
end

function get_SDF_array(naca, AoA, Re, x_coordinates, y_coordinates;
    maxdist::Float64 = 1/512, maxiter::Int = 10000)

    dir = "Airfoils/naca/naca$naca.dat"
    airfoil = Airfoil(dir)

    SDF = zeros(size(x_coordinates, 1), size(y_coordinates, 1))

    for (i, x) in enumerate(x_coordinates), (j, y) in enumerate(y_coordinates)
        SDF[i,j] = sdf(airfoil, [x,y]; maxdist, maxiter)
    end

    return SDF
end

function testrun(;n=256, m=256)
    foil = Airfoil("Airfoils/naca/naca4411.dat")
    # foil = setflap_angle(foil, [0.8, 0], π/2)
    # foil = setflap_angle(foil, 0.7, π/2*(0.2))

    # plot(foil.upperspline[:, 1], 
    #     foil.upperspline[:,2])
    # plot!(foil.lowerspline[:, 1], 
    #     foil.lowerspline[:,2])
    S = sdf_array(foil; n, m)
    xs = LinRange(-0.5, 1.5, n)
    ys = LinRange(-1, 1, m)

    contour(xs, ys, S'; size=(1000,800), fill=(true, :jet))
    plot!((foil.points[:,1]), (foil.points[:,2]); color=:white, linewidth=2, label="NACA 6622")
    # savefig("Paper/figures/sdf_naca6622.png")
end


function get_parameters(dir)
    vars = split(dir, ",_")

    naca = split(vars[1], "naca")[2]
    naca = String(naca)

    AoA = split(vars[2], "=")[2]
    AoA = parse(Float64, AoA)

    Re = split(vars[3], "=")[2]
    Re = parse(Float64, Re)

    return naca, AoA, Re
end

function get_PUV_array(dir; inside_obj_num = 0)
    df = CSV.File(dir) |> DataFrame
    rename!(df, ["Velocity: Magnitude (m/s)", "Velocity[i] (m/s)", "Velocity[j] (m/s)", 
        "Pressure (Pa)", "Vorticity: Magnitude (/s)", "X (m)", "Y (m)", "Z (m)"] .=> 
        ["Vm", "Vx", "Vy", "P", "Vorticity", "x", "y", "z"])

    x_min = minimum(df[!, :x])
    x_max = maximum(df[!, :x])
    y_min = minimum(df[!, :y])
    y_max = maximum(df[!, :y])

    df[!, :x] = Int.((df[!, :x] .- x_min) / (x_max - x_min) * 128)
    df[!, :y] = Int.((df[!, :y] .- y_min) / (y_max - y_min) * 128)

    P = ones(129,129)*inside_obj_num
    U = ones(129,129)*inside_obj_num
    V = ones(129,129)*inside_obj_num

    for row in eachrow(df)
        i = row[:x] + 1
        j = row[:y] + 1
        global P[i,j] = row[:P]
        global U[i,j] = row[:Vx]
        global V[i,j] = row[:Vy]
    end

    return [P;;; U;;; V]
end

function get_data_arrays(;zero_value=0)
    did_not_converge = CSV.File("Output/_convergence.csv") |> DataFrame
    did_not_converge = replace.(did_not_converge, ",_res_03000" => "")
    did_not_converge = replace.(did_not_converge, "Residuals_" => "")
    did_not_converge = did_not_converge[!,1]

    output_dirs = readdir("Output/")

    PUV_complete = zeros(129, 129, 3, 0)
    SDF_complete = zeros(129, 129, 3, 0)

    naca_complete = []
    AoA_complete::Vector{Float64} = []
    Re_complete::Vector{Float64} = []

    x_coordinates = ones(129)*zero_value
    y_coordinates = ones(129)*zero_value

    for i in axes(x_coordinates, 1)
        x_coordinates[i] = (i - 1)/128*2 - 0.5
        y_coordinates[i] = (i - 1)/128 - 0.5
    end

    for dir in output_dirs
        if occursin("Residuals", dir) || occursin("_all_", dir) || 
            !iszero(occursin.(did_not_converge, dir)) || occursin("_convergence", dir)
            continue
        end
        naca, AoA, Re = get_parameters(dir)

        PUV = get_PUV_array("Output/$dir")
        PUV_complete = [PUV_complete ;;;; PUV]

        SDF = get_SDF_array(naca, AoA, Re, x_coordinates, y_coordinates)
        SDF_bit = SDF .> 0
        μ = 1.85508e-5 # Pa-s
        ρ = 1.18415 # kg/m^3
        
        U_ini = SDF_bit * Re*μ/ρ * cos(AoA*pi/180)
        V_ini = SDF_bit * Re*μ/ρ * sin(AoA*pi/180)
        SDF_UV = [SDF ;;; U_ini ;;; V_ini]

        SDF_complete = [SDF_complete ;;;; SDF_UV]

        append!(naca_complete, naca)
        append!(AoA_complete, AoA)
        append!(Re_complete, Re)

        println("Finished naca$naca, AoA = $AoA, Re = $Re")
    end

    return SDF_complete, PUV_complete, naca_complete, AoA_complete, Re_complete
end

function save_training_data()
    SDF_UV, PUV, Naca, AoA, Re = get_data_arrays()
    __, PUV_Inf, __, __, __ = get_data_arrays(zero_value=Inf)

    SDF_UV = Float32.(SDF_UV[1:end-1, 1:end-1, :, :])
    PUV = Float32.(PUV[1:end-1, 1:end-1, :, :])
    AoA = Float32.(AoA)
    Re = Float32.(Re)
    samples = size(Re,1)

    jldsave("training_data/sdf_puv_params_128x128x$samples, normed.jld2"; SDF_UV, PUV, PUV_Inf, Naca, AoA, Re)
end

# save_training_data()