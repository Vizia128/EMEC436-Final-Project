using LinearAlgebra, Plots, Interpolations

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


"""
    Airfoil(address::String)

Finds the upper and lower splines of an airfoil from a dat file. 

Assumptions: The points range from x=1 to x=0 for the upper section and 
then range from x=0 to x=1 on the lower section.
"""
function Airfoil(address::String)
    file = open(address)
    lines = readlines(file)
    name = popfirst!(lines)
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

    upperspline = Interpolations.scale(interpolate(upperpoints, 
        (BSpline(Cubic(Natural(OnGrid()))), NoInterp())), 1:size(upperpoints, 1), 1:2)
    lowerspline = Interpolations.scale(interpolate(lowerpoints, 
        (BSpline(Cubic(Natural(OnGrid()))), NoInterp())), 1:size(lowerpoints, 1), 1:2)

    leftmost = 0.0
    rightmost = 1.0

    return Airfoil(name, points, upperpoints, lowerpoints, upperspline, lowerspline, 
        leftmost, rightmost)
end

function Airfoil(name::String, upperpoints::Matrix, lowerpoints::Matrix)
    leftmost = lowerpoints[1,1]
    rightmost = upperpoints[1,1]

    upperspline = Interpolations.scale(interpolate(upperpoints, 
        (BSpline(Cubic(Natural(OnGrid()))), NoInterp())), 1:size(upperpoints, 1), 1:2)
    lowerspline = Interpolations.scale(interpolate(lowerpoints, 
        (BSpline(Cubic(Natural(OnGrid()))), NoInterp())), 1:size(lowerpoints, 1), 1:2)

    return Airfoil(name, [upperpoints; lowerpoints], upperpoints, lowerpoints, 
        upperspline, lowerspline, leftmost, rightmost)
end

"Rotates a 2D point about a pivot point."
function rotate(pivot::Vector, point::Vector, θ::Real)
    s, c = sin(θ), cos(θ)
    point -= pivot
    point = [c -s; s c] * point
    point += pivot
    return point
end

function findclosest2points(points, pos)
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

    pivot[2] = foil.upperpoints[ufb, 2] - foil.lowerpoints[lfb, 2]

    for (i, uf) in enumerate(eachrow(upperbody))
        upperbody[i,:] = rotate(pivot, upperbody[i,:], angle)
    end
    for (i, lf) in enumerate(eachrow(lowerflap))
        lowerflap[i,:] = rotate(pivot, lowerflap[i,:], angle)
    end
    return Airfoil(foil.name, [upperbody; upperflap], [lowerbody; lowerflap])
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
function findvertdist(spline, pos::Vector, maxdist, maxiter)
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
            vdist = findvertdist(foil.upperspline, pos, maxdist, maxiter)
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
            vdist = findvertdist(foil.lowerspline, pos, maxdist, maxiter)
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

function testrun(;n=256, m=256)
    foil = Airfoil("Airfoils/naca/naca6621.dat")
    foil = setflap_angle(foil, [0.8, 0], -π/12)

    # plot(foil.upperspline[:, 1], 
    #     foil.upperspline[:,2])
    # plot!(foil.lowerspline[:, 1], 
    #     foil.lowerspline[:,2])
    S = sdf_array(foil; n, m)

    contour(S'; size=(1000,800), fill=(true))
    plot!((foil.points[:,1].+0.5)*n/2, (foil.points[:,2].+1)*m/2; color=:white, linewidth=2)
end