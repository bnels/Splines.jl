module Splines

export BSplineEval

# Spline definition
type Spline{T}
    t::Vector{Float64}
    u::Vector{T}
    function Spline(t::Vector{Float64}, u::Vector{T})
        if( length(t) == length(u) )
            p = sortperm(t) # make sure we're sorted in time
            t = t[p]
            u = u[p]
            new(t,u)
        else
            throw(DimensionMismatch())
        end    
    end
end

# constructors
Spline{T}(t::Vector{Float64}, u::Vector{T}) = Spline{T}(copy(t), copy(u))
Spline{T}(t::Vector{Int}, u::Vector{T}) = Spline{T}(float(copy(t)), copy(u))
Spline{T}(t::Range, u::Vector{T}) = Spline{T}(Array(float(t)), copy(u))
Spline{T}(s::Spline{T}) = Spline{T}(copy(s.t), copy(s.u))

# find which knot a point x lies in
# returns first first index of s.t that is > x
# if x is > than all of s.t, returns max index + 1
#function get_knot(x::Float64, s::Spline)
#    i = 1
#    imax = length(s.t)
#    found = false
#    while !found
#        x > s.t[i] || found = true
#        found || i += 1 # only increment if we go to the next iteration
#        i <= imax || found = true # quit if we've gone too far
#    end
#    return i
#end

# Evaluate BSpline Basis defined by knot sequence t at point x
function BSplineEval(t::Vector{Float64}, ord::Int, knot::Int, x::Float64)
    B = zeros(ord+1, ord+1)
    for i = 1:ord+1
        (t[knot + i - 1] <= x < t[knot + i]) ? B[i,1] = 1 : B[i,1] = 0
    end
    for k = 1:ord
        for i = 1:1+ord-k
            r1 = (t[knot + i + k - 1] - t[knot + i - 1])
            r2 = (t[knot + i + k] - t[knot + i])
            r1 == 0 ? c1 = 0 : c1 = (x - t[knot + i - 1]) / r1
            r2 == 0 ? c2 = 0 : c2 = (t[knot + i + k] - x) / r2
            B[i,k+1] = c1*B[i,k] + c2*B[i+1,k]
        end
    end
    return B[1,ord+1]
end

function BSplineEval(t::Vector{Float64}, ord::Int, knot::Int, x::Vector{Float64})
    n = size(x,1)
    f = zeros(n)
    for i = 1:n
        f[i] = BSplineEval(t, ord, knot, x[i])
    end
    return f
end

# create matrix that has BSpline basis function value evaluated at each knot
function SplineMatrix(t::Vector{Float64}, ord::Int)
end

#function DeBoor(s::Spline, ord::Int, x::Float64, k::Int, i::Int)
#    α = (x - s.u[i]) / (s.u[i + ord + 1 - k] - s.u[i])
#    return (1 - α) * DeBoor(s, ord, x, k-1, i-1) + α * DeBoor(s, ord, x, k-1, i)
#end

# De Boor's Algorithm for evaluating Spline at a point x
#function SplineEval(s::Spline, ord::Int, x::Float64)
#    knot = get_knot(x, s)
#    return DeBoor(s, ord, x, ord, knot)
#end


end #module
