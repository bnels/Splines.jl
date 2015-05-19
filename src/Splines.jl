module Splines
export BasisSpline, BasisEval, SplineCoeffMatrix, SplineEvalMatrix #BSplineEval,

# TODO: implement more than zero BC


# Definition of collection of basis splines
type BasisSpline
    # The knot sequence vector on which this collection is to be defined
    t::Vector{Float64}
    # Order of splines (note: linear = 2 (NOT 1), quadratic = 3, etc.)
    m::Int

    # Redefine standard constructor
   BasisSpline(t::Vector{Float64}, m::Int=4) = new(sort(t), m)
end




# Evaluate the k-th basis function at point X
# We will not do this recursively
# d is the order of the derivative to evaluate
function BasisEval(B::BasisSpline, k::Int, x::Float64, d::Int=0)
    # Pad on m-1 knots to the beginning / end
    t_pad = [ B.t[1] * ones(B.m-1); B.t; B.t[end] * ones(B.m-1) ]
    # k relative to padded vector
    k_pad = k + B.m #(note that k starts at one)

    # Tree values
    T = zeros(B.m, B.m)

    # Check if x is in the support of each leaf-level basis function
    for i = 1:B.m
        T[i, 1] =  (t_pad[k_pad+i-1] <= x < t_pad[k_pad+i]) ? 1 : 0
    end

    # Recurse via iteration
    for j = 1:B.m-1
        for i = 1:B.m-j
            r1 = (t_pad[k_pad+i+j-1] - t_pad[k_pad+i-1])
            r2 = (t_pad[k_pad+i+j] - t_pad[k_pad+i])
            if B.m-1 - j >= d
                c1 = (r1 == 0) ? 0 : (x - t_pad[k_pad+i-1]) / r1
                c2 = (r2 == 0) ? 0 : (t_pad[k_pad+i+j] - x) / r2
            else
                c1 = (r1 == 0) ? 0 : j / r1
                c2 = (r2 == 0) ? 0 : -j / r2
            end

            T[i, j+1] = c1 * T[i, j] + c2 * T[i+1, j]
        end
    end
    return T[1,B.m]
end

# Evaluate for a whole bunch of points
function BasisEval(B::BasisSpline, k::Int, x::Vector{Float64}, d::Int=0)
    n = size(x,1)
    f = zeros(n)
    for i = 1:n
        f[i] = BasisEval(B, k, x[i], d)
    end
    return f
end

# # d is the order of the derivative to be evaluated
# function DerivEval(B::BasisSpline, k::Int, x::Vector{Float64}, d::Int)
#     if d == 0
#         return BasisEval(B, k, x)
#     end
# end

# Definition of Interpolatory spline
type Spline{T}

    # Function values sampled at the knot sequence
    v::Vector{T}
    B::BasisSpline

    # Redefine standard constructor
    # TODO: don't hold values, hold coefficients!
    function Spline(v::Vector{T}, B::BasisSpline, )
        if( length(B.t) == length(v) )
            new(v, B)
        else
            throw(DimensionMismatch())
        end
    end
end

# Utility constructor
Spline{T}(v::Vector{T}, t::Vector{Float64}, m::Int=4) = Spline(v, BasisSpline(t,m))

# Construct system matrix to interpolate with, assuming zero boundary conditions (i.e., spline gets flat)
# TODO: Derivative matching
# HACK: currently we eval everything slightly to the left of the knots because it isn't defined at the other location
function SplineCoeffMatrix(B::BasisSpline)
    # number of points to eval at
    m = size(B.t, 1)
    # Why is n so weird?
    n = B.m + length(B.t) - 2
    M = n-m
    # Why is offset so weird?
    offset = B.m

    A = zeros(n,n)
    for i = 1:m-1
        for j=1:n
            A[i,j] = BasisEval(B, j-offset, B.t[i])
        end
    end
    A[m,n] = 1.

    # Now we add derivative conditions at the ends
    for i = m+1:2:n
        d = Int(ceil( float(i-m)/ 2.))
        for j = 1:n
            A[i,j] = BasisEval(B, j-offset, B.t[1], d)
            if i+1 <= n
                A[i+1,j] = BasisEval(B, j-offset, B.t[m]-eps(B.t[m]), d)
            end
        end
    end
    return A
end

# Construct system matrix to interpolate with, assuming zero boundary conditions (i.e., spline gets flat)
# TODO: Derivative matching
function SplineEvalMatrix(B::BasisSpline, x::Vector{Float64})

    # number of points to eval at
    m = size(x, 1)
    # Why is n so weird?
    n = B.m + length(B.t) - 2
    # Why is offset so weird?
    offset = B.m

    A = zeros(m,n)
    for i = 1:m
        for j=1:n
            A[i,j] = BasisEval(B, j-offset, x[i])
        end
    end
    return A
end









# # Spline definition
# type Spline{T}
#     t::Vector{Float64}
#     u::Vector{T}
#     function Spline(t::Vector{Float64}, u::Vector{T})
#         if( length(t) == length(u) )
#             p = sortperm(t) # make sure we're sorted in time
#             t = t[p]
#             u = u[p]
#             new(t,u)
#         else
#             throw(DimensionMismatch())
#         end
#     end
# end

# # constructors
# Spline{T}(t::Vector{Float64}, u::Vector{T}) = Spline{T}(copy(t), copy(u))
# Spline{T}(t::Vector{Int}, u::Vector{T}) = Spline{T}(float(copy(t)), copy(u))
# Spline{T}(t::Range, u::Vector{T}) = Spline{T}(Array(float(t)), copy(u))
# Spline{T}(s::Spline{T}) = Spline{T}(copy(s.t), copy(s.u))

# # find which knot a point x lies in
# # returns first first index of s.t that is > x
# # if x is > than all of s.t, returns max index + 1
# #function get_knot(x::Float64, s::Spline)
# #    i = 1
# #    imax = length(s.t)
# #    found = false
# #    while !found
# #        x > s.t[i] || found = true
# #        found || i += 1 # only increment if we go to the next iteration
# #        i <= imax || found = true # quit if we've gone too far
# #    end
# #    return i
# #end

# # Evaluate BSpline Basis defined by knot sequence t at point x
# function BSplineEval(t::Vector{Float64}, ord::Int, knot::Int, x::Float64)
#     app  = t[end] * ones(ord+1)
#     bapp = t[1] * ones(ord+1)
#     knot = knot + ord + 1
#     tapp = [bapp; t; app]
#     B = zeros(ord+1, ord+1)

#     # Check if x is in the support of this basis function
#     for i = 1:ord+1
#         B[i, 1] =  (tapp[knot+i-1] <= x < tapp[knot+i]) ? 1 : 0
#     end

#     for k = 1:ord
#         for i = 1:1+ord-k
#             r1 = (tapp[knot+i+k-1] - tapp[knot+i-1])
#             r2 = (tapp[knot+i+k] - tapp[knot+i])
#             c1 = (r1 == 0) ? 0 : (x - tapp[knot+i-1]) / r1
#             c2 = (r2 == 0) ? 0 : (tapp[knot+i+k] - x) / r2
#             B[i, k+1] = c1 * B[i, k] + c2 * B[i+1, k]
#         end
#     end
#     return B[1,ord+1]
# end

# # Evaluate the BS

# function BSplineEval(t::Vector{Float64}, ord::Int, knot::Int, x::Vector{Float64})
#     n = size(x,1)
#     f = zeros(n)
#     for i = 1:n
#         f[i] = BSplineEval(t, ord, knot, x[i])
#     end
#     return f
# end

# create matrix that has BSpline basis function value evaluated at each knot
# function SplineCoeffMatrix(t::Vector{Float64}, ord::Int)
#     # m = size(t,1)
#     # n = m + ord - 1
#     # A = zeros(n,m)
#     # for i = 1:m-1
#     #     for j = 1:n-1
#     #         A[j,i] = BSplineEval(t, ord, j-ord, t[i])
#     #     end
#     # end
#     # A[n,m] = 1.0
#     # return A


#     # Number of equations: 1 for each knot plus zero padd
#     # Number of knots without repeats
#     n_knot = size(t,1)

#     # Number of basis functions is total number of knots (including repeats)
#     n = n_knot + 2*(ord-2)

#     # number of equations is the same, we will use n_knot of them to be interpolation conditions
#     # and m - n_knot to set derivatives to zero at the end

#     m = n

#     A = zeros(m,n)
#     # for i = 1:n_knot-1
#     #     for j = 1:n-1
#     #         A[i,j] = BSplineEval(t, ord, i-ord, t[j])
#     #     end
#     # end
#     return A
# end

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
