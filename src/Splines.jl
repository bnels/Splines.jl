module Splines
export BasisSpline, Spline, BasisEval, SplineCoeffMatrix, SplineEvalMatrix #BSplineEval,

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
function BasisEval(B::BasisSpline, k::Int, x::Float64, d::Int=0, hilbert::Bool=false)
    # Pad on m-1 knots to the beginning / end
    t_pad = [ B.t[1] * ones(B.m-1); B.t; B.t[end] * ones(B.m-1) ]
    # k relative to padded vector
    k_pad = k + B.m #(note that k starts at one)

    # Tree values
    T = zeros(B.m, B.m)

    # Check if x is in the support of each leaf-level basis function
    if !hilbert
        for i = 1:B.m
            T[i, 1] =  (t_pad[k_pad+i-1] <= x < t_pad[k_pad+i]) ? 1 : 0
        end
    else #What do I do here if it goes to infinity?
        for i = 1:B.m
            T[i, 1] = log(abs((x-t_pad[k_pad+i-1]) / (x-t_pad[k_pad+i])))
        end
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
function BasisEval(B::BasisSpline, k::Int, x::Vector{Float64}, d::Int=0, hilbert::Bool=false)
    n = size(x,1)
    f = zeros(n)
    for i = 1:n
        f[i] = BasisEval(B, k, x[i], d, hilbert)
    end
    return f
end





# Construct system matrix to interpolate with, assuming zero boundary conditions (i.e., spline gets flat)
# TODO: Derivative matching, make sparse
# HACK: currently we eval the last guy slightly to the left of the knots because it isn't defined at the other location
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
    return sparse(A)
end

# Construct system matrix to interpolate with, assuming zero boundary conditions (i.e., spline gets flat)
function SplineEvalMatrix(B::BasisSpline, x::Vector{Float64}, hilbert::Bool)

    # number of points to eval at
    m = size(x, 1)
    # Why is n so weird?
    n = B.m + length(B.t) - 2
    # Why is offset so weird?
    offset = B.m

    A = zeros(m,n)
    for i = 1:m
        for j=1:n
            A[i,j] = BasisEval(B, j-offset, x[i], 0, hilbert)
        end
    end
    return sparse(A)
end




# Definition of Interpolatory spline
type Spline{T}
    # Function values sampled at the knot sequence
    B::BasisSpline
    alpha::Vector{T}

    Spline{T} (B::BasisSpline, alpha::Vector{T}) = new(B, alpha)
end

# Construct a set of B spline coefficients from values
function Spline{T}(v::Vector{T}, B::BasisSpline)
    if( length(B.t) == length(v) )
        A = SplineCoeffMatrix(B)
        m_mat = size(A,1)
        v1 = [ v; zeros(m_mat - length(v)); ]
        alpha = A \ v1
        return Spline{T}(B, alpha)
    else
        throw(DimensionMismatch())
    end
end

# Utility constructor
Spline{T}(v::Vector{T}, t::Vector{Float64}, m::Int=4) = Spline(v, BasisSpline(t,m))

function call{T}(S::Spline{T}, x::Vector{Float64}, hilbert::Bool=false)
    A = SplineEvalMatrix(S.B, x, hilbert)
    return A * S.alpha
end


function call{T}(S::Spline{T}, x::Float64, hilbert::Bool=false)
    return S([x], hilbert)[1]
end

function *( scalar::Number, S::Spline )
    alpha1 = [ promote(i,scalar)[1] for i in S.alpha ]
    T = typeof(alpha1[1])
    return Spline{T}( S.B, scalar * S.alpha )
end

function *( S::Spline, scalar::Number )
    return scalar * S
end

function +( S1::Spline, S2::Spline )
    knots = unique(sort([ S1.B.t; S2.B.t ]))
    vals  = [ S1(knots[i]) + S2(knots[i]) for i in 1:length(knots) ]
    m = max( S1.B.m, S2.B.m )
    return Spline(vals, knots, m)
end

function -( S1::Spline, S2::Spline )
    return S1 + ((-1.) * S2)
end

function /( S::Spline, scalar::Number )
    return (1.0/scalar) * S
end


end #module
