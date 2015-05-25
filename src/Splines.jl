module Splines
export BasisSpline, Spline, BasisEval, SplineCoeffMatrix, SplineEvalMatrix, PadKnots

# TODO: implement more than zero BC


# Definition of collection of basis splines
type BasisSpline
    # The knot sequence vector on which this collection is to be defined
    t::Vector{Float64}
    # length of knot sequence
    n::Int
    # Order of splines (note: linear = 2 (NOT 1), quadratic = 3, etc.)
    m::Int

    # Redefine standard constructor
   function BasisSpline(t::Vector{Float64}, m::Int=4)
        # sort and pad the knot sequence.  Default padding is repeated knots
        sort!(t)
        difft = diff(t)
        t = t + 0.5*[ difft; difft[end]; ]
        n = length(t)
        t_pad = [t[1] * ones(m-1); t; t[end]*ones(m-1)]
        return new(t_pad, n, m)
   end
end

function PadKnots(B::BasisSpline, scheme::ASCIIString="repeat")
    if scheme == "reflect"
        for i = 1:(B.m - 1)
            B.t[i] = 2.0*B.t[B.m] - B.t[2*B.m - 1 - i]
            B.t[B.m + B.n - 1 + i] = 2.0*B.t[B.m + B.n - 1] + B.t[B.m + B.n - 1 - i]
        end
    elseif scheme == "periodic"
        for i = 1:(B.m - 1)
            B.t[i] = B.t[B.n + i]
            B.t[B.n + B.m - 1 + i] = B.t[B.m - 1 + i]
        end
    elseif scheme == "repeat"
        for i = 1:(B.m - 1)
            B.t[i] = B.t[B.m]
            B.t[B.n + B.m - 1 + i] = B.t[B.n + B.m - 1]
        end
    elseif scheme == "extend"
        k = knots(B)

        delta_start = k[2] - k[1]
        delta_end   = k[end] - k[end-1]

        for i = 1:B.m-1
            B.t[B.m-i] = k[1] - delta_start * i
            B.t[end - B.m + i + 1] = k[end] + delta_end * i
        end

    end

end


# return original knot sequence (unpadded)
function knots(B::BasisSpline)
    n = size(B.t,1)
    return sub(B.t, (B.m):n-B.m+1)
end


# Evaluate the k-th basis function at point X
# We will not do this recursively
# d is the order of the derivative to evaluate
function BasisEval(B::BasisSpline, k::Int, x::Float64, d::Int=0, hilbert::Bool=false)
    # Pad on m-1 knots to the beginning / end
    t_pad = B.t
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
            T[i, 1] = log(abs((x-t_pad[k_pad+i-1]) / (x-t_pad[k_pad+i])))/π
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
function SplineCoeffMatrix(B::BasisSpline)
    # number of points to eval at
    t = knots(B)
    n_interior_knots = B.n
    # Why is n so weird?
    n_basis_functions = B.m + n_interior_knots - 2

    A = zeros(n_interior_knots, n_interior_knots)

    #offset = B.m
    idxs = [Int(j) for j in -B.m/2:n_interior_knots-B.m/2-1]

    # TODO: don't evaluate outside of support of spline

    for i = 1:n_interior_knots
        lb = max(i-B.m, 1)
        ub = min(i+B.m, n_interior_knots)
        for j = lb:ub
            A[i,j] = BasisEval(B, idxs[j], t[i])
        end
    end

    return sparse(A)
end

# Construct system matrix to interpolate with, assuming zero boundary conditions (i.e., spline gets flat)
function SplineEvalMatrix(B::BasisSpline, x::Vector{Float64}, derivs::Int=0, hilbert::Bool=false)

    n_interior_knots = B.n
    # number of points to eval at
    m = size(x, 1)
    A = zeros(m,n_interior_knots)

    idxs = [Int(j) for j in -B.m/2:n_interior_knots-B.m/2-1]

    for i = 1:m
        for j = 1:n_interior_knots
            A[i,j] = BasisEval(B, idxs[j], x[i], derivs, hilbert)
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
    # TODO: add flag everywhere to allow changing pad type
    PadKnots(B, "extend")
    if( B.n == length(v) )
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

function call{T}(S::Spline{T}, x::Vector{Float64}, derivs::Int=0, hilbert::Bool=false)
    A = SplineEvalMatrix(S.B, x, derivs, hilbert)
    return A * S.alpha
end


function call{T}(S::Spline{T}, x::Float64, derivs::Int=0, hilbert::Bool=false)
    return S([x], derivs, hilbert)[1]
end


function *( scalar::Number, S::Spline )
    alpha1 = [ promote(i,scalar)[1] for i in S.alpha ]
    T = typeof(alpha1[1])
    return Spline{T}( S.B, scalar * S.alpha )
end

function *( S::Spline, scalar::Number )
    return scalar * S
end

# function +( S::Spline, scalar::Number )
#     return Spline(S.)
# end

# function +( S::Spline, scalar::Number )
#     return S + scalar
# end

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
