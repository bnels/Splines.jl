module Splines

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
function get_knot(x::Float64, s::Spline)
    
end


# De Boor's Algorithm for evaluating Spline at a point x
function SplineEval(s::Spline, ord::Int, x::Float64)


end



end #module
