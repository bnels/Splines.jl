# Splines.jl
A B-Spline interpolation package for Julia

## Getting Started
You can install the Splines module in Julia
~~~julia
Pkg.checkout("git@github.com:bnels/Splines.jl.git")
Pkg.build("Splines")
~~~
You can update the package using the usual
~~~julia
Pkg.update()
~~~
Now you're ready to start using Splines!


## Basics
Splines.jl uses B-Splines as a basis for constructing Spline interpolations.  This is all under the hood, so for basic spline manipulations, you only need to provide a knot sequence, function values at knots, and what order of spline you would like to use (e.g. 4th order splines are piecewise cubic).
~~~julia
using Splines

ts = [linspace(-10,10,50);] # knot sequence
vs = cos(t) # signal values at knots
m = 4 # cubic splines

S = Spline(vs, ts, m)
~~~
You can now treat your spline as a function, and can evaluate it at any points, or a vector of points
~~~julia
using PyPlot

xs = [linspace(-10,10,150);]
ys = S(xs)
plot(xs, ys)
~~~
![cos(x) spline example](./doc/figs/cos_spline.png)

## Examples

## Hilbert Transforms
