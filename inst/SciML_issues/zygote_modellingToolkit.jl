using ModelingToolkit, OrdinaryDiffEq
using DiffEqBase

@parameters α β δ γ
@variables t x(t) y(t) dx(t)
D = Differential(t)
eqs = [
  dx ~ α*x - β*x*y,  # testing observed variables
  D(x) ~ dx,
  D(y) ~ -δ*y + γ*x*y
]
@named lv = ODESystem(eqs)
syss = structural_simplify(lv) 
parms = [α => 1.5, β => 1.0, δ => 3.0, γ => 1.0]
x0 = [x => 1.0, y => 1.0]
tsteps = 0.0:0.1:10.0
prob = ODEProblem(syss, x0, extrema(tsteps), parms, jac = true)
soltrue = solve(prob,  Tsit5(), saveat = tsteps);
popt0 = [1.1]

using ChainRulesCore
function ChainRulesCore.rrule(::typeof(getindex), VA::ODESolution, sym) 
  function ODESolution_getindex_pullback(Δ)
    @show Δ
    @show length(VA)
    @show VA
    @show VA.u
    # convert symbol to index
    i = issymbollike(sym) ? sym_to_index(sym, VA) : sym
    @show i
    # similar to VectorOfArray: return zero for non-matching indices
    Δ′ = [ (i == j ? Δ : zero(x)) for (x,j) in zip(VA.u, 1:length(VA))]
    (NO_FIELDS, Δ′)
    # TODO: care for observed
  end  
  VA[sym], ODESolution_getindex_pullback(Δ)
end

f1(p) = soltrue[x][1] * p[1] # note the indexing by [x]
f1(popt0)
#using Zygote
gr = Zygote.gradient(f1, popt0) # calls the failing rule for VectorOfArrays instead of above rule

