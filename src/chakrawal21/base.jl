function chak21_problem(syss, x0, parms; ti=range(0,48,step = 1/4))
    prob = ODEProblem(syss, x0, extrema(ti), parms, saveat=ti, jac=true)
end

indexof(sym,syms) = findfirst(isequal(sym),syms)

"""
    unum(sym::Symbol, system)
    pnum(sym::Symbol, system)

Get the `Num` associated with given symbol of state or parameter respectively.

# Examples
```jldoctest am; output = false, setup = :()
true
# output
true
```
"""
unum(sym::Symbol,system) = states(system.syss)[system.statespos[sym]],
pnum(sym::Symbol,system) = parameters(system.syss)[system.parmspos[sym]]

"""
    constructparsetter(system, paramnums = keys(system.searchranges_p), statenums = keys(system.searchranges_u0))

Optimization parameter vector includes a subset of model parameters and a subset 
of initial conditions. This function helps to translate between a vector of
parameters to optimize and vectors of parameters and initial states.

It returns a tuple of two functions
- `setpu(popt, prob)`: returns a tuple(p,u0) of corresponding vectors of Problem `prob` 
  overridden by values of `popt` and converted to eltype of `popt`.
- `getpopt(prob)`: returns a vector of parameters to optimized extracted from fields 
  p and u0 of Problem `prob`.

# Arguments
- `system`: tuple with Dicitonary fields `parmspos` and `statespos` that map 
  symbols to positions in the vectors.
- `paramnums`: a generator of `Num`s corresponding to parameters optimized
- `statenums`: a generator of `Num`s corresponding to states optimized

# Examples
```jldoctest am; output = false, setup = :()
true
# output
true
```
"""
function constructparsetter(system, paramnums = keys(system.searchranges_p), statenums = keys(system.searchranges_u0))
    ppos = system.parmspos
    ipar = collect(1:length(parameters(system.syss)))
    ispopt = falses(length(ppos))   
    upos = system.statespos
    iu = collect(1:length(states(system.syss)))
    isuopt = falses(length(upos))
    i,pname = first(enumerate(paramnums))
    for (i,pname) = enumerate(paramnums)
        pos = ppos[Symbol(pname)]
        ispopt[pos] = true
        ipar[pos] = i
    end         
    i,uname = first(enumerate(statenums))
    for (i,uname) = enumerate(statenums)
        pos = upos[Symbol(uname.val.f.name)]
        isuopt[pos] = true
        iu[pos] = length(paramnums) + i
    end         
    function setpu(popt, prob)
        ([ispopt[i] ? popt[ipar[i]] : prob.p[ipar[i]] for i in axes(prob.p,1)],
        [isuopt[i] ? popt[iu[i]] : prob.u0[iu[i]] for i in axes(prob.u0,1)])
    end
    function getpopt(prob)
        vcat(
            prob.p[getindex.(Ref(system.parmspos), (num.val.name for num in paramnums))],
            prob.u0[getindex.(Ref(system.statespos), (num.val.f.name for num in statenums))]
            )
    end
    poptnums = [paramnums..., statenums...]
    setpu, getpopt, poptnums
end

include("mod_simp.jl")
include("mod_phys.jl")
include("mod_growth.jl")

