function chak21_problem(syss, x0, parms; ti=range(0,48,step = 1/4))
    prob = ODEProblem{false}(syss, x0, extrema(ti), parms, saveat=ti, jac=true)
    prob = ODEProblem(syss, x0, extrema(ti), parms, saveat=ti, jac=true)
end
# function chak21_problem_static(syss, x0, ::Val{NU}, parms, Val{NP}; ti=range(0,48,step = 1/4)) where {NU, NP}
#   prob = chak21_problem(syss, x0, parms; ti)
# end

indexof(sym,syms) = findfirst(isequal(sym),syms)

"""
    constructparsetter(system, paropt = keys(system.searchranges_p), stateopt = keys(system.searchranges_u0))

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
- `paropt`: a generator of `Num`s corresponding to parameters optimized
- `stateopt`: a generator of `Num`s corresponding to states optimized

# Examples
```jldoctest am; output = false, setup = :()
true
# output
true
```
"""
# deprecated, superseded by struct ParSetter in parsetter.jl
# function constructparsetter(system::ModelingToolkit.AbstractSystem, paropt, stateopt)
#     parsys = convert.(Num,parameters(system))::Vector{Num}
#     statesys = convert.(Num,states(system))::Vector{Num}
#     ipar = collect(1:length(parsys))
#     ispopt = falses(length(parsys))   
#     iu = collect(1:length(statesys))
#     isuopt = falses(length(statesys))
#     posparopt = zeros(Int8, length(paropt))
#     posstateopt = zeros(Int8, length(stateopt))
#     #i,num = first(enumerate(paropt))
#     for (i,num) = enumerate(paropt)
#         pos = findfirst(isequal(num), parsys)
#         ispopt[pos] = true
#         ipar[pos] = i
#         posparopt[i] = pos
#     end         
#     #i,num = first(enumerate(stateopt))
#     for (i,num) = enumerate(stateopt)
#         pos = findfirst(isequal(num), statesys)
#         isuopt[pos] = true
#         iu[pos] = length(paropt) + i
#         posstateopt[i] = pos
#     end         
#     function setpu(popt, prob)
#         ([ispopt[i] ? popt[ipar[i]] : prob.p[ipar[i]] for i in axes(prob.p,1)],
#         [isuopt[i] ? popt[iu[i]] : prob.u0[iu[i]] for i in axes(prob.u0,1)])
#     end
#     function getpopt(prob)
#         vcat(prob.p[posparopt], prob.u0[posstateopt])
#     end
#     states_map = OrderedDict(Symbol(x.val.f.name) => x for x in statesys)
#     parameters_map = OrderedDict(Symbol(x.val.name) => x for x in parsys)
#     observed_map = OrderedDict(Symbol(x.lhs.f.name) => x.lhs for x in observed(system))
#     num = merge(states_map, parameters_map, observed_map)
#     (setpu = setpu, getpopt = getpopt, posparopt = posparopt, posstateopt = posstateopt, 
#     num = num)
# end

# function constructparsetter(systemt::NamedTuple, 
#     paropt = keys(systemt.searchranges_p), 
#     stateopt = keys(systemt.searchranges_u0)
#     )
#     constructparsetter(systemt.syss, paropt, stateopt)
# end


include("mod_simp.jl")
include("mod_phys.jl")
include("mod_growth.jl")
include("mod_growth_closed.jl")
include("mod_phys_fixedr.jl")
include("fit_phases_lim.jl")
include("fit_phases_growth.jl")
