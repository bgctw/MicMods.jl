"""
    ParSetter 

Translates between a vector of optimization parameters
and the two vectors of a system: its parameters and initial states.
"""
struct ParSetter{NP,NS,NPO,NSO} 
    # the keys of the parameters optimized
    paropt::SVector{NPO,Num}
    # the keys of the states optimized
    stateopt::SVector{NSO,Num}
    #
    # the keys of the parameters in the system
    parsys::SVector{NP,Num}
    # the keys of the states in the system
    statesys::SVector{NS,Num}
    #
    # is the position in system parameter vector optimized?
    ispopt::SVector{NP,Bool}
    # mapping index system-parameter -> index of optim-parameter
    ipopt_for_ipsys::SVector{NP,Int8}
    # mapping index of optim-parameter -> index of system-parameter
    ipsys_for_ipopt::SVector{NPO,Int8}
    #
    # similar for system states
    isuopt::SVector{NS,Bool}
    iuopt_for_iusys::SVector{NS,Int8}
    iusys_for_iuopt::SVector{NSO,Int8}
    # mapping index system-state -> index of optim-vector (state after parms)
    # = iuopt_for_iusys + length(paropt)
    iopt_for_iusys::SVector{NS,Int8}
    # 
    # mapping symbol to Num
    num::OrderedCollections.OrderedDict{Symbol, Num}
end

""" 
    ParSetter(
        systemt, 
        paropt = keys(systemt.searchranges_p), 
        stateopt = keys(systemt.searchranges_u0))

Construct a Parameter-Setter based on the subset of parameters and initial states.        
"""        
function ParSetter(systemt::NamedTuple; 
    paropt = keys(systemt.searchranges_p), 
    stateopt = collect(keys(systemt.searchranges_u0))
    )
    ParSetter(systemt.syss, paropt, stateopt)
end

function ParSetter(system::ModelingToolkit.AbstractSystem, paropt, stateopt)
    parsys = convert.(Num,parameters(system))::Vector{Num}
    statesys = convert.(Num,states(system))::Vector{Num}
    ParSetter(parsys,statesys,paropt,stateopt; observedvars=observed(system))
end

""" 
    ParSetter(parsys,statesys,paropt,stateopt)

Construct a Parameter-Setter based on the subset of parameters and initial states.        

## values        
All arguments are of type AbstractVector{Num}
- `parsys,statesys`: parameters and states of the system
- `paropt,stateopt`: parameters and states in optimization vector
"""        
function ParSetter(parsys::AbstractVector{Num},statesys::AbstractVector{Num}, 
    paropt, stateopt; observedvars=SVector{0,Num}())
    states_map = OrderedDict(Symbol(x.val.f.name) => x for x in statesys)
    parameters_map = OrderedDict(Symbol(x.val.name) => x for x in parsys)
    observed_map = OrderedDict(Symbol(x.lhs.f.name) => x.lhs for x in observedvars)
    num = merge(states_map, parameters_map, observed_map)
    #
    
    if eltype(paropt) <: Symbol 
        paropt = getindex.(Ref(parameters_map), paropt)
    else
        paropt = collect(paropt) #collect to proper vector, make a copy
    end
    if eltype(stateopt) <: Symbol 
        stateopt = getindex.(Ref(states_map), stateopt)
    else
        stateopt = collect(stateopt) #collect to proper vector, make a copy
    end
    ispopt = falses(length(parsys))   
    ipopt_for_ipsys = zeros(length(parsys))
    iuopt_for_iusys = zeros(length(statesys))
    isuopt = falses(length(statesys))
    ipsys_for_ipopt = zeros(Int8, length(paropt))
    iusys_for_iuopt = zeros(Int8, length(stateopt))
    #i,num = first(enumerate(paropt))
    for (ipopt,num) = enumerate(paropt)
        ipsys = findfirst(isequal(num), parsys)
        isnothing(ipsys)  && error("parameter $num not found in system")
        ispopt[ipsys] = true
        ipopt_for_ipsys[ipsys] = ipopt
        ipsys_for_ipopt[ipopt] = ipsys
    end         
    #i,num = first(enumerate(stateopt))
    for (iuopt,num) = enumerate(stateopt)
        iusys = findfirst(isequal(num), statesys)
        # may be given a symbols instead of num
        if isnothing(iusys) 
            iusys = findfirst(isequal(num), map(x -> x.val.f.name, statesys))
        end
        isnothing(iusys) && error("initial value $num not found in system")
        isuopt[iusys] = true
        iuopt_for_iusys[iusys] = iuopt
        iusys_for_iuopt[iuopt] = iusys
    end         
    # assume parameter vector vcat(p,u): position in entire parameter vector
    iopt_for_iusys = iuopt_for_iusys .+ length(paropt)
    #
    ParSetter{length(parsys),length(statesys),length(paropt),length(stateopt)}(
        paropt, stateopt, 
        parsys, statesys,
        ispopt, ipopt_for_ipsys, ipsys_for_ipopt,
        isuopt, iuopt_for_iusys, iusys_for_iuopt,
        iopt_for_iusys,
        num)
end

getpoptnames(ps::ParSetter) = vcat(
    map(x -> Symbol(x.val.name), ps.paropt), 
    map(x -> Symbol(x.val.f.name), ps.stateopt))



"""
    setpu(ps::ParSetter, popt, prob)
    setpu(ps::ParSetter, popt, p, u0)

Set values in tuple (prob.p, prob.u0) to popt.

Eltype of new p and u0 should correspond to type of popt, while vector
type should correspond to original p and u0. If `p` is of type `SVector{3,Float64}` and
`eltype(popt) == Dual128`, then new `p` is of type `SVector{3,Dual128}`.

Value: Tuple (p, u0)
"""
function setpu(ps::ParSetter, popt, prob, ::Val{Labeled} = Val(false)) where Labeled
    setpu(ps, popt, prob.p, prob.u0, Val(Labeled))
end
# function setpu(ps::ParSetter, popt, p, u0)
#     # return type of p and u0 instead of Vector, need constructor on generator
#     # generator only works in special cases
#     # pnew = typeof(p)( (ps.ispopt[i] ? 
#     #         popt[ps.ipopt_for_ipsys[i]] : p[i] for i in axes(p,1)) )
#     # # pnew::Vector{eltype(popt)} = [ps.ispopt[i] ? 
#     # #     popt[ps.ipopt_for_ipsys[i]] : p[i] for i in axes(p,1)]
#     # u0new = typeof(u0)( (ps.isuopt[i] ? 
#     #         #note usage of popt[iopt] , indexing into the full optim vector with u after p
#     #         popt[ps.iopt_for_iusys[i]] : u0[i] for i in axes(u0,1)) )
#     #     # u0new::Vector{eltype(popt)} = [ps.isuopt[i] ? 
#     #     # #note usage of popt[iopt] , indexing into the full optim vector with u after p
#     #     # popt[ps.iopt_for_iusys[i]] : u0[i] for i in axes(u0,1)]
#     pnew = p .* zero(eltype(p)) .+ (ps.ispopt[i] ? 
#              popt[ps.ipopt_for_ipsys[i]] : p[i] for i in axes(p,1)) 
#     u0new = u0 .* zero(eltype(u0)) .+ (ps.isuopt[i] ? 
#             #note usage of popt[iopt] , indexing into the full optim vector with u after p
#             popt[ps.iopt_for_iusys[i]] : u0[i] for i in axes(u0,1))
#     (pnew, u0new)
# end
function setpu(ps::ParSetter, popt, p, u0, ::Val{Labeled} = Val(false)) where Labeled
    # eltype of new p aund u0 should correspond to eltype of popt
    pnew = p .* zero(eltype(popt)) .+ (ps.ispopt[i] ? 
             popt[ps.ipopt_for_ipsys[i]] : p[i] for i in axes(p,1)) 
    u0new = u0 .* zero(eltype(popt)) .+ (ps.isuopt[i] ? 
            #note usage of popt[iopt] , indexing into the full optim vector with u after p
            popt[ps.iopt_for_iusys[i]] : u0[i] for i in axes(u0,1))
    !Labeled && return((pnew, u0new))
    (
        label_parsys(ps, pnew),
        label_statesys(ps, u0new),
    )
end

"""
    label_parsys(ps::ParSetter, popt)

Create an LVector from parameter vector
"""
function label_parsys(ps::ParSetter, p)
    LVector(NamedTuple{map(x -> x.val.name, ps.parsys).data}(p))
end

"""
    label_parsys(ps::ParSetter, popt)

Create an LVector from parameter vector
"""
function label_statesys(ps::ParSetter, u)
    LVector(NamedTuple{map(x -> x.val.f.name, ps.statesys).data}(u))
end



"""
    getpopt(ps::ParSetter, prob)

Extract the parameter vector from prob.p and prob.u0.

Value: Numeric array.
"""
function getpopt(ps::ParSetter, prob, ::Val{Labeled} = Val(false)) where Labeled
    getpopt(ps, prob.p, prob.u0, Val(Labeled))
end
function getpopt(ps::ParSetter, p, u0, ::Val{Labeled} = Val(false)) where Labeled
    if Labeled 
        vcat(
            LVector(NamedTuple{map(x -> x.val.name, ps.paropt).data}(
                p[ps.ipsys_for_ipopt])), 
            LVector(NamedTuple{map(x -> x.val.f.name, ps.stateopt).data}(
                    u0[ps.iusys_for_iuopt]))
            )
    else
        vcat(p[ps.ipsys_for_ipopt], u0[ps.iusys_for_iuopt])
    end
end

"""
    label_popt(ps::ParSetter, popt)

Create an LVector from parameter vector
"""
function label_popt(ps::ParSetter, popt)
    vcat(
        LVector(NamedTuple{map(x -> x.val.name, ps.paropt).data}(
            popt[1:length(ps.paropt)])), 
        LVector(NamedTuple{map(x -> x.val.f.name, ps.stateopt).data}(
            popt[length(ps.paropt)+1:end]))
        )
end

"""
    parindex(ps::ParSetter, sym::Symbol)

Get the index corresponding to symbol in vector of system parameters.

Value: Integer or nothing
"""
function parindex(ps::ParSetter, sym::Symbol)
    i = findfirst(map(p -> p.val.name, ps.parsys) .== sym)
end

"""
    stateindex(ps::ParSetter, sym::Symbol)

Get the index corresponding to symbol in vector of system states.

Value: Integer or nothing
"""
function stateindex(ps::ParSetter, sym::Symbol)
    i = findfirst(map(u -> u.val.f.name, ps.statesys) .== sym)
end




