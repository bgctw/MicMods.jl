"""
    LabeledParSetter 

Translates between a vector of optimization parameters
and the two vectors of a system: its parameters and initial states.
"""
struct LabeledParSetter{PAR,ST,POPT,SOPT, NP,NS,NPO, NSO} 
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
    # mapping symbol to thing indexing a solution
    num::OrderedCollections.OrderedDict{Symbol, Any}
end

""" 
    LabeledParSetter(
        systemt, 
        paropt = keys(systemt.searchranges_p), 
        stateopt = keys(systemt.searchranges_u0))

Construct a Parameter-Setter based on the subset of parameters and initial states.        
"""        
function LabeledParSetter(systemt::NamedTuple; 
    paropt = Tuple(Symbol.(keys(systemt.searchranges_p))), 
    stateopt = Tuple(Symbol(x.val.f) for x in keys(systemt.searchranges_u0))
    )
    LabeledParSetter(systemt.syss, paropt, stateopt)
end

function LabeledParSetter(system::ModelingToolkit.AbstractSystem, paropt, stateopt)
    parameters_map = OrderedDict(Symbol(x.name) => x for x in parameters(system))
    states_map = OrderedDict(Symbol(x.f.name) => x for x in states(system))
    observed_map = OrderedDict(Symbol(x.lhs.f.name) => x.lhs for x in observed(system))
    LabeledParSetter(
        Tuple(keys(parameters_map)), Tuple(keys(states_map)),
        paropt,stateopt; 
        num = merge(parameters_map, states_map, observed_map)
        )
end

""" 
    LabeledParSetter(parsys,statesys,paropt,stateopt)

Construct a Parameter-Setter based on the subset of parameters and initial states.        

## values        
All arguments are of type AbstractVector{Num}
- `parsys,statesys`: parameters and states of the system
- `paropt,stateopt`: parameters and states in optimization vector
"""        
# function LabeledParSetter(parsys::AbstractVector{Symbol},statesys::AbstractVector{Symbol}, 
#     paropt::AbstractVector{Symbol}, stateopt::AbstractVector{Symbol}; num = OrderedDict{Symbol, Any}())

function LabeledParSetter(parsys::NTuple{NP,Symbol}, statesys::NTuple{NS,Symbol}, 
    paropt, stateopt; num = OrderedDict{Symbol, Any}()) where {NP, NS}
    
    ispopt = falses(NP)   
    ipopt_for_ipsys = zeros(NP)
    iuopt_for_iusys = zeros(NS)
    isuopt = falses(NS)
    ipsys_for_ipopt = zeros(Int8, length(paropt))
    iusys_for_iuopt = zeros(Int8, length(stateopt))
    #i,sym = first(enumerate(paropt))
    for (ipopt,sym) = enumerate(paropt)
        ipsys = findfirst(isequal(sym), parsys)
        isnothing(ipsys)  && error("parameter $sym not found in system")
        ispopt[ipsys] = true
        ipopt_for_ipsys[ipsys] = ipopt
        ipsys_for_ipopt[ipopt] = ipsys
    end         
    #i,sym = first(enumerate(stateopt))
    for (iuopt,sym) = enumerate(stateopt)
        iusys = findfirst(isequal(sym), statesys)
        # may be given a symbols instead of num
        if isnothing(iusys) 
            iusys = findfirst(isequal(sym), map(x -> x.val.f.name, statesys))
        end
        isnothing(iusys) && error("initial value $sym not found in system")
        isuopt[iusys] = true
        iuopt_for_iusys[iusys] = iuopt
        iusys_for_iuopt[iuopt] = iusys
    end         
    # assume parameter vector vcat(p,u): position in entire parameter vector
    iopt_for_iusys = iuopt_for_iusys .+ length(paropt)
    #
    LabeledParSetter{parsys,statesys,paropt,stateopt,NP, NS,length(paropt),length(stateopt)}(
        ispopt, ipopt_for_ipsys, ipsys_for_ipopt,
        isuopt, iuopt_for_iusys, iusys_for_iuopt,
        iopt_for_iusys, num)
end

getpoptnames(ps::LabeledParSetter{PAR,ST,POPT,SOPT}) where {PAR,ST,POPT,SOPT} = POPT
getparsys(ps::LabeledParSetter{PAR,ST,POPT,SOPT}) where {PAR,ST,POPT,SOPT} = PAR
getstatesys(ps::LabeledParSetter{PAR,ST,POPT,SOPT}) where {PAR,ST,POPT,SOPT} = ST
getpopt(ps::LabeledParSetter{PAR,ST,POPT,SOPT}) where {PAR,ST,POPT,SOPT} = POPT
getsopt(ps::LabeledParSetter{PAR,ST,POPT,SOPT}) where {PAR,ST,POPT,SOPT} = SOPT


"""
    setpu(ps::LabeledParSetter, popt, prob)
    setpu(ps::LabeledParSetter, popt, p, u0)

Set values in tuple (prob.p, prob.u0) to popt.

Eltype of new p and u0 should correspond to type of popt, while vector
type should correspond to original p and u0. If `p` is of type `SVector{3,Float64}` and
`eltype(popt) == Dual128`, then new `p` is of type `SVector{3,Dual128}`.

Value: Tuple (p, u0)
"""
function setpu(ps::LabeledParSetter, popt, prob, ::Val{Labeled} = Val(false)) where Labeled
    setpu(ps, popt, prob.p, prob.u0, Val(Labeled))
end
# function setpu(ps::LabeledParSetter, popt, p, u0)
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

function setpu(ps::LabeledParSetter, popt, p, u0, ::Val{Labeled} = Val(true)) where Labeled
    # superseded: eltype of new p aund u0 should correspond to eltype of popt
    # type will be widened between popt and p and u0
    # try to make this type-stable - requires to make copies of p and u0
    #    
    # tmp = map(x -> x * zero(eltype(popt)), p)
    # tmp2 = (ps.ispopt[i] ? 
    # popt[ps.ipopt_for_ipsys[i]] : p[i] for i in axes(p,1))
    # tmp3 = convert(typeof(tmp), tmp .+ tmp)::typeof(tmp) # Vector{Float64}
    # pnew = map(*(zero(eltype(popt))), p) .+ (ps.ispopt[i] ? 
    #          popt[ps.ipopt_for_ipsys[i]] : p[i] for i in axes(p,1)) 
    p0 = map(x -> x * zero(eltype(popt)), p)
    pg = (ps.ispopt[i] ? popt[ps.ipopt_for_ipsys[i]] : p[i] for i in axes(p,1))
    pnew = convert(typeof(p0), p0 .+ pg)::typeof(p0)
    # u0new = u0 .* zero(eltype(u0)) .+ (ps.isuopt[i] ? 
    #         #note usage of popt[iopt] , indexing into the full optim vector with u after p
    #         popt[ps.iopt_for_iusys[i]] : u0[i] for i in axes(u0,1))
    u00 = map(x -> x * zero(eltype(popt)), u0)
    u0g = (ps.isuopt[i] ? 
            #note usage of popt[iopt] , indexing into the full optim vector with u after p
            popt[ps.iopt_for_iusys[i]] : u0[i] for i in axes(u0,1))
    u0new = convert(typeof(u00), u00 .+ u0g)::typeof(u00)
    if Labeled 
        (
            label_parsys(ps, pnew),
            label_statesys(ps, u0new),
        )
    else
      (pnew, u0new)
    end
end

# function setpu_labeled(ps::LabeledParSetter, popt, p, u0)
#     # eltype of new p aund u0 should correspond to eltype of popt
#     pnew, u0new = setpu(ps, popt, p, u0)
#         (
#             label_parsys(ps, pnew),
#             label_statesys(ps, u0new),
#         )
# end


"""
    label_parsys(ps::LabeledParSetter, popt)

Create an LVector from parameter vector
"""
function label_parsys(ps::LabeledParSetter{PAR,ST}, p) where {PAR, ST}
    #LVector(NamedTuple{PAR}(p)) # not typestable for vector
    # @LArray reuses the underlying storage 
    # this allows changing the underlying unnamed vector by names
    @LArray p PAR
    # ans = @LVector eltype(p) PAR
    # ans .= p
    # ans
    #NamedTuple{PAR}(p) # not type stable for p Vector
    #(; (PAR .=> p)...) # not type stable for p Vector
end

"""
    label_parsys(ps::LabeledParSetter, popt)

Create an LVector from parameter vector
"""
function label_statesys(ps::LabeledParSetter{PAR,ST}, u) where {PAR, ST}
    #LVector(NamedTuple{ST}(u))
    # ans = @LVector eltype(u) ST
    # ans .= u
    # ans
    @LArray u ST
end


"""
    getpopt(ps::LabeledParSetter, prob)

Extract the parameter vector from prob.p and prob.u0.

Value: Numeric array.
"""
function getpopt(ps::LabeledParSetter, prob, ::Val{Labeled} = Val(true)) where Labeled
    getpopt(ps, prob.p, prob.u0, Val(Labeled))
end
function getpopt(ps::LabeledParSetter{PAR,ST,POPT,SOPT}, p, u0, ::Val{Labeled} = Val(true)) where {PAR,ST,POPT,SOPT,Labeled}
    if Labeled 
        vcat(
            LVector(NamedTuple{POPT}(
                p[ps.ipsys_for_ipopt])), 
            LVector(NamedTuple{SOPT}(
                    u0[ps.iusys_for_iuopt]))
            )
    else
        vcat(p[ps.ipsys_for_ipopt], u0[ps.iusys_for_iuopt])
    end
end
function getpopt(ps::LabeledParSetter{PAR,ST,POPT,SOPT}, p::StaticVector, u0::StaticVector, ::Val{Labeled} = Val(true)) where {PAR,ST,POPT,SOPT,Labeled}
    if Labeled 
        vcat(
            SLVector(NamedTuple{POPT}(
                p[ps.ipsys_for_ipopt])), 
            SLVector(NamedTuple{SOPT}(
                    u0[ps.iusys_for_iuopt]))
            )
    else
        vcat(p[ps.ipsys_for_ipopt], u0[ps.iusys_for_iuopt])
    end
end


# """
#     label_popt(ps::LabeledParSetter, popt)

# Create an LVector from parameter vector
# """
# function label_popt(ps::LabeledParSetter, popt)
#     vcat(
#         LVector(NamedTuple{map(x -> x.val.name, ps.paropt).data}(
#             popt[1:length(ps.paropt)])), 
#         LVector(NamedTuple{map(x -> x.val.f.name, ps.stateopt).data}(
#             popt[length(ps.paropt)+1:end]))
#         )
# end

# """
#     parindex(ps::LabeledParSetter, sym::Symbol)

# Get the index corresponding to symbol in vector of system parameters.

# Value: Integer or nothing
# """
# function parindex(ps::LabeledParSetter, sym::Symbol)
#     i = findfirst(map(p -> p.val.name, ps.parsys) .== sym)
# end

# """
#     stateindex(ps::LabeledParSetter, sym::Symbol)

# Get the index corresponding to symbol in vector of system states.

# Value: Integer or nothing
# """
# function stateindex(ps::LabeledParSetter, sym::Symbol)
#     i = findfirst(map(u -> u.val.f.name, ps.statesys) .== sym)
# end




