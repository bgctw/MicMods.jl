module MicMods

using Base: thread_notifiers
using ModelingToolkit, OrderedCollections
using Distributions, StatsFuns
using MCMCChains
using StaticArrays, LabelledArrays
using Loess, Polynomials

export chak21_simp_system, chak21_phys_system, chak21_growth_system,
    chak21_growth_closed, chak21_fixedr_system,
    indexof, ParSetter, getpopt, setpu, parindex, stateindex, 
    label_parsys, label_statesys, label_popt,
    gettruncdist,
    kinresp_mic, kinresp_exp,
    micfromcoef,
    find_inflection

include("util.jl")
include("findinflection.jl")
include("parsetter.jl")
include("distfit.jl")
include("chakrawal21/base.jl")
include("kinresp/base.jl")

# Write your package code here.

end
