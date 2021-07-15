# need to source these plotting functions, 
# because we do not want to import them in MicMods
# moved to inst/utils_plot.jl


"""
    replace "(t)" in sting vector by "_0"

To create parameter names, that MCMCChains.get can handle.
`get(chn, Symbol("b(t)")` is not extracted.
"""
function replacetby0(s)
    replace.(string.(s),"(t)" => "_0") 
end

"""
    extendrange(lower,upper, frac = 0.1)

Modify lower and upper to so that the included range is
extended by frac.

This is useful for settting plotting limits based on quantiles.
"""    
function extendrange(lower,upper, frac = 0.1)
    ext = (upper - lower)*frac/2
    lower = lower - ext
    upper = upper + ext
    (lower,upper)
end




function best_estimate(chn::MCMCChains.Chains)
    logden = first(get(chn, :log_density))
    imin = argmax(logden)
    chn[imin[1],:, imin[2]]
    vec(Array(chn[imin[1],:, imin[2]]))
end


