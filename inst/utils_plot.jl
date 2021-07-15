"""
    plotchains(getchn, vars = keys(getchn))

Density plot of get(chain, [vars])
"""
function plotchain(getchn, vars = keys(getchn); kwargs...)
    map(vars) do var 
        x = getchn[var]
        density(x, ylab="density", xlab = string(var),
            label = label = permutedims(["chain " * string(i) for i in axes(x,2)]);
            kwargs...
            )
    end
end

"""
    plot_post_and_priors(systemt, ps, chn)

## Arguments
- systemt: Tuple with fields `systemt.searchranges_p`, and `systemt.searchranges_u0`    
- ps: Parameter setter with field `paropt` and `stateopt`
- chn: MCMCChains.Chains object
"""
function plot_post_and_priors(systemt, ps, chn)
    #using StatsPlots
    #@parameters ks km kd Y  # use ps.num[:symbol]
    srk = vcat(
        [key => systemt.searchranges_p[key] for key in ps.paropt],
        [key => systemt.searchranges_u0[key] for key in ps.stateopt]
    )
    #i=2; key, r = collect(srk)[i]
    plts = map(enumerate(srk)) do (i,(key, r))
        #(i,(key, r)) = p
        #r = systemt.searchranges_p[ps.num[:ks]]
        d = gettruncdist(r...)
        s = first(get(chn,[Symbol(MicMods.replacetby0(key))]))
        #s = Array(chn[:,i,:])
        xlim = MicMods.extendrange(quantile(vec(s), (0.025,0.975))..., 0.5)
        #plot(d, label = key)
        plot(d, label = "prior", xlim = xlim, ylab="density", xlab=string(key))
        density!(s, label = permutedims(["chain " * string(i) for i in axes(s,2)]))
    end
    # put all in one plot
    plot(plts..., layout = length(plts))
end

