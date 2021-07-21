function fit_initial_growth(tinfl, popt0g; systemtg = chak21_phys_system(), chak21syn = chak21syn, solver = Rodas5())
    ilimg = findall(chak21syn.solsyn.t .< tinfl)
    tlimgo = chak21syn.solsyn.t[ilimg]
    tlimg = vcat(tlimgo, tinfl)
    obsgr_tot = chak21syn.obsr_tot[ilimg]
    obsgq = chak21syn.obsq[ilimg]    
    
    psg = LabeledParSetter(systemtg, paropt=(), stateopt=(:b,:r))
    probp = remake(systemtg.prob, tspan = extrema(tlimgo))
    flossg = function(popt::AbstractVector{T}, psg=psg, par0 = systemtg.prob.p, u00 = systemtg.prob.u0) where T
        local pp, u0p = setpu(psg, popt, par0, u00)
        #pp, u0p = setpu(psg, popt, systemtg.prob)
        local solpg = solve(systemtg.prob, solver, u0 = u0p.__x, p = pp.__x, saveat=tlimgo);
        if solpg.retcode != :Success 
            predgq = zeros(T, length(obsgq))
            predgr_tot = zeros(T, length(obsgr_tot))
            lossval = convert(T,Inf)::T
        else
            predgq = solpg[psg.num[:q]]::Vector{T}
            # lossval = sum(abs2, predq .- obsq)
            predgr_tot = solpg[psg.num[:r_tot]]::Vector{T}
            # predgr_tot = predg[psg.num[:r_tot]];
            lossval = sum(abs2, predgr_tot - obsgr_tot)
        end
        return (lossval=lossval, sol=solpg, predq=predgq, predr_tot=predgr_tot)
    end
    #return flossg

    #using Optim 
    flossg1(p) = first(flossg(p))
    ansopt = optimize(flossg1, popt0g.__x, LBFGS(); autodiff = :forward)
    popto = label_popt(psg, Optim.minimizer(ansopt))
    return (popt = popto, popt0 = popt0g, floss = flossg, ps=psg, prob = probp, ansopt = ansopt)
end

function fit_growth_and_lim()
end