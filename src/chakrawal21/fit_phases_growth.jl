function fit_initial_growth(tinfl, popt0g; systemtg = chak21_phys_system(), chak21syn = chak21syn, solver = Rodas5())
    ilimg = findall(chak21syn.solsyn.t .< tinfl)
    tlimgo = chak21syn.solsyn.t[ilimg]
    tlimg = vcat(tlimgo, tinfl)
    obsgr_tot = chak21syn.obsr_tot[ilimg]
    obsgq = chak21syn.obsq[ilimg]    
    
    psg = LabeledParSetter(systemtg, paropt=(), stateopt=(:b,:r))
    probpg = remake(systemtg.prob, tspan = (0, maximum(tlimgo)))
    flossg = function(popt::AbstractVector{T}, psg=psg, par0 = systemtg.prob.p, u00 = systemtg.prob.u0) where T
        local pp, u0p = setpu(psg, popt, par0, u00)
        #pp, u0p = setpu(psg, popt, systemtg.prob)
        local solpg = solve(probpg, solver, u0 = u0p.__x, p = pp.__x, saveat=tlimgo);
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
    return (popt = popto, popt0 = popt0g, floss = flossg, ps=psg, prob = probpg, ansopt = ansopt)
end

function fit_growth_and_lim(popt0f; systemtf = chak21_phys_system(), chak21syn = chak21syn, solver = Rodas5())
    tlimfo = chak21syn.solsyn.t
    obsfr_tot = chak21syn.obsr_tot
    obsfq = chak21syn.obsq   

    psf = LabeledParSetter(systemtf, paropt=(:ks, :km, :kd, :Y), stateopt=(:b,:r)) 
    probpf = remake(systemtf.prob, tspan = (0, maximum(tlimfo)))
    #popt = popt0f
    flossf = let probpf = probpf, psf=psf, tlimfo=tlimfo, obsfq=obsfq, obsfr_tot = obsfr_tot
        function(popt::AbstractVector{T}) where T
            local pp, u0p = setpu(psf, popt, probpf)
            #pp, u0p = setpu(psf, popt, systemtf.prob)
            local solpf = solve(probpf, solver, u0 = u0p.__x, p = pp.__x, saveat=tlimfo);
            if solpf.retcode != :Success 
                predfq = zeros(T, length(obsfq))
                predfr_tot = zeros(T, length(obsfr_tot))
                lossval = convert(T,Inf)::T
            else
                predfq = solpf[psf.num[:q]]::Vector{T}
                # lossval = sum(abs2, predq .- obsq)
                predfr_tot = solpf[psf.num[:r_tot]]::Vector{T}
                # predgr_tot = predg[psf.num[:r_tot]];
                lossval = sum(abs2, predfr_tot - obsfr_tot)
            end
            return (lossval=lossval, sol=solpf, predq=predfq, predr_tot=predfr_tot)
        end
    end
    return (floss = flossf, )

    #using Turing, MCMCChains
    sr = label_popt(psf, vcat(
        [systemtf.searchranges_p[psf.num[key]] for key in getpopt(psf)],
        [systemtf.searchranges_u0[psf.num[key]] for key in getsopt(psf)]
    ))
    # uses: sr, floss, σ_obsq, σ_obsr_tot
    σ_obsq = chak21syn.σ_obsq
    σ_obsr_tot = chak21syn.σ_obsr_tot
    @model function fitqf(qobs, r_totobs,::Type{T} = Float64) where {T}
        p = Vector{T}(undef, length(sr))
        
        for (i,r) = enumerate(sr)
             p[i] ~ gettruncdist(r...)
        end
        if !isa(_context, Turing.PriorContext)
            lossval, solp, predq, predr_tot = flossf(p);
            if !isfinite(lossval) 
                Turing.@addlogprob! -Inf; return
            end
            # for (i, qp) = enumerate(predq)
            #     qobsl[i] ~ Normal(qp, σ_obsq)
            # end
            for (i, rp) = enumerate(predr_tot)
                r_totobs[i] ~ Normal(rp, σ_obsr_tot)
            end
        end
    end
    model_func = fitqf(obsfq, obsfr_tot)
    #tmp = sample(model_func, MH(), 10)
    #tmp = @suppress_err sample(model_func, NUTS(), 10)
    chn = sample(model_func, NUTS(),  MCMCThreads(), 10, 3, 
        theta_init = repeat(popt0f.__x, 1, 3)
        );
    chn = replacenames(chn, Dict("p[$i]" => pname for (i,pname) in 
        enumerate(replacetby0(getoptnames(psf)))))
    #popt1 = MicMods.label_popt(psf, MicMods.best_estimate(chn))
    poptof = label_popt(psf, best_estimate(chn))
    return (popt = poptof, popt0 = popt0f, floss = flossf, ps=psf, prob = probpf, chn = chn)
end