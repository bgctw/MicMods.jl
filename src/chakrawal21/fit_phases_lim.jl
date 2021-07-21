function fit_initial_lim(tinfl, popt0; systemt = chak21_fixedr_system(), chak21syn = chak21syn, solver = Rodas5())
    iinfl = findfirst(chak21syn.solsyn.t .> tinfl) - 1

    ilim = (iinfl+1):length(chak21syn.solsyn.t)
    tlimo = chak21syn.solsyn.t[ilim]
    tlim = vcat(tinfl, tlimo)  # include inflection point
    obsr_tot = chak21syn.obsr_tot[ilim];
    obsq = chak21syn.obsq[ilim];
    #
    ps = LabeledParSetter(systemt)
    p0, u0 = setpu(ps, popt0, systemt.prob)
    #u0g, cr_tot  = MicMods.estimate_u0_from_cumresp(tinfl, popt0, systemt, chak21syn; ps)
    # u0g, cr_tot = estimate_u0_from_cumresp(tinfl, popt0, systemt, chak21syn; ps)
    # u0g, cr_tot = MicMods.estimate_u0_from_cumresp(tinfl, u0, u0.b, chak21syn.solsyn.t, chak21syn.obsr_tot)
    u0g, cr_tot = estimate_u0_from_cumresp(tinfl, u0, u0.b, chak21syn.solsyn.t, chak21syn.obsr_tot)
    !(sum(u0g) + cr_tot ≈ sum(u0)) && error("mass balance u at tinfl: $(sum(u0g) + cr_tot) != $(sum(u0))")
    s_plus_b = u0g.s + u0g.b
    #p0.b = u0g.b
    #probp = remake(systemt.prob; tspan = extrema(tlim), u0 = u0g.__x)
    probp = remake(systemt.prob; tspan = extrema(tlim))

    popt0.b = u0g.b # update popt0 initial b - changes popt also outside

    #p = getpopt(ps, systemt.prob, Val(true)) 
    # uses systemt, ps, probp, solver, obsq, obsr_tot, s_plus_b
    # for typestability use let (https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-captured)
    floss_lim = let probp = probp, systemt = systemt, ps=ps, solver=solver, obsq = obsq, obsr_tot = obsr_tot, s_plus_b=s_plus_b
        function(p::AbstractVector{T}, ps=ps, par0=probp.p, u00=probp.u0) where T
            # pp = vcat(p[1:4], systemt.prob.p[5:end])
            # #u0p = convert.(eltype(p),systemt.prob.u0)
            # u0p = vcat(systemt.prob.u0[1], p[5], systemt.prob.u0[3])
            #pp, u0p = ps.setpu(p, systemt.prob)
            # change the optimized subset of parameter and initial state
            local pp, u0p = setpu(ps, p, par0, u00, Val(true))
            # adjustr_tot initial substrate for estimated biomass, so that sum 
            # of bimass +substret is kept constant
            u0p.s = s_plus_b - u0p.b
            local solp = solve(probp, 
                solver, 
                maxiters = 1000, # stop early if cannot determine 
                # use unlabelled parameters to avoid many compilations
                p = pp.__x, u0 = u0p.__x,
                #isoutofdomain = (u, p, t) -> (any(ui < zero(ui) for ui in u)),
                saveat = tlim,
            );
            local lossval, predq, predr_tot
            if solp.retcode != :Success 
                predq = zeros(T, length(obsq))
                predr_tot = zeros(T, length(obsr_tot))
                lossval = convert(T,Inf)::T
            else
                # fixed problems in gradient but is not type-stable
                predq = solp[ps.num[:q]][2:end]::Vector{T}
    #            lossval = sum(abs2, predq .- obsq)
                predr_tot = solp[ps.num[:r_tot]][2:end]::Vector{T}
                # local prednew = systemt.predout.(solp.u, Ref(pp.__x))
                # predq = getindex.(prednew, :q)[2:end]
                # predr_tot = getindex.(prednew, :r_tot)[2:end]
                lossval = sum(abs2, predr_tot .- obsr_tot)
                # @show T
                # @show typeof(lossval)
            end
            return (lossval=lossval, sol=solp, predq=predq, predr_tot=predr_tot)
        end
    end
    # return (popt = popt0, popt0 = popt0, s_plus_b, floss = floss_lim, ps = ps, prob = probp)

    # #using Optim # fails to converge
    # floss_lim1(popt) = first(floss_lim(popt))
    # #return floss_lim1(popt0)
    # @show popt0
    # ansopt = optimize(floss_lim1, popt0, LBFGS(); autodiff = :forward)
    # popto = Optim.minimizer(ansopt)
    # return (popt = popto, floss = floss_lim, prob = probp, ansopt = ansopt)

    #using Turing, MCMCChains
    sr = label_popt(ps, vcat(
        [systemt.searchranges_p[ps.num[key]] for key in getpopt(ps)],
        [systemt.searchranges_u0[ps.num[key]] for key in getsopt(ps)]
    ))
    # uses: sr, floss, σ_obsq, σ_obsr_tot
    σ_obsq = chak21syn.σ_obsq
    σ_obsr_tot = chak21syn.σ_obsr_tot
    @model function fitql(qobsl, r_totobsl,::Type{T} = Float64) where {T}
        p = Vector{T}(undef, length(sr))
        for (i,r) = enumerate(sr)
             p[i] ~ gettruncdist(r...)
        end
        if !isa(_context, Turing.PriorContext)
            lossvall, solpl, predql, predr_totl = floss_lim(p);
            if !isfinite(lossvall) 
                Turing.@addlogprob! -Inf; return
            end
            # for (i, qp) = enumerate(predq)
            #     qobsl[i] ~ Normal(qp, σ_obsq)
            # end
            for (i, rp) = enumerate(predr_totl)
                r_totobsl[i] ~ Normal(rp, σ_obsr_tot)
            end
        end
    end
    model_func = fitql(obsq, obsr_tot)
    #tmp = sample(model_func, MH(), 10)
    #tmp = @suppress_err sample(model_func, NUTS(), 10)
    #Random.seed!(0815)
    # chn = chn0 = sample(model_func, NUTS(),  MCMCThreads(), 40, 3)
    # chn = replacenames(chn0, Dict("p[$i]" => pname for (i,pname) in 
    #      enumerate(MicMods.replacetby0(getpoptnames(ps)))
    # drop the 1st chain
    #chn = chn[:,:,[2,3]]
    # resample using 
    chn = chn1 = sample(model_func, NUTS(),  MCMCThreads(), 100, 3, 
    #chn = chn1 = @suppress_err sample(model_func, NUTS(),  MCMCThreads(), 400, 3, 
        #theta_init = MicMods.best_estimate(chn0)
        );
    chn = replacenames(chn1, Dict("p[$i]" => pname for (i,pname) in 
        enumerate(replacetby0(getoptnames(ps)))))
    popt1 = MicMods.label_popt(ps, MicMods.best_estimate(chn))

#        theta_init = popt_true);
    popto = label_popt(ps, best_estimate(chn1))
    return (popt = popto, popt0 = popt0, s_plus_b = s_plus_b, floss = floss_lim, ps=ps, prob = probp, chn = chn, obsr_tot=obsr_tot)
end

function tmpf_p28370()
    chn = replacenames(chn1, Dict("p[$i]" => pname for (i,pname) in 
          enumerate(MicMods.replacetby0(getpoptnames(ps)))))
          popt1 = MicMods.label_popt(ps, MicMods.best_estimate(chn1))

    tmp = Array(chn[:,:,[2,3]]);
    # tmp = Array(chn[:,:,[3]]);
    # tmp = Array(chn[:,:,[1]]);
    # tmp = Array(chn[:,:,[1,3]]);
    #chn = Chains(reshape(tmp, (size(tmp)...,1)), chn.name_map.parameters) 
    chn = Chains(reshape(tmp, (size(tmp)...,1)), collect(getpoptnames(ps))) 
    function interactive_plotchains()
        popt1 = MicMods.label_popt(ps, MicMods.best_estimate(chn1))
        plot(chn)
        corner(chn)
        map(display, plotchain(get(chn1, :log_density)));
        plot_post_and_priors(systemt, ps, chn)
        soltrue[ps.num[:s]]
    end
end

# function iplot_29348()
#     #+using StatsPlots
#     ans_optlim = fit_initial_lim(tinfl, popt0; systemt, chak21syn, solver);
#     popt1 = first(ans_optlim)

# end

# function tmpf_check_type_stability()
#     poptl0 = copy(popt0);  poptl0.b = 40.0 # start with higher biomass
#     floss_lim = fit_initial_lim(tinfl, popt0; systemt, chak21syn, solver);
#     tmp = floss_lim(poptl0); first(tmp)
#     floss1(p) = first(floss_lim(p))
#     floss1(poptl0)
#     #@code_warntype floss_simp(poptl0)

#     #using ForwardDiff
#     ForwardDiff.gradient(floss1, poptl0)
# end

"""
    estimate_u0_from_cumresp(tinfl, popt0, systemt, chak21syn; ps = LabelledParSetter(systemt))

Obtain initial estimate of substrate s and biomass u by computing uptake from microbial parameters and respiration.

# Arguments
- `tinfl`: (numeric) start-time, after amendment
- `u0`: labelled array of initial state with entries s and b
- `b0`: microbial biomass at amendmend time 0. The default is obtained from chak21syn problem
- `t`, `r_tot`: respiration rate measured at time t
- `r0`, `Y`: representative phyisological state and yield for the period 0..tinfl.

# Value
- Tuple `(u0, cr_tot_gr)` of parameters and cumulated respiration.
"""
function estimate_u0_from_cumresp(tinfl, u0, b0, t, r_tot, ::Val{replace_u0} = Val(false);
    r0 = 0.6, 
    Y0 = 0.7,
    )  where replace_u0
    # cumulative respiration: integrating respiration rate (r_tot) from 0 to tinfl
    cr_tot_gr = first(integrate_smoother(t, r_tot, tinfl))    
    # estimate biomass increment from cumulative respiration
    # by computing with cumulative numbers -> can avoid rates mm_s, ks
    # cumulative resp due to both growth ((1-)u1) and maintenance u2
    #   cr = ((1-Y) u1 + u2) * t  
    # ratio of these two assumed fixed by ratio of uptake rates
    #    u1/u2 = ks r / (ksm + (1-r)) and assume ksm ~ 0.11 ks
    #    -> u2 = u1 * 0.11 (1-r)/r
    # factor out u1
    #    cr = u1 * [(1-y) + 0.11 (1-r)/r] * t
    # assume: bgr = u1*t
    #    cr = bgr * up     with up = (1-y) + 0.11 (1-r)/r
    #    bgr = cr/up
    u0t = replace_u0 ? u0 : copy(u0)
    up = (1 - Y0) + 0.11 * (1 - r0) / r0
    bgr = cr_tot_gr/up
    u0t.s = u0.s - cr_tot_gr - bgr
    u0t.b = b0 + bgr
    u0t, cr_tot_gr
end



function fit_growth_and_lim()
end