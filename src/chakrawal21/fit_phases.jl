function fit_initial_lim(tinfl, popt0; systemt = chak21_fixedr_system(), chak21syn = chak21syn, solver = Rodas5())
    iinfl = findfirst(chak21syn.solsyn.t .> tinfl) - 1

    ilim = (iinfl+1):length(chak21syn.solsyn.t)
    tlimo = chak21syn.solsyn.t[ilim]
    tlim = vcat(tinfl, tlimo)  # include inflection point
    obsr_tot = chak21syn.obsr_tot[ilim];
    obsq = chak21syn.obsq[ilim];
    #
    ps = LabeledParSetter(systemt)
    #u0g, cr_tot  = MicMods.estimate_u0_from_cumresp(tinfl, popt0, systemt, chak21syn; ps)
    u0g, cr_tot = estimate_u0_from_cumresp(tinfl, popt0, systemt, chak21syn; ps)
    s_plus_b = u0g.s + u0g.b
    @show sum(systemt.prob.u0)
    @show sum(u0g) + cr_tot
    probp = remake(systemt.prob; tspan = extrema(tlim), u0 = u0g)
    # @parameters t
    # @variables q(t) b(t)

    #p = getpopt(ps, systemt.prob, Val(true)) 
    # uses systemt, ps, probp, solver, obsq, obsr_tot, s_plus_b
    function floss_lim(p::AbstractVector{T}, ps=ps, par0 = probp.p, u00 = probp.u0) where T
        # pp = vcat(p[1:4], systemt.prob.p[5:end])
        # #u0p = convert.(eltype(p),systemt.prob.u0)
        # u0p = vcat(systemt.prob.u0[1], p[5], systemt.prob.u0[3])
        #pp, u0p = ps.setpu(p, systemt.prob)
        # change the optimized subset of parameter and initial state
        local pp, u0p = setpu(ps, p, par0, u00, Val(true))
        # adjust initial substrate for estimated biomass, so that sum 
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
        return lossval, solp, predq, predr_tot
    end
    floss_lim
end

function tmpf()
    poptl0 = copy(popt0);  poptl0.b = 41.0 # start with higher biomass
    floss_simp = fit_initial_lim(tinfl, popt0; systemt, chak21syn, solver);
    tmp = floss_simp(poptl0); first(tmp)
    floss1(p) = first(floss_simp(p))
    floss1(popt0)
    #@code_warntype floss_simp(poptl0)

    #using ForwardDiff
    ForwardDiff.gradient(floss1, poptl0)

    p = popt0; par0 = probp.p; u00 = probp.u0
    pp, u0p = setpu(ps, p, par0, u00, Val(true))
    #@code_warntype setpu(ps, p, par0, u00, Val(true))
    tmpf(ps, p, par0, u00) = setpu(ps, p, par0, u00, Val(true))
    #@code_warntype tmpf(ps, p, par0, u00)

#@code_warntype systemt.predout(solp.u[1], pp) # finally with let ok
    #@code_warntype setpu(ps, p, systemt.prob) # ok
    #@code_warntype floss_simp(p) # solp of type Any?glo
    #@code_warntype(ps.num[:q]) # ok
    lossval_true, solp_true, predq_true = floss_simp(popt0); lossval_true

    probp = remake(systemt.prob; tspan = extrema(tlim))
    probp.tspan
    systemt.prob.tspan # not changed

end

"""
    estimate_u0_from_cumresp(tinfl, popt0, systemt, chak21syn; ps = LabelledParSetter(systemt))

Obtain initial estimate of substrate s and biomass u by computing uptake from microbial parameters and respiration.

# Arguments
- `tinfl`: (numeric) start-time, after amendment
- `systemt`: tuple with ODE problem `prob`, and searchranges for subset of optimized parameters and initial states: `searchranges_p`, `searchranges_u0`
- `popt0`: vector of parameters and initial states to be optimized
- `chak21syn`: Tuple of respiration vector `obsr_tot` observed at time `solsyn.t`
- `b0`: microbial biomass at amendmend time 0. The default is obtained from chak21syn problem
- `r0`, `mm_s0`: representative phyisological state and  for the period 0..tinfl.

# Value
- Tuple `(u0, cr_tot_gr)` of parameters and cumulated respiration.
"""
function estimate_u0_from_cumresp(tinfl, popt0, systemt, chak21syn; 
    ps = LabelledParSetter(systemt), 
    b0 = label_statesys(chak21syn.pss, chak21syn.systemt.prob.u0).b,
    r0 = 0.6, # mm_s0 = 0.9,
    )
    # cumulative respiration: integrating respiration rate (r_tot) from 0 to tinfl
    cr_tot_gr = first(integrate_smoother(chak21syn.solsyn.t, chak21syn.obsr_tot, tinfl))    
    pp, u0p = setpu(ps, popt0, systemt.prob, Val(true))
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
    up = (1 - pp.Y) + 0.11 * (1 - r0) / r0
    bgr = cr_tot_gr/up
    u0p.s = u0p.s - cr_tot_gr - bgr
    u0p.b = b0 + bgr
    u0p, cr_tot_gr
end


function fit_initial_growth()
end

function fit_growth_and_lim()
end