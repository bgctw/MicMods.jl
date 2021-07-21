function fit_initial_lim(tinfl, popt0; systemt = chak21_fixedr_system(), chak21syn = chak21syn, solver = Rodas5())
    iinfl = findfirst(chak21syn.solsyn.t .> tinfl) - 1

    ilim = (iinfl+1):length(chak21syn.solsyn.t)
    tlimo = chak21syn.solsyn.t[ilim]
    tlim = vcat(tinfl, tlimo)  # include inflection point
    obsr_tot = chak21syn.obsr_tot[ilim];
    obsq = chak21syn.obsq[ilim];
    #
    ps = LabeledParSetter(systemt)
    probp = remake(systemt.prob; tspan = extrema(tlim))
    #pg, u0g = MicMods.estimate_u0_from_cumresp(tinfl, popt0, systemt, chak21syn; ps)
    pg, u0g = estimate_u0_from_cumresp(tinfl, popt0, systemt, chak21syn; ps)
    probp.u0[:] .= u0g
    # @parameters t
    # @variables q(t) b(t)

    #p = getpopt(ps, systemt.prob, Val(true)) 
    # uses systemt, ps, probp, solver, obsq, obsr_tot)
    function floss_lim(p::AbstractVector{T}, ps=ps, par0 = probp.p, u00 = probp.u0) where T
        # pp = vcat(p[1:4], systemt.prob.p[5:end])
        # #u0p = convert.(eltype(p),systemt.prob.u0)
        # u0p = vcat(systemt.prob.u0[1], p[5], systemt.prob.u0[3])
        #pp, u0p = ps.setpu(p, systemt.prob)
        local pp, u0p = setpu(ps, p, par0, u00, Val(true))
        local pt, u0t = systemt.adjust_p0u0(pp, u0p)
        local solp = solve(probp, 
            solver, 
            maxiters = 1000, # stop early if cannot determine 
            # use unlabelled parameters to avoid many compilations
            p = pt.__x, u0 = u0t.__x,
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
            # local prednew = systemt.predout.(solp.u, Ref(pt.__x))
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
    floss_simp = fit_initial_lim(tinfl, popt0; systemt, chak21syn, solver);
    tmp = floss_simp(popt0); first(tmp)
    floss1(p) = first(floss_simp(p))
    floss1(popt0)
    #@code_warntype floss_simp(popt0)

    #using ForwardDiff
    ForwardDiff.gradient(floss1, popt0)

    p = popt0; par0 = probp.p; u00 = probp.u0
    pp, u0p = setpu(ps, p, par0, u00, Val(true))
    #@code_warntype setpu(ps, p, par0, u00, Val(true))
    tmpf(ps, p, par0, u00) = setpu(ps, p, par0, u00, Val(true))
    #@code_warntype tmpf(ps, p, par0, u00)

#@code_warntype systemt.predout(solp.u[1], pt) # finally with let ok
    #@code_warntype setpu(ps, p, systemt.prob) # ok
    #@code_warntype floss_simp(p) # solp of type Any?glo
    #@code_warntype(ps.num[:q]) # ok
    lossval_true, solp_true, predq_true = floss_simp(popt0); lossval_true

    probp = remake(systemt.prob; tspan = extrema(tlim))
    probp.tspan
    systemt.prob.tspan # not changed

end


function estimate_u0_from_cumresp(tinfl, popt0, systemt, chak21syn; ps = ParSetter(systemt))
    cr_tot_gr = first(integrate_smoother(chak21syn.solsyn.t, chak21syn.obsr_tot, tinfl))    
    pp, u0p = setpu(ps, popt0, systemt.prob, Val(true))
    # TODO: in real applications, estimate initial biomass by kinresp
    _b00 = label_statesys(chak21syn.pss, chak21syn.systemt.prob.u0).b
    _r = 0.8
    _mm_s = 0.9
    # estimate biomass increment from cumulative respiration
    # r = (Y u1 + u2) * b * mm_s
    #     assume mm_s = 0.9
    # u1/u2 = ks r / (ksm + (1-r)) and assume ksm ~ 0.11 ks
    #    u2 = 0.11 (1-r)/r
    bgr = cr_tot_gr / (pp.Y + 0.11 * (1-_r)/_r) / _mm_s
    u0p.s = u0p.s + _b00 - cr_tot_gr
    u0p.b = _b00 + bgr
    pp, u0p
end


function fit_initial_growth()
end

function fit_growth_and_lim()
end