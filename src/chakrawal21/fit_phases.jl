function fit_initial_lim(tinfl; systemt = chak21_fixedr_system(), chak21syn = chak21syn, solver = Rodas5())
    iinfl = findfirst(chak21syn.solsyn.t .> tinfl) - 1

    ilim = (iinfl+1):length(chak21syn.solsyn.t)
    tlimo = chak21syn.solsyn.t[ilim]
    tlim = vcat(tinfl, tlimo)  # include inflection point
    obsr_tot = chak21syn.obsr_tot[ilim];
    obsq = chak21syn.obsq[ilim];
    #
    pp, u0p = estimate_u0_from_cumresp(tinfl, systemt, chak21syn)
    ps = ParSetter(systemt)
    systemt.prob.u0[:] .= u0p

    p = popt0
    function floss_lim(p, ps, par0 = systemt.prob.p, u00 = systemt.prob.u0)
        # pp = vcat(p[1:4], systemt.prob.p[5:end])
        # #u0p = convert.(eltype(p),systemt.prob.u0)
        # u0p = vcat(systemt.prob.u0[1], p[5], systemt.prob.u0[3])
        #pp, u0p = ps.setpu(p, systemt.prob)
        pp, u0p = setpu(ps, p, par0, u00)
        pt, u0t = systemt.adjust_p0u0(pp, u0p)
        #@info typeof(p)
        probp = remake(systemt.prob; p = pt, u0 = u0t)
        #@show pt
        solp = solve(probp, 
            solver, 
            maxiters = 1000, # stop early if cannot determine sol
            #isoutofdomain = (u, p, t) -> (any(ui < zero(ui) for ui in u)),
        );
        if solp.retcode != :Success 
            lossval = convert(eltype(p),Inf)::eltype(p) 
            predq = zeros(eltype(pp), length(obsq))
            predr_tot = zeros(eltype(pp), length(obsr_tot))
        else
            # fixed problems in gradient and is not type-stable
            predq = solp[ps.num[:q]][2:end]
            #lossval = sum(abs2, predq .- obsq)
            predr_tot = solp[ps.num[:r_tot]][2:end]
            lossval = sum(abs2, predr_tot .- obsr_tot)
        end
        return lossval, solp, predq, predr_tot
    end

    floss_lim
end

function tmpf()
    #@code_warntype systemt.predout(solp.u[1], pt) # finally with let ok
    #@code_warntype setpu(ps, p, systemt.prob) # ok
    floss_simp = getfloss_simp(systemt, solver, ps, obsq, obsr_tot);
    #@code_warntype floss_simp(p) # solp of type Any?glo
    #@code_warntype(ps.num[:q]) # ok
    lossval_true, solp_true, predq_true = floss_simp(popt0); lossval_true
    lossval_true, solp_true, predq_true = floss_simp(popt0, systemt.prob.p, systemt.prob.u0); lossval_true

end


function estimate_u0_from_cumresp(tinfl, systemt, chak21syn)
    cr_tot_gr = first(integrate_smoother(chak21syn.solsyn.t, chak21syn.obsr_tot, tinfl))    
    ps = ParSetter(systemt)
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