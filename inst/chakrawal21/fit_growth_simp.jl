# based onf fitsimpmod.jl
#
# alternates between
# -fitting the growh phase 
#    for initial b r - given Y, ksm
# - and simplified model 
#    for ks, km, kd, Y - given initial b,s and assumring r = const = 1 

using MicMods
using ModelingToolkit, DifferentialEquations
using Plots
using RecursiveArrayTools
using Statistics

function gen_synthetic()
    #using StableRNGs, Serialization
    rng = StableRNG(123)
    systemts = chak21_phys_system();
    pss = ParSetter(systemts)
    solsyn = solve(systemts.prob, Rodas5(), p=systemts.prob.p, u0 = systemts.prob.u0);
    t = solsyn.t;
    #plot(solsyn)
    obsqsyn = solsyn[pss.num[:q]];
    σ_obsq = (t -> t[2] - t[1])(extrema(obsqsyn))/10
    obsq = obsqsyn .+ σ_obsq*randn(rng,size(obsqsyn));
    obsr_totsyn = solsyn[pss.num[:r_tot]];
    σ_obsr_tot = (t -> t[2] - t[1])(extrema(obsr_totsyn))/10
    obsr_tot = obsr_totsyn .+ σ_obsr_tot*randn(rng,size(obsr_totsyn));
    function iplot()
        plot(solsyn)
        plot(t, solsyn[pss.num[:r]], label = "r")
        plot(t, obsqsyn, label = "q_true", xlab = "time (day)", ylab = "q")
        scatter!(t, obsq, label = "q")
        plot(t, obsr_totsyn, label = "r_tot_true", xlab = "time (day)", ylab = "resp (gC/m2/day)")
        scatter!(t, obsr_tot, label = "r_tot")
        plot(t, solsyn[pss.num[:r]])
    end
    chak21syn = (
        solsyn = solsyn, pss = pss, systemt = systemts,
        obsqsyn=obsqsyn, σ_obsq=σ_obsq, obsq=obsq,
        obsr_totsyn=obsr_totsyn, σ_obsr_tot=σ_obsr_tot, obsr_tot=obsr_tot);
    #serialize("inst/chak21_syn.jls", chak21syn); # does not work with closure functions
    chak21syn
end

chak21syn = gen_synthetic() #deserialize("inst/chak21_syn.jls");

# find time of unlimited growth
#tinfl = find_inflection(chak21syn.solsyn.t, chak21syn.obsr_tot)
tinfl0, ml = find_max(chak21syn.solsyn.t, chak21syn.obsr_tot);
tinfl = 0.9 * tinfl0 # start shortly before maximum

solver = Rodas5()
systemt = chak21_fixedr_system();
psl0 = LabeledParSetter(systemt)
popt0 = getpopt(psl0, systemt.prob)
tmp = fit_initial_lim(tinfl, popt0; systemt, chak21syn, solver)

function tmpf()
    #tmp = systemt.predout(systemt.prob.u0, systemt.prob.p)
    tmp = systemt.predout(first(solp.u), systemt.prob.p)
    @code_warntype systemt.predout(first(solp.u), systemt.prob.p)
    tmp.r_tot
    getindex(tmp, :r_tot)
    solp = solve(systemt.prob, solver)
    prednew = systemt.predout.(solp.u, Ref(systemt.prob.p))
    predr_tot = getindex.(prednew, :r_tot)
end


function tmpf()
    cr_tot_gr = first(integrate_smoother(chak21syn.solsyn.t, chak21syn.obsr_tot, tinfl))
    iinfl = findfirst(chak21syn.solsyn.t .> tinfl) - 1

    ilim = (iinfl+1):length(chak21syn.solsyn.t)
    tlimo = chak21syn.solsyn.t[ilim]
    tlim = vcat(tinfl, tlimo)
    obsr_tot = chak21syn.obsr_tot[ilim];
    obsq = chak21syn.obsq[ilim];
    label_statesys(chak21syn.pss, chak21syn.solsyn(tinfl))

    # fit the latter part using fixed pysiological state r
    #systemt = chak21_simp_system();
    systemt = chak21_fixedr_system();
    #ps = ParSetter(systemt; stateopt=[]) # omit b0 from fitting
    ps = ParSetter(systemt)
    popt0 = getpopt(ps, systemt.prob, Val(true)) 
    # here, adjust r parameter - for obs this is unknown
    systemt.prob.p[parindex(ps, :r)] = 
        chak21syn.solsyn(tinfl)[stateindex(chak21syn.pss, :r)]
    #popt_names = [ps.keys(systemt.searchranges_p)..., keys(systemt.searchranges_u0)...]
    popt_names = [ps.paropt..., ps.stateopt...]


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

    systemt.prob.u0[:] .= u0p
    popt0.b = u0p.b # update initial biomass in parameter vector (will be optimized)
    pt, u0t  = systemt.adjust_p0u0(pp, u0p)

    # adjust time in the problem to start after unlimited growth phase
    systemt = tmp = merge(systemt, (prob = remake(systemt.prob, 
        tspan = extrema(tlim), saveat=tlim),));
    systemt.prob.tspan
    soltrue = solve(systemt.prob, solver, p=pt, u0 = u0t);
    #plot(soltrue)
    soltrue[1] == u0t
    soltrue.t[1] == tinfl

    function i_plot()
        plot(tlim, soltrue[ps.num[:r_tot]])
        #scatter!(chak21syn.solsyn.t, chak21syn.obsr_totsyn, label = "r_tot")
        scatter!(tlimo, obsr_tot)
    end

    # obsqsyn = soltrue[ps.num[:r_tot]];
    # σ_obsq = (t -> t[2] - t[1])(extrema(obsqsyn))/10
    # obsq = obsqsyn .+ σ_obsq*randn(size(obsqsyn));
    #pvec_tmp = copy(pvec_true); x0vec_tmp = copy(x0vec_true);
    p = popt0
    function getfloss_simp(systemt, solver, ps, obsq, obsr_tot)
        # put into closure so that types can be inferred
        function floss_(p, par0 = systemt.prob.p, u00 = systemt.prob.u0)
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
    end
    #@code_warntype systemt.predout(solp.u[1], pt) # finally with let ok
    #@code_warntype setpu(ps, p, systemt.prob) # ok
    floss_simp = getfloss_simp(systemt, solver, ps, obsq, obsr_tot);
    #@code_warntype floss_simp(p) # solp of type Any?glo
    #@code_warntype(ps.num[:q]) # ok
    lossval_true, solp_true, predq_true = floss_simp(popt0); lossval_true
    lossval_true, solp_true, predq_true = floss_simp(popt0, systemt.prob.p, systemt.prob.u0); lossval_true
end

function plot_soltrue()
    plot(solp_true)
    plot!(tlim, solp_true[ps.num[:r_tot]])

end

function test_forwardDiff()
    # ForwardDiff
    ftest = function(p)
        lossval, solp, predq = floss_simp(p)
        lossval
    end
    tmp = ForwardDiff.gradient(ftest, p)
        
    p = popt_true
    #using DualNumbers
    pd = map(Dual, p)
    floss_simp(p)
    floss_simp(pd)
    pp, u0p = setpu(ps, p, systemt.prob)
    pp, u0p = setpu(ps, pd, systemt.prob)
    # #@show pp
    ftest = function(pp, u0p = systemt.prob.u0)
        probp = remake(systemt.prob; p = pp, u0 = u0p)
        # #@show pp
        solp = solve(probp, #solver, 
            maxiters = 1000, # stop early if cannot determine sol
            #isoutofdomain = (u, p, t) -> (any(ui < zero(ui) for ui in u)),
        );
        solp[end][1]
    end
    ftest(p)
    ftest(pd)
    pr = collect(p)
    tmp3 = ForwardDiff.gradient(ftest, pr)
end

function optim_turing_initialdecay()
    using Turing, FillArrays, MCMCChains
    using StatsPlots
    import Random
    using Suppressor
    #Turing.setadbackend(:zygote)
    Turing.setadbackend(:forwarddiff)
    # dn = Normal(0.5, 0.3)
    # pobs = rand(dn, 20)
    # @model function testpopt(pobs, ::Type{T} = Float64) where {T}
    #     # test sampling the mean of on p1
    #     # m ~ Normal(0, 2)
    #     #pobs ~ MvNormal(Fill(m, length(pobs)), dn.σ)
    #     p = Vector{T}(undef, 1)
    #     for i = 1:1
    #          p[i] ~ Normal(0, 2)
    #     end
    #     #p = [Normal(0,2) for i in [1]]
    #     pobs ~ MvNormal(Fill(p[1], length(pobs)), dn.σ)
    # end
    # model_func = testpopt(pobs)
    # chn = sample(model_func, NUTS(0.65), 1000, )
    #plot(chn)
    sr = vcat(
        [systemt.searchranges_p[key] for key in ps.paropt],
        [systemt.searchranges_u0[key] for key in ps.stateopt]
    )
    # uses: sr, floss, σ_obsq, σ_obsr_tot
    σ_obsq = chak21syn.σ_obsq
    σ_obsr_tot = chak21syn.σ_obsr_tot
    @model function fitql(qobsl, r_totobsl,::Type{T} = Float64) where {T}
        p = Vector{T}(undef, length(sr))
        for (i,r) = enumerate(sr)
             p[i] ~ gettruncdist(r...)
        end
        if !isa(_context, Turing.PriorContext)
            lossvall, solpl, predql, predr_totl = floss_simp(p);
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
    tmp = @suppress_err sample(model_func, NUTS(), 10)
    Random.seed!(0815)
    # chn = chn0 = sample(model_func, NUTS(),  MCMCThreads(), 40, 3)
    # chn = replacenames(chn0, Dict("p[$i]" => pname for (i,pname) in 
    #      enumerate(MicMods.replacetby0(getpoptnames(ps)))
    # drop the 1st chain
    #chn = chn[:,:,[2,3]]
    # resample using 
    chn = chn1 = sample(model_func, NUTS(),  MCMCThreads(), 400, 3, 
        #theta_init = MicMods.best_estimate(chn0)
        );
#        theta_init = popt_true);
    chn = replacenames(chn1, Dict("p[$i]" => pname for (i,pname) in 
          enumerate(MicMods.replacetby0(getpoptnames(ps)))))
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

function depr_optim_initial_growth_closed()
    using LabelledArrays
     ilimg = findall(chak21syn.solsyn.t .< tinfl)
    tlimgo = chak21syn.solsyn.t[ilimg]
    tlimg = vcat(tlimgo, tinfl)
    obsgr_tot = chak21syn.obsr_tot[ilimg]
    obsgq = chak21syn.obsq[ilimg]    
    #
    # growth system is not enough - need full system with substrate lim
    # for experimenter: enough substrate to estimate kmr
    systemtg = chak21_growth_closed(tend = tinfl)
    # labelled Array from Parm (OrderedDict Num -> val)
    pg0 = pg = LVector((; (k.val.name => v for (k,v) in systemtg.parms)...))
    pg.ks = popt1.ks
    pg.Y = popt1.Y
    u0g0 = collect(Float64, values(systemtg.x0))
    predg = systemtg.predcum.(tlimg, Ref(u0g0), Ref(pg));
    predga = VectorOfArray(predg);
    predgr_tot = predga[5,:]
    function i_plot()
        plot(chak21syn.solsyn.t, chak21syn.solsyn[chak21syn.pss.num[:r_tot]])
        plot!(chak21syn.solsyn.t[ilimg], chak21syn.solsyn[chak21syn.pss.num[:r_tot]][ilimg])
        scatter!(tlimgo, obsgr_tot)
        #plot(tlimgo, obsgq)
        #plot(tlimgo, chak21syn.solsyn[chak21syn.pss.num[:r]][ilimg])
        plot!(tlimg, predgr_tot)
        plot!(predg2.t, predg2r_tot)
        label_parsys(chak21syn.pss, chak21syn.systemt.prob.p)
        label_statesys(chak21syn.pss, chak21syn.systemt.prob.u0)
        u0gopt
    end
    systemtg2 = chak21_growth_system(;tend = tinfl);
    psg2 = ParSetter(systemtg2)
    pg20 = pg2 = label_parsys(psg2, systemtg2.prob.p)
    pg2[parindex.(Ref(psg2), SA[:ks, :Y, :kd])] .= popt1[[:ks, :Y, :kd]]
    u0g2 = label_statesys(psg2, systemtg2.prob.u0)
    predg2 = solve(systemtg2.prob, solver, u0 = u0g2, p = pg2);
    predg2r_tot = predg2[psg2.num[:r_tot]]

    probsyn = chak21syn.systemt.prob;
    solsyn = chak21syn.solsyn; pss = chak21syn.pss
    pg2[parindex.(Ref(psg2), SA[:ks, :Y, :kd])] .= 
        label_parsys(chak21syn.pss, probsyn.p )[[:ks, :Y, :kd]]
    # same initial state
    u0g2
    label_statesys(chak21syn.pss, probsyn.u0 )    
    # check parameters
    pg2
    label_parsys(chak21syn.pss, probsyn.p )    
    label_statesys(psg2, systemtg2.prob.f(u0g2, pg2, 0))
    label_statesys(chak21syn.pss, probsyn.f( probsyn.u0, probsyn.p, 0))
    predg2[psg2.num[:u1]][1], solsyn[pss.num[:u1]][1] # differ
    solsyn[pss.num[:mm_s]][1] # valid to neglect?

    using Optim
    u0g = copy(u0g0)
    uopt = u0g0[2:3]  # b,r
    function flossg(uopt)
        #u0g[2:3] .= uopt  # b, r
        local u0p = u0g .* zero(eltype(uopt)) .+ (i ∈ SA[2,3] ? 
            uopt[i-1] : u0g[i] for i in axes(u0g,1)) 
        predg = systemtg.predcum.(tlimg, Ref(u0p), Ref(pg));
        pred_infl = predg[end]
        predga = VectorOfArray(predg[1:(end -1)]);
        predgr_tot = predga[5,:];
        sum(abs2, predgr_tot - obsgr_tot)
    end
    flossg(uopt)
    flossg(uopt .* 1.01)
    using ForwardDiff
    ForwardDiff.gradient(flossg, uopt)
    ansopt = optimize(flossg, uopt, LBFGS(); autodiff = :forward)
    Optim.minimizer(ansopt)

    u0gopt = copy(u0g)
    u0gopt[2:3] .= Optim.minimizer(ansopt)
    predg = systemtg.predcum.(tlimg, Ref(u0gopt), Ref(pg));
    predga = VectorOfArray(predg);
    predgr_tot = predga[5,:]
    # infers much too small initial r
end

function optim_initial_growth()
    using LabelledArrays, StaticArrays
    popt1 = MicMods.label_popt(ps, MicMods.best_estimate(chn1))
    ilimg = findall(chak21syn.solsyn.t .< tinfl)
    tlimgo = chak21syn.solsyn.t[ilimg]
    tlimg = vcat(tlimgo, tinfl)
    obsgr_tot = chak21syn.obsr_tot[ilimg]
    obsgq = chak21syn.obsq[ilimg]    
    #
    # growth system is not enough - need full system with substrate lim
    # for experimenter: enough substrate to estimate kmr
    systemtg = chak21_phys_system();
    # adjust tspan in problem to be faster? or does saveat suffice?
    systemtg = merge(systemtg, (prob = remake(systemtg.prob, tspan = (0, tinfl)),))
    # labelled Array from Parm (OrderedDict Num -> val)
    # assume initial alues of cr_tot=0  and s given
    psg = ParSetter(systemtg, paropt = SA[], stateopt = SA[:b,:r])
    u0g = label_statesys(psg, systemtg.prob.u0)
    pg = label_parsys(psg, systemtg.prob.p)
    pg[parindex.(Ref(psg), SA[:ks, :Y, :kd, :km])] .= popt1[[:ks, :Y, :kd, :km]];
    pg
    predg = solve(systemtg.prob, solver, u0 = u0g, p = pg, saveat=tlimgo);
    predgr_tot = predg[psg.num[:r_tot]];
    function i_plot()
        plot(chak21syn.solsyn.t, chak21syn.solsyn[chak21syn.pss.num[:r_tot]])
        plot!(chak21syn.solsyn.t[ilimg], chak21syn.solsyn[chak21syn.pss.num[:r_tot]][ilimg])
        scatter!(tlimgo, obsgr_tot)
        #plot!(predg.t, predgr_tot)
        plot!(soloptg1.t, soloptg1[psg.num[:r_tot]])
        label_parsys(chak21syn.pss, chak21syn.systemt.prob.p)
        label_statesys(chak21syn.pss, chak21syn.systemt.prob.u0)
        u0gopt
    end

    using Optim
    popt0 = label_popt(psg, getpopt(psg, systemtg.prob))
    popt0 = LVector(b = _b00, r = 0.1)
    popt = copy(popt0)
    function getfloss_g(systemtg, solver, psg, obsgq, obsgr_tot, tlimgo = tlimgo)
        # put into closure so that types can be inferred
        function flossg(popt, par0 = systemtg.prob.p, u00 = systemtg.prob.u0)
            pp, u0p = setpu(psg, popt, par0, u00)
            #pp, u0p = setpu(psg, popt, systemtg.prob)
            predg = solve(systemtg.prob, solver, u0 = u0p, p = pp, saveat=tlimgo);
            if predg.retcode != :Success 
                lossval = convert(eltype(pp),Inf)::eltype(pp) 
                # todo check type in forwarddiff
                predgq = zeros(eltype(pp), length(obsgq))
                predgr_tot = zeros(eltype(pp), length(obsgr_tot))
            else
                #pred_infl = predg[end]
                predgq = predg[psg.num[:q]];
                predgr_tot = predg[psg.num[:r_tot]];
                lossval = sum(abs2, predgr_tot - obsgr_tot)
            end
            lossval, predg, predgq, predgr_tot
        end
    end
    flossg = getfloss_g(systemtg, solver, psg, obsgq, obsgr_tot)
    flossg1(popt) = flossg(popt)[1]
    flossg1(popt0)
    flossg1(popt0 .* 1.01)
    using ForwardDiff
    ForwardDiff.gradient(flossg1, popt)
    ansopt = optimize(flossg1, popt, LBFGS(); autodiff = :forward)
    Optim.minimizer(ansopt)
    poptg2 = label_popt(psg, Optim.minimizer(ansopt))

    pp, u0p = setpu(psg, poptg2, systemtg.prob)
    predg = solve(systemtg.prob, solver, u0 = u0p, p = pp, saveat=tlimg);
    pred_infl = predg[end]
    predgr_tot = predg[psg.num[:r_tot]];

    # explore uncertainties in growth curve    
    srg = vcat(
        [systemtg.searchranges_p[key] for key in psg.paropt],
        [systemtg.searchranges_u0[key] for key in psg.stateopt]
    )
    @model function fitqg(qobsg, r_totobsg,::Type{T} = Float64) where {T}
        p = Vector{T}(undef, length(srg))
        for (i,r) = enumerate(srg)
             p[i] ~ gettruncdist(r...)
        end
        if !isa(_context, Turing.PriorContext)
            # note flossg instead of floss
            lossvalg, solpg, predqg, predr_totg = flossg(p);
            if !isfinite(lossvalg) 
                Turing.@addlogprob! -Inf; return
            end
            #predq = Fill(sum(p), length(obsqg))
            # for (i, qp) = enumerate(predq)
            #     qobs[i] ~ Normal(qp, σ_obsq)
            # end
            for (i, rp) = enumerate(predr_totg)
                r_totobsg[i] ~ Normal(rp, σ_obsr_tot)
            end
        end
    end
    model_funcg = fitqg(obsgq, obsgr_tot)
    #tmp = sample(model_funcg, MH(), 2)
    chng = sample(model_funcg, NUTS(), 10)
    Random.seed!(0815)
    # resample using 
    chng = chng1 = sample(model_funcg, NUTS(),  MCMCThreads(), 400, 3, 
        theta_init = MicMods.best_estimate(chng));
#        theta_init = popt_true);
    chng = replacenames(chng1, Dict("p[$i]" => pname for (i,pname) in 
        enumerate(MicMods.replacetby0(getpoptnames(psg)))))
    tmp = Array(chng[:,:,[2,3]]);
    tmp = Array(chng[:,:,[3]]);
    tmp = Array(chng[:,:,[1]]);
    tmp = Array(chng[:,:,[1,3]]);
    #chng = Chains(reshape(tmp, (size(tmp)...,1)), chng.name_map.parameters) 
    chng = Chains(reshape(tmp, (size(tmp)...,1)), popt_names) 
    function interactive_plotchains()
        popt1 = MicMods.label_popt(ps, MicMods.best_estimate(chng1))
        plot(chng)
        corner(chng)
        map(display, plotchain(get(chng, :log_density)));
        plot_post_and_priors(systemt, ps, chng)
        soltrue[ps.num[:s]]
    end

    poptg1 = MicMods.label_popt(psg, MicMods.best_estimate(chng1))
    pg1, u0g1 = setpu(psg, poptg1, systemtg.prob)
    soloptg1 = solve(systemtg.prob, u0=u0g1, p=pg1, saveat = tlimgo);
    soloptg1(tinfl)
    predgr_tot = soloptg1[psg.num[:r_tot]];

end


function Turing_entire_Gibbs()
    #Turing.setadbackend(:zygote)
    Turing.setadbackend(:forwarddiff)

    # variable parameters for growth phase now include 
    # parameters for limited phase now include b,r
    psl = ParSetter(systemt, stateopt = SA[])
    srl2 = vcat(
        [systemt.searchranges_p[key] for key in psl.paropt],
        [systemt.searchranges_u0[key] for key in psl.stateopt]
    )
    # need loss function working with psl
    floss_simp2 = getfloss_simp(systemt, solver, psl, obsq, obsr_tot);


    popt1 = MicMods.label_popt(ps, MicMods.best_estimate(chn1))
    pl = label_popt(psl, [popt1[pn] for pn in getpoptnames(psl)])
    # 
    pg = poptg2

    # p = poptg1
    pfroml = SA[:Y, :km, :ks, :kd]                
    @model function fitqb(qobsg, r_totobsg, qobsl, r_totobsl, ::Type{T} = Float64) where {T}
        pg = label_popt(psg, Vector{T}(undef, length(srg)))
        #pg = Vector{T}(undef, length(srg))
        for (i,r) = enumerate(srg)
             pg[i] ~ gettruncdist(r...)
        end
        pl = label_popt(psl, Vector{T}(undef, length(srl2)))
        #pl =  Vector{T}(undef, length(srl2))
        for (i,r) = enumerate(srl2)
             pl[i] ~ gettruncdist(r...)
        end
        if !isa(_context, Turing.PriorContext)
            # other parameters in parg
            parg, u0g = setpu(psg, pg, systemtg.prob, Val(true))
            #parg, u0g = setpu(psg, pg, systemtg.prob)
            parg[parindex.(Ref(psg), pfroml)] .= pl[pfroml] # may only work with forwarddiff
            #  parga = parg .* zero(eltype(pl)) .+ (i ∈ SA[1,3,4,5] ?
            #     )
            lossvalg, solpg, predqg, predr_totg = flossg(pg, parg, u0g);
            #@show parg[4], lossvalg
            if !isfinite(lossvalg) 
                Turing.@addlogprob! -Inf; return
            end
            # set intial state b and r in parms and u0 for lim-floss
            ut = label_statesys(psg, solpg(tinfl))
            parl, u0l = setpu(psl, pl, systemt.prob, Val(true))
            parl.r = ut.r
            u0l.b = ut.b
            lossvall, solpl, predql, predr_totl = floss_simp2(pl, parl, u0l);
            if !isfinite(lossvall) 
                Turing.@addlogprob! -Inf; return
            end
            # # for (i, qp) = enumerate(predq)
            # #     qobs[i] ~ Normal(qp, σ_obsq)
            # # end
            # @show typeof(predr_totg), typeof(r_totobsg)

            for (i, rp) = enumerate(predr_totg)
                r_totobsg[i] ~ Normal(rp, σ_obsr_tot)
            end

            # # for (i, qp) = enumerate(predq)
            # #     qobsl[i] ~ Normal(qp, σ_obsq)
            # # end
            for (i, rp) = enumerate(predr_totl)
                r_totobsl[i] ~ Normal(rp, σ_obsr_tot)
            end
        end
    end
    model_funcb = fitqb(obsgq, obsgr_tot, obsq, obsr_tot)
    #tmp = sample(model_funcb, MH(), 2)
    tmp = sample(model_funcb, NUTS(), 2, theta_init = vcat(pg,pl))
    Random.seed!(0815)
    # resample using 
    chng = chng1 = sample(model_funcb, NUTS(),  MCMCThreads(), 400, 3, 
#        theta_init = MicMods.best_estimate(chng)
#        theta_init = popt_true
        theta_init = collect(vcat(pg,pl))
    );
    chng = replacenames(chng1, merge(
        Dict("pg[$i]" => pname for (i,pname) in 
        enumerate(getpoptnames(psg))),
        Dict("pl[$i]" => pname for (i,pname) in 
        enumerate(getpoptnames(psl)))
        ))
    tmp = Array(chng[:,:,[2,3]]);
    tmp = Array(chng[:,:,[3]]);
    tmp = Array(chng[:,:,[1]]);
    tmp = Array(chng[:,:,[1,3]]);
    #chng = Chains(reshape(tmp, (size(tmp)...,1)), chng.name_map.parameters) 
    chng = Chains(reshape(tmp, (size(tmp)...,1)), popt_names) 
    function interactive_plotchains()
        popt1 = MicMods.label_popt(ps, MicMods.best_estimate(chng1))
        plot(chng)
        corner(chng)
        map(display, plotchain(get(chng, :log_density)));
        plot_post_and_priors(systemt, ps, chng)
        soltrue[ps.num[:s]]
    end

    poptg1 = MicMods.label_popt(psg, MicMods.best_estimate(chng1))
    pg1, u0g1 = setpu(psg, poptg1, systemtg.prob)
    soloptg1 = solve(systemtg.prob, u0=u0g1, p=pg1, saveat = tlimgo);
    soloptg1(tinfl)
    predgr_tot = soloptg1[psg.num[:r_tot]];

end


