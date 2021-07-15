using LsqFit: ForwardDiff
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
    serialize("inst/chak21_syn.jls", chak21syn);
end

using Serialization
chak21syn = deserialize("inst/chak21_syn.jls");

# find time of unlimited growth
tinfl = find_inflection(chak21syn.solsyn.t, chak21syn.obsr_tot)
iinfl = findfirst(chak21syn.solsyn.t .> tinfl)


# compute cumulated respiration:
# fit exponential model 
using LsqFit, GLM, StaticArrays
model(t, p) = p[1] * exp.(p[2] * t)
tdata = chak21syn.solsyn.t[1:iinfl]
ydata = chak21syn.obsr_tot[1:iinfl]
reg = lm(hcat(ones(length(tdata)),tdata), log.(ydata))
popt0 = [exp(coef(reg)[1]), coef(reg)[2]]
pf = coef(curve_fit(model, tdata, ydata, popt0))
ypred = model(tdata, pf)
function i_plot()
    scatter(tdata, ydata)
    plot!(tdata, ypred)
end
#pf[1]/pf[2] .* (exp.(pf[2] .* tdata) .- 1)
cr_tot_unlim = pf[1]/pf[2] * (exp(pf[2] * tinfl) - 1)

ilim = findall(chak21syn.solsyn.t .> tinfl)
tlimo = chak21syn.solsyn.t[ilim]
tlim = vcat(tinfl, tlimo)
obsr_tot = chak21syn.obsr_tot[ilim]
obsq = chak21syn.obsq[ilim]
chak21syn.solsyn(tinfl)

# fit the latter part using fixed pysiological state r
#systemt = chak21_simp_system();
systemt = chak21_fixedr_system();
#ps = ParSetter(systemt; stateopt=[]) # omit b0 from fitting
ps = ParSetter(systemt)
popt0 = getpopt(ps, systemt.prob, Val(true)) 
b00 = popt0.b
# here, adjust r parameter - for obs this is unknown
systemt.prob.p[parindex(ps, :r)] = 
    chak21syn.solsyn(tinfl)[stateindex(chak21syn.pss, :r)]
#popt_names = [ps.keys(systemt.searchranges_p)..., keys(systemt.searchranges_u0)...]
popt_names = [ps.paropt..., ps.stateopt...]


pp, u0p = setpu(ps, popt0, systemt.prob, Val(true))
# _ks, _r = map(sym -> systemt.prob.p[parindex(ps,sym)], (:ks, :r))
# _ksm = MicMods.estimate_ksm(_ks)
_r = 0.1
uptake = cr_tot_unlim/(1-pp.Y)
u0p.s = u0p.s + b00 - cr_tot_unlim
u0p.b = b00 + pp.Y * uptake * _r/(_r + 0.11*(1-_r))  # assume biomass yield of resp
systemt.prob.u0[:] .= u0p
popt0.b = u0p.b # update initial biomass in parameter vector
pt, u0t  = systemt.adjust_p0u0(pp, u0p)

solver = Rodas5()
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
    chak21syn.solsyn(tinfl)
    u0t
    chak21syn.solsyn(0)
    sum(chak21syn.solsyn(0)[1:2])
    chak21syn.solsyn(tinfl)[1:2]
    soltrue(tinfl)
    sum(soltrue(tinfl))
end

# obsqsyn = soltrue[ps.num[:r_tot]];
# σ_obsq = (t -> t[2] - t[1])(extrema(obsqsyn))/10
# obsq = obsqsyn .+ σ_obsq*randn(size(obsqsyn));
#pvec_tmp = copy(pvec_true); x0vec_tmp = copy(x0vec_true);
p = popt0
function getfloss_simp(systemt, solver, ps, obsq, obsr_tot)
    # put into closure so that types can be inferred
    function floss_(p)
        # pp = vcat(p[1:4], systemt.prob.p[5:end])
        # #u0p = convert.(eltype(p),systemt.prob.u0)
        # u0p = vcat(systemt.prob.u0[1], p[5], systemt.prob.u0[3])
        #pp, u0p = ps.setpu(p, systemt.prob)
        pp, u0p = setpu(ps, p, systemt.prob)
        pt, u0t = systemt.adjust_p0u0(pp, u0p)
        #@info typeof(p)
        probp = remake(systemt.prob; p = pt, u0 = u0t)
        #@show pt
        solp = solve(probp, 
            solver, 
            maxiters = 1000, # stop early if cannot determine sol
            #isoutofdomain = (u, p, t) -> (any(ui < zero(ui) for ui in u)),
        );
        if any(solp.retcode != :Success) 
            lossval = convert(eltype(p),Inf)::eltype(p) 
            predq = zeros(length(obsq))
            predr_tot = zeros(length(obsr_tot))
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
    @model function fitq(qobs, r_totobs,::Type{T} = Float64) where {T}
        p = Vector{T}(undef, length(sr))
        for (i,r) = enumerate(sr)
             p[i] ~ gettruncdist(r...)
        end
        if !isa(_context, Turing.PriorContext)
            lossval, solp, predq, predr_tot = floss_simp(p);
            if !isfinite(lossval) 
                Turing.@addlogprob! -Inf; return
            end
            #predq = Fill(sum(p), length(obsq))
            # for (i, qp) = enumerate(predq)
            #     qobs[i] ~ Normal(qp, σ_obsq)
            # end
            for (i, rp) = enumerate(predr_tot)
                r_totobs[i] ~ Normal(rp, σ_obsr_tot)
            end
        end
    end
    model_func = fitq(obsq, obsr_tot)
    #chn = chn0 = sample(model_func, MH(), 10)
    chn = chn0 = sample(model_func, NUTS(), 10)
    Random.seed!(0815)
    chn = chn0 = sample(model_func, NUTS(),  MCMCThreads(), 40, 3)
    chn = replacenames(chn0, Dict("p[$i]" => pname for (i,pname) in 
         enumerate(MicMods.replacetby0(popt_names))))
    # drop the 1st chain
    #chn = chn[:,:,[2,3]]
    # resample using 
    chn = chn1 = sample(model_func, NUTS(),  MCMCThreads(), 400, 3, 
        theta_init = MicMods.best_estimate(chn0));
#        theta_init = popt_true);
    chn = replacenames(chn1, Dict("p[$i]" => pname for (i,pname) in 
        enumerate(MicMods.replacetby0(popt_names))))
    tmp = Array(chn[:,:,[2,3]]);
    tmp = Array(chn[:,:,[3]]);
    tmp = Array(chn[:,:,[1]]);
    #chn = Chains(reshape(tmp, (size(tmp)...,1)), chn.name_map.parameters) 
    chn = Chains(reshape(tmp, (size(tmp)...,1)), popt_names) 
    function interactive_plotchains()
        popt1 = MicMods.label_popt(ps, MicMods.best_estimate(chn1))
        plot(chn)
        corner(chn)
        map(display, plotchain(get(chn, :log_density)));
        plot_post_and_priors(systemt, ps, chn)
        soltrue[ps.num[:s]]
    end
end

function depr_optim_initial_growth_closed()
    using LabelledArrays
    popt1 = MicMods.label_popt(ps, MicMods.best_estimate(chn1))
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
    using LabelledArrays
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
        plot!(predg.t, predgr_tot)
        label_parsys(chak21syn.pss, chak21syn.systemt.prob.p)
        label_statesys(chak21syn.pss, chak21syn.systemt.prob.u0)
        u0gopt
    end

    using Optim
    popt0 = label_popt(psg, getpopt(psg, systemtg.prob))
    popt = copy(popt0)
    function flossg(popt)
        pp, u0p = setpu(psg, popt, systemtg.prob)
        predg = solve(systemtg.prob, solver, u0 = u0p, p = pp, saveat=tlimgo);
        #pred_infl = predg[end]
        predgr_tot = predg[psg.num[:r_tot]];
        sum(abs2, predgr_tot - obsgr_tot)
    end
    flossg(popt0)
    flossg(popt0 .* 1.01)
    using ForwardDiff
    ForwardDiff.gradient(flossg, uopt)
    ansopt = optimize(flossg, uopt, LBFGS(); autodiff = :forward)
    Optim.minimizer(ansopt)

    popt = Optim.minimizer(ansopt)
    pp, u0p = setpu(psg, popt, systemtg.prob)
    predg = solve(systemtg.prob, solver, u0 = u0p, p = pp, saveat=tlimg);
    pred_infl = predg[end]
    predgr_tot = predg[psg.num[:r_tot]];
end




using StaticArrays
mod_growth = chak21_growth_closed();
x10 = mod_growth.predcum(10, values(mod_growth.x0), values(mod_growth.parms))
ts = 0:0.5:50
xg = map(ts) do t
    x1 = mod_growth.predcum(t, values(mod_growth.x0), values(mod_growth.parms))
end;
xga = VectorOfArray(xg);
plot(ts, xga[1,:])
plot(ts, xga[2,:])
plot(ts, xga[3,:], label = string(LabelledArrays.symnames(typeof(x10))[3]))
solsyn[:t]


function getfloss_comb(systemt, solver, ps, obsq)
    # put into closure so that types can be inferred
    function floss_(p)
        chak21_growth_closed
        # pp = vcat(p[1:4], systemt.prob.p[5:end])
        # #u0p = convert.(eltype(p),systemt.prob.u0)
        # u0p = vcat(systemt.prob.u0[1], p[5], systemt.prob.u0[3])
        #pp, u0p = ps.setpu(p, systemt.prob)

        pp, u0p = setpu(ps, p, systemt.prob)
        #@info typeof(p)
        probp = remake(systemt.prob; p = pp, u0 = u0p)
        #@show pp
        solp = solve(probp, 
            # solver, 
            maxiters = 1000, # stop early if cannot determine sol
            #isoutofdomain = (u, p, t) -> (any(ui < zero(ui) for ui in u)),
        );
        if any(solp.retcode != :Success) 
            lossval = convert(eltype(p),Inf)::eltype(p) 
            predq = zeros(length(obsq))
        else
            # fixed problems in gradient and is not type-stable
            predq = solp[ps.num[:q]]
            # predq = solp[ps.num[:q]]::typeof(obsq) 
            # predobs = systemt.predout.(solp.u, Ref(pp));
            # predq = [obs.q for obs in predobs];
            lossval = sum(abs2, predq .- obsq)
        end
        return lossval, solp, predq
    end
end
#@code_warntype systemt.predout(solp.u[1], pp) # finally with let ok
#@code_warntype setpu(ps, p, systemt.prob) # ok
floss_simp = getfloss_simp(systemt, solver, ps, obsq);
#@code_warntype floss_simp(p) # solp of type Any?glo
#@code_warntype(ps.num[:q]) # ok

p = getpopt(ps, systemt.prob); p .*= 0.8; lossval, solp, predq = floss_simp(p); lossval
popt_true = getpopt(ps, systemt.prob)
lossval_true, solp_true, predq_true = floss_simp(popt_true); lossval_true

