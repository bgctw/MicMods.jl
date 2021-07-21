# based onf fitsimpmod.jl
#
# alternates between
# -fitting the growh phase 
#    for initial b r - given Y, ksm
# - and simplified model 
#    for ks, km, kd, Y - given initial b,s and assumring r = const = 1 

using MicMods
using ModelingToolkit, DifferentialEquations
using StatsPlots
using RecursiveArrayTools
using Statistics
using Random

gen_synthetic() = gen_synthetic(GLOBAL_RNG)
function gen_synthetic(rng)
    systemts = chak21_phys_system();
    systemts.prob.f
    pss = LabeledParSetter(systemts)
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

using StableRNGs
rng = StableRNG(123)
chak21syn = gen_synthetic(rng) #deserialize("inst/chak21_syn.jls");

# find time of unlimited growth
#tinfl = find_inflection(chak21syn.solsyn.t, chak21syn.obsr_tot)
tinfl0, ml = find_max(chak21syn.solsyn.t, chak21syn.obsr_tot);
tinfl = 0.9 * tinfl0 # start shortly before maximum

#------------- optimized initial decay
solver = Rodas5()
systemt = chak21_fixedr_system();
psl = LabeledParSetter(systemt)
popt0 = getpopt(psl, systemt.prob)
#poptg = copy(popt0); poptg.b = 40.0
ans_optlim = fit_initial_lim(tinfl, copy(popt0); systemt, chak21syn, solver);
floss_lim = ans_optlim.floss;
p = popt1 = ans_optlim.popt

function tmpf()
    pg = ans_optlim.popt0
    first(floss_lim(p))
    first(@inferred floss_lim(ans_optlim.popt0))
    @code_warntype floss_lim(p, ans_optlim.ps, ans_optlim.prob.p, ans_optlim.prob.u0)
    @code_warntype floss_lim(p)
    popt1 = first(ans_optlim)
    getpopt(chak21syn.pss, chak21syn.systemt.prob)
end

function iplot()
    # need to use loss-function instead of prob, because s must be adjusted
    soll = floss_lim(p).sol;
    label_statesys(psl, soll(tinfl))
    
    scatter(soll.t, ans_optlim.obsr_tot)
    plot!( soll.t, soll[psl.num[:r_tot]], label= "r_tot opt")
    plot!(chak21syn.solsyn.t, chak21syn.solsyn[psl.num[:r_tot]], label = "r_tot syn" )
    #scatter!(chak21syn.solsyn.t, chak21syn.obsr_tot, label = "r_tot")

    solg = flossg(poptg).sol;
    plot!( solg.t, solg[psl.num[:r_tot]], label= "r_tot opt gr")

    using StatsPlots
    plot(ans_optlim.chn)
    map(display, plotchain(get(chn, :log_density)));
    tmp = Array(chn[:,:,[1,2]]);

    ptrue = label_parsys(chak21syn.pss, chak21syn.systemt.prob.p)
    u0true = label_statesys(chak21syn.pss, chak21syn.solsyn(tinfl))
    poptt = label_popt(psl, vcat( ptrue[collect(getpopt(psl))], u0true[collect(getsopt(psl))] ))
    p = poptt
    u0true.s, u0true.b
    popt1.b, 
end

#------------- optimize initial growth
systemtg = chak21_phys_system();
psg = LabeledParSetter(systemtg, paropt=(), stateopt=(:b,:r))
pg[parindex.(Ref(psg), SA[:ks, :Y, :kd, :km])] .= popt1[[:ks, :Y, :kd, :km]];
systemtg.prob.p[[findfirst(==(sym), getparsys(psg)) for sym in getpopt(psl)]] .= popt1[1:length(getpopt(psl))]
popt0g = getpopt(psg, systemtg.prob) .* 1.1
#label_parsys(psg, systemtg.prob.p)
ans_g = fit_initial_growth(tinfl, copy(popt0g); systemtg, chak21syn, solver);
poptg = ans_g.popt
ans_g.ansopt
flossg = ans_g.floss

function tmp()
    # flossg = fit_initial_growth(tinfl, copy(popt0g); systemtg, chak21syn, solver)
    tmp = flossg(popt0g)
    @code_warntype flossg(popt0g)
    flossg1(p) = first(flossg(p))
    flossg1(popt0g.__x)
    using ForwardDiff
    ForwardDiff.gradient(flossg1, popt0g)
    #tmpf(p) = ForwardDiff.gradient(flossg1, popt0g)
    #@code_warntype tmpf(popt0g)
end

#----------------- optimize entire growth phase by Turing























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


