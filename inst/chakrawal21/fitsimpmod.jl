using MicMods
using ModelingToolkit, DifferentialEquations
using Plots
using RecursiveArrayTools

# solve physiological and simple system with straw-N parameters
# for kmr small r->1 they should yield the same prediction 

# systemt = chak21_simp_system();
# systemt = chak21_growth_system();
# systemt = chak21_simp_system();

fsystems = Dict(
    #:growth => chak21_growth_system, 
    :simp => chak21_simp_system, 
    #:phys => chak21_phys_system,
    )
systems = Dict(map(zip(keys(fsystems),values(fsystems))) do (key, f)
    Pair(key, f())
end); 
sysid = :simp
sysid = :growth
sysid = first(systems).first
systemt = systemt = systems[sysid];


#@variables t q(t) r_tot(t)
#@variables s(t) b(t) cr_tot(t) dec_s(t) r(t)
#@variables dec_s(t) tvr_b(t)

#mapdict(f,d) = Dict(pair.first => f(pair.second) for pair in d)
# d2 = Dict(:a => (1,2), :b => (3,4))
# mapdict(d2) do x
#     2 .* x
# end
    

function compute_sols()
    function solve_chak21(prob, x0vec, pvec)
        #sol = solve(prob, Rodas5(), p=p, x0 = x0, isoutofdomain=hasnegativestate);
        sol = solve(prob, AutoTsit5(Rosenbrock23()), p = pvec, x0 = x0vec);
        @assert all(sol[systemt.states[:s]] .+ sol[systemt.states[:b]] .+ 
            sol[systemt.states[:cr_tot]] .≈ sum(x0vec[1:2]))
        sol
    end
    sols = mapdict(systems) do systemt
        pvec = ModelingToolkit.varmap_to_vars(systemt.parms, parameters(systemt.syss))
        x0vec = ModelingToolkit.varmap_to_vars(systemt.x0, states(systemt.syss))
        sol = solve_chak21(systemt.prob, x0vec, pvec)
    end;
    @assert all(isapprox.(sols[:phys](48.0)[1:3], sols[:simp](48.0), atol=0.01))
    @assert all(isapprox.(sols[:simp](48.0), sols[:phys](48.0)[1:3], atol=0.01))
end


function tmp_inspectsol()
    #plot(sol)
    #plot(sol,vars=[s, b, cr_tot, dec_s, r_tot])
    sol = soltrue;
    #sol = solp_true;
    plot(sol)
    plot(sol,vars=[ps.num[:r]])
    plot(sol,vars=getindex.(Ref(systemt.states),[:s,:b,:cr_tot]))
    plot!(sol.t, sol[systemt.observed[:r_tot]], label="r_tot")
    plot(sol.t, sol[systemt.observed[:q]], label="q")
    #plot(sol.t, predq_true, label = "predq_true")
    #plot(sol.t, predq, label = "predq")
    scatter!(sol.t, obsq, label="obsq")
    #scatter!(sol.t, predq, label="predq_i")
    plot!(solp.t, solp[ps.num[:q]], label = "q_pred")
    xlim = (17,18)
    bt = (sol.t .>= 17) .& (sol.t .<= 18)
    plot(sol,vars=[s, b, cr_tot, dec_s, r_tot], xlim = xlim)
    plot(sol,vars=[u1,u2,s], xlim = xlim)
    plot(sol[s,16:25],vars=[s])

    # glucose taken up within 16 hours, but stored in microbial biomass
    # instead of respired
    plot(sol,vars=[q], ylim = (0, 800e-6)) 
    plot(sol,vars=[r_tot], ylim = (0,2e-6))
    # in the model resp and q only differ by factor (1-Y)/HG 
end


using Statistics
# ps = constructparsetter(systemt);
# p = ps.getpopt(systemt.prob) 
#ps = ParSetter(systemt);
pss = Dict(
    :phys => ParSetter(systemt; stateopt=[]),
    :growth => ParSetter(systemt; paropt=[])
);
ps = pss[sysid]
p = getpopt(ps, systemt.prob) 
#popt_names = [ps.keys(systemt.searchranges_p)..., keys(systemt.searchranges_u0)...]
popt_names = [ps.paropt..., ps.stateopt...]


#pvec_true = ModelingToolkit.varmap_to_vars(systemt.parms, parameters(systemt.syss))
#x0vec_true = ModelingToolkit.varmap_to_vars(systemt.x0, states(systemt.syss))
solver = AutoTsit5(Rosenbrock23());
#solver = Rodas4();
soltrue = solve(systemt.prob, solver, p=systemt.prob.p, x0 = systemt.prob.u0);
obsqtrue = soltrue[ps.num[:q]];
σ_obsq = (t -> t[2] - t[1])(extrema(obsqtrue))/10
obsq = obsqtrue .+ σ_obsq*randn(size(obsqtrue));
#pvec_tmp = copy(pvec_true); x0vec_tmp = copy(x0vec_true);
function getfloss(systemt, solver, ps, obsq)
    # put into closure so that types can be inferred
    function floss_(p)
        # pp = vcat(p[1:4], systemt.prob.p[5:end])
        # #u0p = convert.(eltype(p),systemt.prob.u0)
        # u0p = vcat(systemt.prob.u0[1], p[5], systemt.prob.u0[3])
        #pp, u0p = ps.setpu(p, systemt.prob)
        pp, u0p = setpu(ps, p, systemt.prob)
        probp = remake(systemt.prob; p = pp, u0 = u0p)
        # #@show pp
        solp = solve(probp, solver, 
            maxiters = 1000, # stop early if cannot determine sol
            #isoutofdomain = (u, p, t) -> (any(ui < zero(ui) for ui in u)),
        );
        if any(solp.retcode != :Success) 
            lossval = convert(eltype(p),Inf)::eltype(p) 
            predq = zeros(length(obsq))
        else
            # fixed problems in gradient and is not type-stable
           predq = solp[ps.num[:q]]::typeof(obsq) 
        #   predobs = systemt.predout.(solp.u, Ref(pp));
        #   predq = [obs.q for obs in predobs];
          lossval = sum(abs2, predq .- obsq)
        end
        return lossval, solp, predq
    end
end
#@code_warntype systemt.predout(solp.u[1], pp) # finally with let ok
#@code_warntype setpu(ps, p, systemt.prob) # ok
floss = getfloss(systemt, solver, ps, obsq);
#@code_warntype floss(p) # solp of type Any?glo
#@code_warntype(ps.num[:q]) # ok
@code_warntype solp_true[]

p = getpopt(ps, systemt.prob); p .*= 0.8; lossval, solp, predq = floss(p); lossval
popt_true = getpopt(ps, systemt.prob)
lossval_true, solp_true, predq_true = floss(popt_true); lossval_true

function optim_Zygote()
    using DiffEqFlux, Flux
    using Zygote
    gr = gradient(x -> floss(x)[1],p)
    gr = gradient(x -> floss(x)[1],pvec_true[[1]])

    callback = function (p, lossval, solp)
        display(lossval)
        # plt = plot(solp.t, solp[q])
        # scatter!(plt, solp.t, qobs)
        # display(plt)
        # Tell sciml_train to not halt the optimization. If return true, then
        # optimization stops.
        return false
    end

    result_ode1 = result_ode = DiffEqFlux.sciml_train(
        floss, p, ADAM(0.1), cb = callback, maxiters = 100)
    lossval, solp, predq = floss(result_ode.minimizer); lossval
    hcat(systemt.prob.p[1:3], result_ode1.minimizer, p )
end

function optim_Optim()
    using Optim, NLSolversBase
    floss1 = OnceDifferentiable(x -> floss(x)[1], p; autodiff=:forward);
    result_ode2 = optimize(floss1, p, Optim.Options(show_trace = true, time_limit = 2*60.0)) # faster than with autodiff here
    #result_ode2 = optimize(floss1, p, LBFGS()) # faster than with autodiff here
    #result_ode2 = optimize(x -> floss(x)[1], p, LBFGS(), Optim.Options(show_trace = true, time_limit = 2*60.0); autodiff = :forward) 
    hcat(systemt.prob.p[1:3], result_ode2.minimizer, p )
    lossval2, solp2, predq = floss(result_ode2.minimizer); lossval2
end

function bboptimize_ctrl(functionOrProblem, parameters::BlackBoxOptim.Parameters = BlackBoxOptim.EMPTY_PARAMS; kwargs...)
    optctrl = BlackBoxOptim.bbsetup(functionOrProblem, parameters; kwargs...)
    res = BlackBoxOptim.run!(optctrl)
    res, optctrl
end

function optim_blackbox()
    using BlackBoxOptim
    # function rosenbrock2d(x)
    #     return (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
    # end
    # res = bboptimize(rosenbrock2d; SearchRange = (-5.0, 5.0), NumDimensions = 2);
    # best_candidate(res)
    DifferentialEquations.parameters(systemt.syss)      
    floss1b(p) = floss(p)[1]
    sr = [systemt.searchranges_p..., systemt.searchranges_u0...]
    res, optctrl = bboptimize_ctrl(floss1b; SearchRange = [p.second for p in sr], MaxTime = 20);
    hcat(popt_true, best_candidate(res), p )
    lossval, solp = floss(best_candidate(res)); lossval

    arch = optctrl.runcontrollers[1].evaluator.archive;
    cands = arch.candidates;
    [fitness(x) for x in cands]
    pa = VectorOfArray([x.params for x in cands]);
    #size(pa)
    pa[1,:] # 1: km almost all the same - no spread
    pa[3,:] # 3: Y
    # to explore the model, use Bayesian, look into Turing.jl and Tpapp
end

function test_distfit()
    r = systemt.searchranges_p[ps.num[:Y]]
    dn = truncnormalfromrange(r[1:2]...)
    hcat(quantile.(dn, [0.025, 0.975]), collect(r))
    dln = trunclognormalfromrange(r[1:2]...)
    hcat(quantile.(dln, [0.025, 0.975]), collect(r))
    dltn = trunclogitnormalfromrange(r[1:2]...)
    hcat(quantile.(dltn, [0.025, 0.975]), collect(r))
    r = systemt.searchranges_p[ps.num[:ks]]
    r = systemt.searchranges_u0[ps.num[:r]]
    dp = gettruncdist(r...)
    function tmp()
        using StatsPlots
        plot(dn)
        plot!(dln)
        plot!(dltn)
        plot(dp)
    end
end


function optim_turing()
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
    # uses: sr, floss, σ_obsq
    @model function fitq(qobs, ::Type{T} = Float64) where {T}
        p = Vector{T}(undef, length(sr))
        for (i,r) = enumerate(sr)
             p[i] ~ gettruncdist(r...)
        end
        if !isa(_context, Turing.PriorContext)
            lossval, solp, predq = floss(p);
            if !isfinite(lossval) 
                Turing.@addlogprob! -Inf; return
            end
            #predq = Fill(sum(p), length(obsq))
            for (i, qp) = enumerate(predq)
                qobs[i] ~ Normal(qp, σ_obsq)
            end
        end
    end
    model_func = fitq(obsq)
    #chn = chn0 = sample(model_func, NUTS(), 10)
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
        enumerate(replacetby0(popt_names))))
    tmp = Array(chn[:,:,[2,3]]);
    chn = Chains(reshape(tmp, (size(tmp)...,1)), chn.name_map.parameters) 
    function interactive_plotchains()
        plot(chn)
        corner(chn)
        map(display, plotchain(get(chn, :log_density)));
        plot_post_and_priors(systemt, ps, chn)
    end
end

function tmp_predict()
    using AxisArrays
    pa0 = Array(chn)
    pa = AxisArray(Array(chn); row = axes(pa0,1), col = chn.name_map.parameters)

    vec(pa[:,:kd])
    lossval, solp, predq = floss(pa0[1,:])
    plot(solp.t, predq, label = "model fit")
    scatter!(solp.t, obsq, label = "observations", ylab="heat", xlab = "time (hours)")
end


function optim_turing_Gibbs()
    using Turing, FillArrays, MCMCChains
    using StatsPlots
    import Random
    #Turing.setadbackend(:zygote)
    Turing.setadbackend(:forwarddiff)
    sr = vcat(
        [systemt.searchranges_p[key] for key in ps.paropt],
        [systemt.searchranges_u0[key] for key in ps.stateopt]
    )
    # uses: sr, floss, σ_obsq
    @model function fitq(qobs, ::Type{T} = Float64) where {T}
        p = Vector{T}(undef, length(sr))
        for (i,r) = enumerate(sr)
             p[i] ~ gettruncdist(r...)
        end
        if !isa(_context, Turing.PriorContext)
            lossval, solp, predq = floss(p);
            if !isfinite(lossval) 
                Turing.@addlogprob! -Inf; return
            end
            #predq = Fill(sum(p), length(obsq))
            for (i, qp) = enumerate(predq)
                qobs[i] ~ Normal(qp, σ_obsq)
            end
        end
    end
    model_func = fitq(obsq)
    #chn = chn0 = sample(model_func, NUTS(), 10)
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
        enumerate(replacetby0(popt_names))))
    tmp = Array(chn[:,:,[2,3]]);
    chn = Chains(reshape(tmp, (size(tmp)...,1)), chn.name_map.parameters) 
    function interactive_plotchains()
        plot(chn)
        corner(chn)
        map(display, plotchain(get(chn, :log_density)));
        plot_post_and_priors(systemt, ps, chn)
    end
end








