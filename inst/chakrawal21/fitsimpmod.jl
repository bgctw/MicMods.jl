using MicMods
using ModelingToolkit, DifferentialEquations
using Plots
using DiffEqFlux, Flux
using RecursiveArrayTools

# solve physiological and simple system with straw-N parameters
# for kmr small r->1 they should yield the same prediction 

# system = chak21_simp_system();
# system = chak21_growth_system();
# system = chak21_simp_system();

fsystems = Dict(
    :growth => chak21_growth_system, 
    :simp => chak21_simp_system, 
    #:phys => chak21_phys_system,
    )
systems = Dict(map(zip(keys(fsystems),values(fsystems))) do (key, f)
    Pair(key, f())
end); 
system = first(systems).second;

@variables t q(t) r_tot(t)
#@variables s(t) b(t) cr_tot(t) dec_s(t) r(t)
#@variables dec_s(t) tvr_b(t)

mapdict(f,d) = Dict(pair.first => f(pair.second) for pair in d)
# d2 = Dict(:a => (1,2), :b => (3,4))
# mapdict(d2) do x
#     2 .* x
# end
    

function compute_sols()
    function solve_chak21(prob, x0vec, pvec)
        #sol = solve(prob, Rodas5(), p=p, x0 = x0, isoutofdomain=hasnegativestate);
        sol = solve(prob, AutoTsit5(Rosenbrock23()), p = pvec, x0 = x0vec);
        @assert all(sol[system.states[:s]] .+ sol[system.states[:b]] .+ 
            sol[system.states[:cr_tot]] .≈ sum(x0vec[1:2]))
        sol
    end
    sols = mapdict(systems) do system
        pvec = ModelingToolkit.varmap_to_vars(system.parms, parameters(system.syss))
        x0vec = ModelingToolkit.varmap_to_vars(system.x0, states(system.syss))
        sol = solve_chak21(system.prob, x0vec, pvec)
    end;
    @assert all(isapprox.(sols[:phys](48.0)[1:3], sols[:simp](48.0), atol=0.01))
    @assert all(isapprox.(sols[:simp](48.0), sols[:phys](48.0)[1:3], atol=0.01))
end

sysid = :simp
sysid = :growth
#sol = sols[sysid]; 
system = systems[sysid];


function tmp_inspectsol()
    #plot(sol)
    #plot(sol,vars=[s, b, cr_tot, dec_s, r_tot])
    plot(sol,vars=getindex.(Ref(system.states),[:s,:b,:cr_tot]))
    plot!(sol.t, sol[system.observed[:r_tot]], label="r_tot")
    plot(sol.t, sol[system.observed[:q]], label="q")
    scatter!(sol.t, obsq, label="obsq")
    #scatter!(sol.t, predq, label="predq_i")
    plot!(solp.t, solp[system.observed[:q]], label = "q_pred")
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
DifferentialEquations.parameters(system.syss)
states(system.syss)
setpu, getpopt, poptnums = constructparsetter(system);
p = getpopt(system.prob) 
pp, u0p = setpu(p, system.prob)
@assert all(pp .== system.prob.p)
@assert all(u0p .== system.prob.u0)
# todo write macro to generate loss function for specific parameters
#pvec_true = ModelingToolkit.varmap_to_vars(system.parms, parameters(system.syss))
#x0vec_true = ModelingToolkit.varmap_to_vars(system.x0, states(system.syss))
solver = AutoTsit5(Rosenbrock23());
#solver = Rodas4();
soltrue = solve(system.prob, solver, p=system.prob.p, x0 = system.prob.u0);
obsqtrue = soltrue[system.observed[:q]];
obsq = obsqtrue .+ median(obsqtrue)/10*randn(size(obsqtrue));
#pvec_tmp = copy(pvec_true); x0vec_tmp = copy(x0vec_true);
function floss(p)
    # pp = vcat(p[1:4], system.prob.p[5:end])
    # #u0p = convert.(eltype(p),system.prob.u0)
    # u0p = vcat(system.prob.u0[1], p[5], system.prob.u0[3])
    pp, u0p = pset(p, system.prob.p, system.prob.u0)
    probp = remake(system.prob; p = pp, u0 = u0p)
    #@show pp
    #solp = solve(probp; saveat = tsteps)
    solp = solve(probp; solver, 
        maxiters = 1000, # stop early if cannot determine sol
        isoutofdomain = (u, p, t) -> (any(ui < zero(ui) for ui in u)),
    );
    lossval = any(solp.retcode != :Success) ? convert(eltype(p),Inf)::eltype(p) : begin
      predobs = system.predout.(solp.u, Ref(pp), Ref(system));
      predq = [obs.q for obs in predobs];
      #predq = solp[q]
      lossval = sum(abs2, predq .- obsq)
    end
    return lossval, solp
end
#@code_warntype(floss(p)) # solp of type Any?
p[1] = 0.1; p[2:5] .*= 0.8; lossval, solp = floss(p); lossval
popt_true = vcat(system.prob.p[1:4], system.prob.u0[2])
lossval_true, solp_true = floss(popt_true); lossval_true

function optim_Zygote()
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
    lossval, solp = floss(result_ode.minimizer); lossval
    hcat(system.prob.p[1:3], result_ode1.minimizer, p )
end

function optim_Optim()
    using Optim, NLSolversBase
    floss1 = OnceDifferentiable(x -> floss(x)[1], p; autodiff=:forward);
    result_ode2 = optimize(floss1, p, Optim.Options(show_trace = true, time_limit = 2*60.0)) # faster than with autodiff here
    #result_ode2 = optimize(floss1, p, LBFGS()) # faster than with autodiff here
    #result_ode2 = optimize(x -> floss(x)[1], p, LBFGS(), Optim.Options(show_trace = true, time_limit = 2*60.0); autodiff = :forward) 
    hcat(system.prob.p[1:3], result_ode2.minimizer, p )
    lossval2, solp2 = floss(result_ode2.minimizer); lossval2
end

function optim_blackbox()
    using BlackBoxOptim
    # function rosenbrock2d(x)
    #     return (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
    # end
    # res = bboptimize(rosenbrock2d; SearchRange = (-5.0, 5.0), NumDimensions = 2);
    # best_candidate(res)
    DifferentialEquations.parameters(system.syss)      
    floss1b(p) = floss(p)[1]
    res = bboptimize(floss1b; SearchRange = [
        (0.005, 1.0),   #ks
        (0.005, 1.0),   #ksm
        (0.1, 10.0),    #km
        (1/200, 1/10),  #kd
        (0.1, 0.8),  #Y
        #(0.1, 1/10),  #s0 is fixed
        (1.0, 100.0),  #b0
        ]);
    hcat(popt_true, best_candidate(res), p )
    lossval, solp = floss(best_candidate(res)); lossval
  
end








