# https://github.com/SciML/ModelingToolkit.jl/issues/972#issue-859986238
using ModelingToolkit, DifferentialEquations

# lessons
# - test loss function with screwed parameters
# - check solver return code and return Inf misfit
# - in Turing model check for inifinite misfit
# - test Truing-model with single-core MC run

@parameters α β δ γ
@variables t x(t) y(t) dx(t) obs1(t)
D = Differential(t)
# based on example in https://diffeqflux.sciml.ai/stable/examples/optimization_ode/
eqs = [
  dx ~ α*x - β*x*y,  # testing observed variables
  D(x) ~ dx,
  D(y) ~ -δ*y + γ*x*y,
  obs1 ~ 2x + y,
]
@named lv = ODESystem(eqs)
syss = structural_simplify(lv) 
parms = [α => 1.5, β => 1.0, δ => 3.0, γ => 1.0]
x0 = [x => 1.0, y => 1.0]
tsteps = 0.0:0.1:10.0
prob = ODEProblem(syss, x0, extrema(tsteps), parms, jac = true)

# test solving with original and modified parameters
pvec = ModelingToolkit.varmap_to_vars(parms, parameters(syss))
x0vec = ModelingToolkit.varmap_to_vars(x0, states(syss))
p0 = copy(pvec); p0[1] = 1.1
#p0 = copy(pvec); p0[1] = 5 # stuck in local minimum
# using Rodas5 because its required in the real DAE problem 
soltrue = solve(prob, Rodas5(), saveat = tsteps);
obs_dt_true = soltrue[dx]; 
obs_dt = obs_dt_true .+ 0.8 .* rand(length(tsteps))

#sol0 = solve(prob, Rodas5(), p=p0, saveat = tsteps);
obsdt = soltrue[dx]
@code_warntype soltrue[dx]

using ForwardDiff, StaticArrays
ispopt = @SVector [false, true, false, false]
ipopt_for_ip = @SVector [0, 1, 0, 0]
floss_dt = function(popt)
  p = pvec .* zero(popt) .+ (ispopt[i] ? 
    popt[ipopt_for_ip[i]] : pvec[i] for i in axes(pvec,1)) 
  probp = remake(prob, p = p)
  solp = solve(probp, Rodas5(), saveat = tsteps);
  pred_dt = solp[dx]
  sum(abs2, pred_dt - obs_dt)
end
popt = pvec[[2]] .* 1.1
floss_dt(popt)
ForwardDiff.gradient(floss_dt, popt)

using DualNumbers
pd = map(Dual, prob.p)
u0d = map(Dual, prob.u0)
tmp1 = solve(prob, p = pd, saveat = tsteps);
# ambibuous method error
tmp2 = solve(prob, Rodas5(), p = pd, saveat = tsteps);
tmp3 = solve(prob, Rodas4(), p = pd, saveat = tsteps);

# works
tmp4 = solve(prob, AutoTsit5(Rosenbrock23()), p = pd, saveat = tsteps);
tmp5 = solve(prob, AutoTsit5(Rosenbrock23()), p = pd, u0 = u0d, saveat = tsteps);



obst = soltrue[obs1]
σ_obsq = 1.0
obstn = obst .+ σ_obsq .* randn(length(obst))

using Plots
# error: plot(soltrue, vars=[dx]) 
#ERROR: MethodError: no method matching *(::Vector{Float64}, ::Vector{Float64})
#plot(soltrue.t, soltrue[dx])
plot(soltrue.t, soltrue[obs1])
scatter!(soltrue.t, obstn)

callback = function (p, l, pred)
  display(l)
  # plt = plot(pred, ylim = (0, 6))
  # display(plt)
  return false # true to stop optimization
end

#obs = Array(soltrue)
function loss(p)
  solp = solve(prob, Rodas5(), p=p, saveat = tsteps, maxiters = 1e3)
  solp.retcode == :Success || return Inf, solp, 0.0
  # sol = solve(prob, Rodas5(), p=p, saveat = tsteps,
  #    isoutofdomain=(u,p,t) -> u[1] < zero(u[1]))
  # lossval = sum(abs2, solp .- obs) # this works
  # lossval = sum(abs2, solp[dx] .- obsdt)
  lossval = sum(abs2, solp[obs1] .- obstn)
  return lossval, solp, solp[obs1]
end

loss(pvec)[1]
loss(p0)[1]
loss([11.711229812609849, 1.739834865393254, -14.125420144916909, -6.177638231623366])

result_ode = DiffEqFlux.sciml_train(loss, p0, ADAM(0.1), cb = callback, maxiters = 500);
popt = result_ode.minimizer
hcat(pvec, p0, popt)


sol0 = solve(prob, Rodas5(), p=p0, saveat = tsteps);
solopt = solve(prob, Rodas5(), p=popt, saveat = tsteps);
# plot!( sol0.t, sol0[dx])
# plot!( solopt.t, solopt[dx])
plot!( sol0.t, sol0[obs1])
plot!( solopt.t, solopt[obs1])


function optim_turing_Gibbs()
  using Turing, FillArrays, MCMCChains
  using StatsPlots
  import Random
  #Turing.setadbackend(:zygote)
  Turing.setadbackend(:forwarddiff)
  # uses: floss, p0, obstn, σ_obsq
  @model function fitq(qobs, ::Type{T} = Float64) where {T}
      p = Vector{T}(undef, length(p0))
      for (i,r) = enumerate(p0)
           p[i] ~ Normal(0,10)
      end
      if !isa(_context, Turing.PriorContext)
          lossval, solp, pred_obs = loss(p);
          if !isfinite(lossval) 
              Turing.@addlogprob! -Inf; return
          end
          for (i, qp) = enumerate(pred_obs)
              qobs[i] ~ Normal(qp, σ_obsq)
          end
      end
  end
  model_func = fitq(obstn)
  #chn = chn0 = sample(model_func, MH(), 10)
  Random.seed!(0815)
  chn = chn0 = sample(model_func, NUTS(),  MCMCThreads(), 40, 3)
  # drop the 1st chain
  #chn = chn[:,:,[2,3]]
  # resample using 
  chn = chn1 = sample(model_func, NUTS(),  MCMCThreads(), 400, 3, 
      theta_init = MicMods.best_estimate(chn0));
#        theta_init = popt_true);
  # drop the 2st chain
  chn = chn[:,:,[1,3]]
  function interactive_plotchains()
      plot(chn)
      corner(chn)
      map(display, plotchain(get(chn, :log_density)));
      popt2 = MicMods.best_estimate(chn)
      solopt2 = solve(prob, Rodas5(), p=popt, saveat = tsteps);
      plot( solopt2.t, solopt2[obs1])
      plot!( solopt.t, solopt[obs1])
      plot!( sol0.t, sol0[obs1])
      scatter!(soltrue.t, obstn)
      hcat(pvec, p0, popt2)

    end
end

