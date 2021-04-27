#----------------- ModelingToolkit system
using ModelingToolkit, DifferentialEquations, DiffEqFlux, Flux, Optim

@parameters α β δ γ
@variables t x(t) y(t)
D = Differential(t)
eqs = [
  D(x) ~ α*x - β*x*y,
  D(y) ~ -δ*y + γ*x*y
]
@named lv = ODESystem(eqs)
syss = structural_simplify(lv) 
parms = [α => 1.5, β => 1.0, δ => 3.0, γ => 1.0]
x0 = [x => 1.0, y => 1.0]
tsteps = 0.0:0.1:10.0
prob = ODEProblem(syss, x0, extrema(tsteps), parms, jac = true)

# test solving with original and modified parameters
ptrue = [x.second for x in parms]
p0 = copy(ptrue)
p0[1] = 1.1
soltrue = solve(prob, Rodas5(), saveat = tsteps);
sol0 = solve(prob, Rodas5(), p=p0, saveat = tsteps);

callback = function (p, l, pred)
    display(l)
    # plt = plot(pred, ylim = (0, 6))
    # display(plt)
    return false # true to stop optimization
  end

obs = Array(soltrue) + 0.2*randn(size(sol))
function loss(p)
  solp = solve(prob, Rodas5(), p=p, saveat = tsteps)
  # sol = solve(prob, Rodas5(), p=p, saveat = tsteps,
  #    isoutofdomain=(u,p,t) -> u[1] < zero(u[1]))
  #lossval = sum(abs2, solp.-1)
  lossval = sum(abs2, solp .- obs)
  return lossval, solp
end

result_ode = DiffEqFlux.sciml_train(
  loss, p0, ADAM(0.1), 
  cb = callback, 
  maxiters = 20
  );

# result_ode = DiffEqFlux.sciml_train(
#   loss, result_ode.minimizer, BFGS(), #Descent(), #ADAM(0.1), 
#   cb = callback, 
#   maxiters = 30
#   );
# hcat(ptrue, result_ode.minimizer, p0)

result_ode = DiffEqFlux.sciml_train(
  loss, result_ode.minimizer, ADAM(0.1), #BFGS(), 
  cb = callback, 
  maxiters = 20#, allow_f_increases=true
  );

using Plots
plot(soltrue)
solp = solve(prob, Rodas5(), p=result_ode.minimizer, saveat = tsteps);
scatter(soltrue.t, obs[1,:])
scatter!(soltrue.t, obs[2,:])
plot!(solp)
  
