using MicMods
using ModelingToolkit, DifferentialEquations
using Plots

# solve physiological and simple system with straw-N parameters
# for kmr small r->1 they should yield the same prediction 

fsystems = Dict(
    :simp => chak21_simp_system, 
    :phys => chak21_phys_system,
    )

systems = Dict(map(zip(keys(fsystems),values(fsystems))) do (key, f)
    Pair(key, f())
end) 

@variables t q(t) r_tot(t)
@variables s(t) b(t) cr_tot(t) dec_s(t) r(t)

#is_negative_su(u,p,t) = u[s] < 0 | u[b] < 0 # cannot index s into u
hasnegativestate(u,p,t) = any(x -> x < 0, u)

mapdict(f,d) = Dict(pair.first => f(pair.second) for pair in d)
# d2 = Dict(:a => (1,2), :b => (3,4))
# mapdict(d2) do x
#     2 .* x
# end
    

function solve_system(sys_s,sys,p = p_straw; tspan = (0.0, 48.0))
    prob = ODEProblem(sys_s, p[:x0], tspan, p[:parms] )
    #sol = solve(prob, Rodas5(), isoutofdomain=is_negative_su);
    sol = solve(prob, Rodas5(), isoutofdomain=hasnegativestate);
    @assert all(sol[s] .+ sol[b] .+ sol[cr_tot] .≈ sum(x.second for x in p[:x0][1:2]))
    sol
end
#sol = gen_problem(systems[:simp]...);
sols = mapdict(systems) do system
    solve_system(system...);
end;
@assert all(isapprox.(sols[:phys](48.0)[1:3], sols[:simp](48.0), atol=0.01))


tmp_inspectsol()
    #plot(sol)
    plot(sol,vars=[s, b, cr_tot, dec_s, r_tot])
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




