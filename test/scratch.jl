using DifferentialEquations

intro1 = function()
    f(u,p,t) = 1.01*u
    u0 = 1/2
    tspan = (0.0,1.0)

    prob = ODEProblem(f,u0,tspan)
    using Plots
    plot(sol,linewidth=5,title="Solution to the linear ODE with a thick line",
        xaxis="Time (t)",yaxis="u(t) (in μm)",label="My Thick Line!") # legend=false
    plot!(sol.t, t->0.5*exp(1.01t),lw=3,ls=:dash,label="True Solution!")
end

function mtk1
    using ModelingToolkit

    @variables t x(t) RHS(t)  # independent and dependent variables
    @parameters τ       # parameters
    D = Differential(t) # define an operator for the differentiation w.r.t. time
    
    # your first ODE, consisting of a single equation, indicated by ~
    @named fol_separate = ODESystem([ RHS  ~ (1 - x)/τ,
                                      D(x) ~ RHS ])
    
    using DifferentialEquations: solve
    using Plots: plot
    
    prob = ODEProblem(structural_simplify(fol_separate), [x => 0.0], (0.0,10.0), [τ => 3.0])
    sol = solve(prob)
    plot(sol, vars=[x,RHS])
    #sol[RHS] # get the derived variable as Array
    #hcat(sol.t, (sol[v] for v in (x,RHS))...)
end    

function lorenz!(du,u,p,t)
    σ,ρ,β = p
    du[1] = σ*(u[2]-u[1])
    du[2] = u[1]*(ρ-u[3]) - u[2]
    du[3] = u[1]*u[2] - β*u[3]
end
u0 = [1.0,0.0,0.0]
p = (10,28,8/3) # we could also make this an array, or any other type!
tspan = (0.0,100.0)
prob = ODEProblem(lorenz!,u0,tspan,p)
sol = solve(prob)


function interactive()
    using Plots
    plot(sol)
end

function mtk_example()
    #using ModelingToolkit
    @parameters t σ ρ β
    @variables x(t) y(t) z(t)
    D = Differential(t)
    eqs = [D(x) ~ σ*(y-x),
        D(y) ~ x*(ρ-z)-y,
        D(z) ~ x*y - β*z]
    de = ODESystem(eqs,t,[x,y,z],[σ,ρ,β])
    sol = solve(prob)
    plot(sol)
end

function catalyst_example1()
    #using Catalyst, Latexify
    rs = @reaction_network begin
    c1, X --> 2X
    c2, X --> 0
    c3, 0 --> X
    end c1 c2 c3
    #Graph(rs)
    p = (1.0,2.0,50.) # [c1,c2,c3]
    tspan = (0.,4.)
    u0 = [5.]         # [X]
    #osys = convert(ODESystem,rs)
    #latexify(osys)
    # solve ODEs
    oprob = ODEProblem(rs, u0, tspan, p)
    osol  = solve(oprob, Tsit5())
    #plot(osol)
    equations(osys)

    # solve for Steady-States
    #ssprob = SteadyStateProblem(rs, u0, p)
    ssprob = SteadyStateProblem(osys, u0, p)
    sssol  = solve(ssprob, SSRootfind())
end

function simpmod_example()
    #using ModelingToolkit, DifferentialEquations, Plots
    @parameters t ks km kd Y HS HB
    # t in hours
    # masses per soilmass in mol/g
    @variables s(t) b(t) cr(t) r_tot(t) q(t) dec_s(t) tvr_b(t) 
    D = Differential(t)
    HG = HS - Y * HB # can take out of system because does not involve t
    eqs = [
        dec_s ~ (ks*s*b)/(s + km),
        tvr_b ~ kd*b,
        D(s) ~ -dec_s + tvr_b,
        D(b) ~ Y * dec_s - tvr_b,
        r_tot ~ (1-Y) * dec_s,
        D(cr) ~ r_tot,
        q ~ -HG * dec_s,
        ]
    de = ODESystem(eqs; name=:simpmod)
    des = structural_simplify(de) # omit r,q
    p_straw = Dict(
        :parms => [ks => 0.15, km => 2.62e-6, kd => 1.35e-2, Y => 0.72, 
            HB => -492, HS => -469],
        #:x0 => [s => 4.17e-6, b => 1.14e-5, cr => 0] # twutz: wrong in paper
        :x0 => [s => 4.17e-5, b => 1.14e-5, cr => 0]
    ) # see S5
    prob = ODEProblem(des, 
        #[s => 1.0, b => 1.0, r => 0.0, q => 0.0],
        p_straw[:x0],
        (0.0, 48.0),
        #[ks => 0.5, km => 0.5, kd => 1/20, Y => 0.5, HG => -10]
        p_straw[:parms]
    )
    #sol = solve(prob);
    #sol = solve(prob,AutoTsit5(Rodas5()));
    #sol = solve(prob,lsoda()); # does not work with DAE
    sol = solve(prob,Rodas5());
    #sol[r_tot]
    #plot(sol)
    plot(sol,vars=[s, b, cr, dec_s, r_tot])
    # glucose taken up within 16 hours, but stored in microbial biomass
    # instead of respired
    plot(sol,vars=[q], ylim = (0, 800e-6)) 
    plot(sol,vars=[r_tot], ylim = (0,2e-6))
    # in the model resp and q only differ by factor (1-Y)/HG 
end


function test_min()
    # using minimum function in equation works nicely
    #using ModelingToolkit, DifferentialEquations, Plots
    @parameters t ks km kd Y HG cn_s cn_b imm_N
    @variables s(t) b(t) r(t) rO(t) q(t) dec_s(t) tvr_b(t) 
    @variables syn_c_pot(t) syn_n_pot(t) syn(t)
    D = Differential(t)
    HG = HS - Y * HB
    eqs = [
        dec_s ~ (ks*s*b)/(s + km),
        tvr_b ~ kd*b,
        syn_c_pot ~ Y * dec_s,
        syn_n_pot ~ (dec_s/cn_s + imm_N) * cn_b,
        syn ~ min(syn_c_pot, syn_n_pot), # note the min function
        D(s) ~ -dec_s + tvr_b,
        D(b) ~ +syn - tvr_b,
        r ~ dec_s - syn,
        rO ~ r - (1 - Y) * dec_s,
        q ~ -HG * dec_s,
        ]
    de = ODESystem(eqs)
    des = structural_simplify(de) # omit r,q
    prob = ODEProblem(des, 
        #[s => 1.0, b => 1.0, r => 0.0, q => 0.0],
        [s => 100.0, b => 1.0],
        (0.0, 50.0),
        [ks => 0.5, km => 0.5, kd => 1/10, Y => 0.5, HG => -10,
        cn_b => 8, cn_s => 20,
        imm_N => 0.05];
        #jac=true # takes longer but ok here
    )
    sol = solve(prob);
    #plot(sol)
    plot(sol,vars=[r, rO, syn, syn_c_pot, syn_n_pot]) # r and q are tracked
    plot(sol,vars=[s, b, r, q, dec_s]) # r and q are tracked
end

function test_local()
    x = "bla"
    function inner(arr)
        local x,y
        x, y = arr
        x
    end
    ans = inner( (2,4))
end


