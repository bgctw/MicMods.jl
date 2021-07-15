test_tspan = function()
    #check whether u[1] relates to t=0 or tspan[1]
    using DifferentialEquations
    f(u,p,t) = 1.01*u
    u0 = 1/2
    tspan = (0.0,1.0)
    prob = ODEProblem(f,u0,tspan)
    sol0 = solve(prob)

    prob1 = ODEProblem(f,u0,tspan .+ 0.5)
    sol1 = solve(prob1)

    sol1[1] #t = 0.5
    sol1(0) #t = 0
end

