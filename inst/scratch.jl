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

function ftmp()
    tmp = (Vector{ForwardDiff.Dual{ForwardDiff.Tag{Turing.Core.var"#f#3"{DynamicPPL.TypedVarInfo{NamedTuple{(:pg, :pl), Tuple{DynamicPPL.Metadata{Dict{AbstractPPL.VarName{:pg, Tuple{Tuple{Int64}}}, Int64}, Vector{Distribution{Univariate, Continuous}}, Vector{AbstractPPL.VarName{:pg, Tuple{Tuple{Int64}}}}, Vector{Float64}, Vector{Set{DynamicPPL.Selector}}}, DynamicPPL.Metadata{Dict{AbstractPPL.VarName{:pl, Tuple{Tuple{Int64}}}, Int64}, Vector{Distribution{Univariate, Continuous}}, Vector{AbstractPPL.VarName{:pl, Tuple{Tuple{Int64}}}}, Vector{Float64}, Vector{Set{DynamicPPL.Selector}}}}}, Float64}, DynamicPPL.Model{var"#73#74", (:qobsg, :r_totobsg, :qobsl, :r_totobsl, :T), (:T,), (), Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Type{Float64}}, Tuple{Type{Float64}}}, DynamicPPL.Sampler{NUTS{Turing.Core.ForwardDiffAD{40}, (), AdvancedHMC.DiagEuclideanMetric}}, DynamicPPL.DefaultContext}, Float64}, Float64, 4}}, Vector{Float64})

    tmp2 = (Vector{ForwardDiff.Dual{ForwardDiff.Tag{Turing.Core.var"#f#3"{DynamicPPL.TypedVarInfo{NamedTuple{(:pg, :pl), Tuple{DynamicPPL.Metadata{Dict{AbstractPPL.VarName{:pg, Tuple{Tuple{Int64}}}, Int64}, Vector{Distribution{Univariate, Continuous}}, Vector{AbstractPPL.VarName{:pg, Tuple{Tuple{Int64}}}}, Vector{Float64}, Vector{Set{DynamicPPL.Selector}}}, DynamicPPL.Metadata{Dict{AbstractPPL.VarName{:pl, Tuple{Tuple{Int64}}}, Int64}, Vector{Distribution{Univariate, Continuous}}, Vector{AbstractPPL.VarName{:pl, Tuple{Tuple{Int64}}}}, Vector{Float64}, Vector{Set{DynamicPPL.Selector}}}}}, Float64}, DynamicPPL.Model{var"#75#76", (:qobsg, :r_totobsg, :qobsl, :r_totobsl, :T), (:T,), (), Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Type{Float64}}, Tuple{Type{Float64}}}, DynamicPPL.Sampler{NUTS{Turing.Core.ForwardDiffAD{40}, (), AdvancedHMC.DiagEuclideanMetric}}, DynamicPPL.DefaultContext}, Float64}, Float64, 4}}, Vector{Float64})
end

function tmpf()
    ftmp = function()
        if rand() > 0.5 
            return nothing
        end
        [1,2,3]
    end
    @code_warntype ftmp()
    # bettern return zeros(3) instead of nothing
end

function _negate(f)
    local fc = f
    function negate(args...; kwargs...)
        -fc(args...; kwargs...)
    end
end

nsum = _negate(sum)




