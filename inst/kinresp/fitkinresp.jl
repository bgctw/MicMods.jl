using MicMods
using Turing, Optim
using OrderedCollections

#systemt = kinresp_exp()
systemt = kinresp_mic()
tstep = 1.0:0.25:10.0
nobs = length(tstep)
sr = OrderedDict(Iterators.take(systemt.searchranges_p, 3))
popt_names = keys(sr)
popt_true = [OrderedDict(systemt.parms)[pname] for pname in popt_names]
npar = length(popt_true)
# include initial conditions in parameters to estimate
pall = vcat(popt_true, [x.second for x in systemt.parms[(npar+1):end]])
qobs_true = systemt.predout(tstep,pall)
σ_qobs0 = 0.02*(maximum(qobs_true)-minimum(qobs_true))
δ = 0.5 # assume fixed and known, see discussion in Wutzler11
σ_qobs = σ_qobs0 .* abs.(qobs_true).^δ
qobs = qobs_true .+ σ_qobs .* randn(nobs)

# using Plots
scatter(tstep, qobs, label = "Observations",
    ylab = "Respiration (g/g)", xlab="Time (hr)",
    legend=:topleft,)
plot!(tstep, qobs_true, label = "True process")
#plot!(tstep, qpred, label = "Initial fit")
savefig("tmp.svg")


# initial values of microbial parameters from log(y-b0) ~ log(b1) + b2*t
using GLM: lm, @formula, coef
using DataFrames
b0 = 0.9*minimum(qobs)
df = DataFrame(ly = log.(qobs .- b0), tstep = tstep);
#plot(df.tstep, df.ly)
ols = lm(@formula(ly ~ tstep), df);
b1l, b2 = coef(ols); b1 = exp(b1l)
qpred = @. b0 + b1*exp(b2*tstep)
p0 = micfromcoef([b0,b1,b2])
function tmp.depr()
    # initial values of uncertainty from: log(res) ~ log(σ0) + δ qpred
    df = DataFrame(logres = log.(abs.(qpred .- qobs)), tstep = tstep, qpred = qpred);
    plot(df.tstep, df.logres)
    ols = lm(@formula(logres ~ qpred), df);
    σ00l, δ0 =  coef(ols); σ00 = exp(σ00l)
    plot(df.logres, σ00l .+ δ0 .* qpred)
    plot(exp.(df.logres), σ00 .* qpred.^δ0)
end
#σ00 = median(abs.(qpred .- qobs))
σ00 = median(abs.(qpred .- qobs)./qpred.^δ)

function optim_turing()
    #qobs, popt_true, sr, systemt.parms systemt.predout, step, σ00, δ
    @model function fitq(qobs, ::Type{T} = Float64) where {T}
        npar = length(popt_true)
        p = Vector{T}(undef, npar)
        for (i,r) = enumerate(sr)
             p[i] ~ gettruncdist(r.second...)
        end
        pall = vcat(p, [x.second for x in systemt.parms[(npar+1):end]])
        qpred = systemt.predout(tstep,pall)
        σ0 ~ LogNormal(log(σ00), log(2))
        #δ ~ Normal(0.5, 0.5) # assume given - same across many experiments
        σ_pred = σ0 .* abs.(qpred).^δ
        qobs ~ MvNormal(qpred, σ_qobs)
    end
    #p0 from linear fit above, popt_names
    model_func = fitq(qobs)
    chn = chn0 = sample(model_func, NUTS(), 10, init_theta = vcat(p0, σ00))
    chn = chn0 = sample(model_func, NUTS(),  MCMCThreads(), 400, 3, init_theta = vcat(p0, σ00))
    #chn = chn0 = sample(model_func, NUTS(),  MCMCThreads(), 800, 3, init_theta = Array(chn0)[end,:])
    chn = replacenames(chn0, Dict("p[$i]" => pname for (i,pname) in enumerate(popt_names)))
    function interactive_plotchains()
        using StatsPlots
        plot(chn)
        namesingroup(chn0, :p)
        meanplot(chn)
        corner(chn)

        plotchains(get(chn, [:log_density]))
        
        size(Array(chn))
        plot(get(chn, :acceptance_rate)[1])
        plot(get(chn, :hamiltonian_energy)[1])
        Array(chn0)[end,:]
        
    end
end








