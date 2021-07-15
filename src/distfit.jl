function truncnormalfromrange(lower,upper)
    σ = (upper - lower)/1.96/2
    μ = lower + 1.96σ
    Truncated(Normal(μ,σ),lower,upper)
    #Normal(μ,σ)
end
function trunclognormalfromrange(lower,upper)
    loglower = log(lower)
    logupper = log(upper)
    σ = (logupper - loglower)/1.96/2
    μ = loglower + 1.96σ
    Truncated(LogNormal(μ,σ), lower, upper)
    #LogNormal(μ,σ)
end
function logitnormalfromrange(lower,upper)
    logitlower = logit(lower)
    logitupper = logit(upper)
    σ = (logitupper - logitlower)/1.96/2
    μ = logitlower + 1.96σ
    LogitNormal(μ,σ)
    #LogitNormal(μ,σ)
end
function trunclogitnormalfromrange(lower,upper)
    logitlower = logit(lower)
    logitupper = logit(upper)
    σ = (logitupper - logitlower)/1.96/2
    μ = logitlower + 1.96σ
    Truncated(LogitNormal(μ,σ), lower, upper)
    #LogitNormal(μ,σ)
end
function gettruncdist(lower, upper, dist)
    dist == :normal && return(truncnormalfromrange(lower,upper))
    dist == :lognormal && return(trunclognormalfromrange(lower,upper))
    dist == :logitnormal && return(trunclogitnormalfromrange(lower,upper))
    dist == :uniform && return(Uniform(lower,upper))
    error("unknown distribution label: $dist")
end