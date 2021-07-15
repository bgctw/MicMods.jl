# kinetic respiration model with simplified coefficients

"""
    micfromcoef(params; λ = 0.9, YCO2 = 1.5)

Compute microbial parameters from coefficient regressions β0,β1,β2.
See Wutzler11 eq. (2,3,4)

Value: Numeric array: `[μmax, r0, x0, λ, YCO2]`
"""
function micfromcoef(params; λ = 0.9, YCO2 = 1.5)
    β0,β1,β2 = (p for p in params)    
    μmax = β2
    r0 = β1*(1-λ)/(β0 + β1*(1-λ))
    x0 = β1*λ*YCO2/(r0*μmax)
    [μmax, r0, x0, λ, YCO2]
end

"""
    Fit exponential equation with lognormal priors on coefficients.
"""
function kinresp_exp()
    function predout(t, params)
        β0,β1,β2 = (p for p in params)    
        q =  @. β0 + β1 * exp(β2 * t)
    end
    parms = [
        :β0 => 1.0, 
        :β1 => 1.0,  
        :β2 => 1.0, 
    ]
    searchranges_p = OrderedDict(
        :β0 => (0.1, 10, :lognormal),  
        :β1 => (0.1, 10, :lognormal),  
        :β2 => (0.1, 10, :lognormal),  
    )
    (parms = parms, searchranges_p = searchranges_p, predout=predout)
end

"""
    Fit the kinetic respiraiton model using microbial parameters
"""
function kinresp_mic()
    function predout(t, params)
        x0,μmax,r0,YCO2,λ = params
        #x0l,μmaxl,r0l,YCO2l,λ = (p in params)
        # x0 = exp(x0l)
        # μmax = exp(μmaxl)
        # r0 = logistic(r0l),
        # YCO2 = logistic(YCO2l)
        q = @. x0 * (1 - r0) * (1/λ - 1)*μmax/YCO2 + 
            x0 * r0 * 1/λ * μmax/YCO2 * exp(μmax * t)
    end
    parms = [
        :x0 => 140, 
        :μmax => 0.13,
        :r0 => 0.026,
        :YCO2 => 1.5,
        :λ => 0.9 # ratio between productive (growth associated) part and total of 
        # the specific respiration activity
    ]
    searchranges_p = OrderedDict(
        :x0 => (10, 10000, :lognormal),  
        :μmax => (0.001, 1.5, :lognormal), 
        :r0 => (0.005, 1.0, :uniform), 
        :YCO2 => (0.1, 0.8, :uniform), 
        :λ => (0.88, 0.92, :uniform), 
    )
    (parms = parms, searchranges_p = searchranges_p, predout=predout)
end

