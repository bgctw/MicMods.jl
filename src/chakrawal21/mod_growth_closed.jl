# integrated solution to the unlimited growth phase
# 
# TODO add microbial turnover tvr_b = kd*b to closed form

function chak21_growth_closed(;tend = 10.0)
    #using ModelingToolkit, DifferentialEquations, Plots
    @parameters t ks ksm Y HS HB 
    # t in hours
    @variables s(t) b(t) cr_tot(t) r(t) r_tot(t) q(t) 
    # masses per soilmass in mol/g
    HG = HS - Y * HB # can take out of system because does not involve t
    # see S5 for strawN
    #λ = 0.9
    ksm = estimate_ksm(ks) 
    parms = OrderedDict(
        ks => 0.15, 
        #ksm => 0.0, # updated down by ks and Y
        Y => 0.72, 
        HB => -492, HS => -469,
    )
    #parms[ksm] = (1 - parms[Y]) * parms[ks] * (1/λ - 1) # λ = r_gr/(r_gr + ksm)
    searchranges_p = OrderedDict(
        #ks => (parms[ks]/50, parms[ks]*50, :lognormal),   #ks 
        #ksm => (parms[ksm]/50, parms[ksm]*50, :lognormal),   #ks 
        # Y unconstrained in growth model - fix it
        #Y => (0.1, 0.8, :normal),  #Y
    )
    x0 = OrderedDict(
        s => 4.17e-5*1e6,
        b => 1.14e-5*1e6, 
        r => 0.4,
        cr_tot => 0,
    )
    searchranges_u0 = OrderedDict(
        #s => ((0.1 .+ (-0.01,+0.01))..., :normal),  #s0 is fixed
        b => (1.0, 100.0, :lognormal),  #b0
        r => (1e-8, 1-1e-8, :uniform), 
        # cr_tot is fixed 0
    )
    # parmspos = map(enumerate(parms)) do (i,(k,v)); k => i  end
    # statespos = map(enumerate(x0)) do (i,(k,v)); k => i  end
    # function predout(u,p)
    #     local s,b,cr_tot,r, 
    #         ks,ksm,Y,HB,HS,
    #         u1,u2,dec_s,r_gr,r_m,r_tot,HG,q
    #     # here we know the order in u and p 
    #     s,b,r,cr_tot = u
    #     ks,ksm,Y,HB,HS = p
    #     u1 = b * ks * r
    #     u2 = b * ksm # *(1-r) # uncoupled resp also for growing P components
    #     dec_s = u1 + u2
    #     r_gr = (1-Y) * u1
    #     r_m = u2
    #     r_tot = r_gr + r_m
    #     HG = HS - Y * HB 
    #     q = -HG * u1 - HS * u2
    #     (dec_s = dec_s, r_tot = r_tot, q = q, u1 = u1, u2 = u2,
    #     r_gr = r_gr, r_m = r_m)
    # end
    function predcum(t,u0,p)
        local s0,x0,r0, 
            ks,ksm,Y,HB,HS,
            HG
        s0,x0,r0,_ = u0 # x is biomass b
        ks,Y,HB,HS = p
        ksm = estimate_ksm(ks)
        # all quantities involving t are vectors, others scalars
        μ = ks * Y
        emt = exp.(μ .* t)
        r0_emt = r0 .* emt
        #r = r0 ./ ((1-r0) .* emt .+ r0)
        r = r0_emt ./ ((1 - r0) .+ r0_emt)
        xr = (x0 * r0) .* emt
        x = (x0*(1-r0)) .+ xr
        a = x0*(1-r0)*ksm
        b = (1-Y)*ks + ksm
        r_tot =  a .+  b .* xr
        cr_tot = a .* t + b/μ .* (emt .- 1)
        HG = HS - Y * HB
        q = (-HS * x0 * (1-r0)) .- (HG * ks + HS * ksm) .* emt
        s = (s0+x0) .- x .- cr_tot
        SLVector(
            s = s, b = x, r = r, cr_tot = cr_tot, r_tot = r_tot, q = q)
    end
    (#prob = prob, syss = syss,
    #sys = chak21_growth,
    #x0 = x0, parms = parms, # from prob
    #parmspos = parmspos, statespos = statespos, 
    #states = states_map, observed = observed_map,
    #predout = predout,
    x0 = x0, parms = parms, # here cannot infer from prob
    predcum = predcum,
    searchranges_p = searchranges_p, searchranges_u0 = searchranges_u0,
    )
end

