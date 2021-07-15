# simplified model from Chakrawal21

function chak21_growth_system(;tend = 10.0)
    #using ModelingToolkit, DifferentialEquations, Plots
    @parameters t ks kd kmr Y HS HB #ksm
    # t in hours
    # masses per soilmass in mol/g
    @variables s(t) b(t) cr_tot(t) r(t) r_tot(t) q(t) dec_s(t) tvr_b(t)
    @variables mm_s(t) u1(t) u2(t) r_gr(t) r_m(t) db(t)
    D = Differential(t)
    HG = HS - Y * HB # can take out of system because does not involve t
    #λ = 0.9
    #ksm = (1-Y)*ks * (1/λ -1) # ksm = 0.11 ks
    ksm = estimate_ksm(ks) 
    eqs = [
        D(s) ~ -dec_s + tvr_b,
        db ~ Y * u1 - tvr_b, # store output
        D(b) ~ db,
        #D(r) ~ Y * u1/b * (1 - r),
        D(r) ~ Y * u1/b * (s/(s+kmr) - r),
        D(cr_tot) ~ r_tot,
        u1 ~ b * ks * r, #* s/(s + km),
        u2 ~ b * ksm * (1-r), #* s/(s + km),
        dec_s ~ u1 + u2,
        tvr_b ~ kd*b,
        r_gr ~ (1-Y) * u1,
        r_m ~ u2,
        r_tot ~ r_gr + r_m,
        q ~ -HG * u1 - HS * u2,
    ]
    @named chak21_growth = ODESystem(eqs)
    syss = structural_simplify(chak21_growth) 
    # see S5 for strawN
    parms = [
        ks => 0.15, 
        #ksm => 0.15/10, # uptake for P (maint comp), here maint resp
        #km => 2.62, 
        kd => 1.35e-2, 
        #kd => 0.0, 
        Y => 0.72, 
        HB => -492, HS => -469,
        kmr => 1e-9, # half-saturation low so that optimum r -> 1
    ]
    searchranges_p = OrderedDict(
        #ks and km only appear as um, need to fix one
        ks => (0.05, 1.0, :lognormal),   #ks 
        #ksm => (0.005, 1.0, :lognormal),   #ksm
        #kd => (1/200, 1/10, :lognormal),  #kd
        #kmr => (1, 40, :lognormal),  #kmr
        # Y unconstrained in growth model - fix it
        #Y => (0.1, 0.8, :normal),  #Y
    )
    x0 = [s => 4.17e-5*1e6, b => 1.14e-5*1e6, r => 0.4, cr_tot => 0]
    searchranges_u0 = OrderedDict(
        #s => ((0.1 .+ (-0.01,+0.01))..., :normal),  #s0 is fixed
        b => (1.0, 100.0, :lognormal),  #b0
        r => (1e-8, 1-1e-8, :uniform), 
    )
    parmspos = OrderedDict(Symbol(x) => indexof(x,parameters(syss)) for x in 
    (ks,kd,Y,HB,HS)) #kms from alpha
    statespos = OrderedDict(Symbol(x.val.f.name) => indexof(x,states(syss)) for x in 
    (s,b,cr_tot,r))
    # states_map = OrderedDict(Symbol(x.f.name) => x for x in states(syss))
    # observed_map = OrderedDict(Symbol(x.lhs.f.name) => x.lhs for x in observed(syss))
    prob = chak21_problem(syss, x0, parms; ti=range(0,tend,step = 1/4))
    # uses closure parmspos, λ
    function predout(u,p)
        # # without let or local unstable reassingning in closure
        # let s,b,cr_tot,r, ks,ksm,kd,Y,HB,HS, r, r_tot, q, dec_s, tvr_b,
        #     u1, u2, r_gr, r_m, HG
        local s,b,cr_tot,r, ks,ksm,kd,Y,HB,HS, r, r_tot, q, dec_s, tvr_b,
            u1, u2, r_gr, r_m, HG
            s,b,cr_tot,r = (u[pos] for pos in values(statespos))
            ks,kd,Y,HB,HS = (p[pos] for pos in values(parmspos))
            ksm = ks * (1/λ -1) # ksm = 0.11 ks
            u1 = b * ks * r
            u2 = b * ksm * (1-r)
            dec_s = u1 + u2
            tvr_b = kd*b
            r_gr = (1-Y) * u1
            r_m = u2
            r_tot = r_gr + r_m
            HG = HS - Y * HB 
            q = -HG * u1 - HS * u2
            (dec_s = dec_s, r_tot = r_tot, q = q, u1 = u1, u2 = u2, tvr_b = tvr_b,
            r_gr = r_gr, r_m = r_m)
        #end
    end
    (prob = prob, syss = syss,
    #sys = chak21_growth,
    #x0 = x0, parms = parms, # from prob
    #parmspos = parmspos, statespos = statespos, 
    #states = states_map, observed = observed_map,
    predout = predout,
    searchranges_p = searchranges_p, searchranges_u0 = searchranges_u0,
    )
end

