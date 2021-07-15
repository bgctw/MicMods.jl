# physiological model from Chakrawal21

function chak21_phys_system()
    #using ModelingToolkit, DifferentialEquations, Plots
    @parameters t ks km kd Y HS HB kmr #ksm
    # t in hours
    # masses per soilmass in mol/g
    @variables s(t) b(t) cr_tot(t) r_tot(t) q(t) dec_s(t) tvr_b(t)
    @variables mm_s(t) u1(t) u2(t) r(t) r_gr(t) r_m(t) db(t)
    D = Differential(t)
    HG = HS - Y * HB # can take out of system because does not involve t
    # assume constant λ = ks/(ks+ksm): ratio between productive (growth associated) part and 
    # total of the specific respiration activity
    ksm = estimate_ksm(ks) 
    eqs = [
        D(s) ~ -dec_s + tvr_b,
        db ~ Y * u1 - tvr_b, # store output
        D(b) ~ db,
        D(cr_tot) ~ r_tot,
        D(r) ~ Y * u1/b * (s/(s+kmr) - r),
        mm_s ~ s/(s + km),
        u1 ~ b * ks * r * mm_s,
        u2 ~ b * ksm * (1-r) * mm_s,
        dec_s ~ u1 + u2,
        tvr_b ~ kd*b,
        r_gr ~ (1-Y) * u1,
        r_m ~ u2,
        r_tot ~ r_gr + r_m,
        q ~ -HG * u1 - HS * u2,
        ]
    @named chak21_phys = ODESystem(eqs)
    syss = structural_simplify(chak21_phys) 
    # see S5 for strawN
    parms = [ks => 0.15, km => 2.62, kd => 1.35e-2, 
            Y => 0.72, # yield: 
            HB => -492, HS => -469,
            kmr => 1e-9, # half-saturation low so that optimum r -> 1
            #ksm => 0.10, # uptake for P (maint comp): lower than for U (growth)
            ]
    searchranges_p = OrderedDict(
        ks => (0.05, 1.0, :uniform),   #ks
        km => (0.1, 10.0, :lognormal),    #km
        kd => (1/200, 1/10, :lognormal),  #kd
        Y => (0.1, 0.8, :normal),  #Y
        kmr => (1, 40, :lognormal),  #Y1e-9, # half-saturation low so that optimum r -> 1
    )
    x0 = [s => 4.17e-5*1e6, b => 1.14e-5*1e6, cr_tot => 0, r => 0.4]
    searchranges_u0 = OrderedDict(
        s => (0.1, 1/10, :lognormal),  #s0 is fixed
        b => (1.0, 100.0, :lognormal),  #b0 assume fixed by growth model
        r => (1e-5, 1-1e-5, :uniform),  #r0 assume fixed by growth model
    )
    parmspos = OrderedDict(Symbol(x) => indexof(x,parameters(syss)) for x in 
        (ks,km,kd,Y,HB,HS,kmr))
    #states_map = OrderedDict(x.f.name => x for x in states(syss))
    prob = chak21_problem(syss, x0, parms)
    # parmspos and λ from closure
    function predout(u,p)
        local s,b,r,cr, 
            ks, km, kd, Y, HB, HS, kmr,
            u1,u2,ksm,
            dec_s, r_gr, r_m, r_tot, HG, q
        s,b,cr,r = u
        ks, km, kd, Y, HB, HS, kmr = (p[pos] for pos in values(parmspos))
        ksm = ks * (1/λ -1) # ksm = 0.11 ks
        u1 = b * ks * r * s/(s + km)
        u2 = b * ksm * (1-r) * s/(s + km)
        dec_s = u1 + u2
        r_gr = (1-Y) * u1
        r_m = u2
        r_tot = r_gr + r_m
        HG = HS - Y * HB # can take out of system because does not involve t
        q = -HG * u1 - HS * u2
        (dec_s = dec_s, r_tot = r_tot, q = q)
    end
    (prob = prob, syss = syss, 
        # sys = chak21_simp, x0 = x0, parms = parms, 
        # parmspos = parmspos, states = states_map, observed = observed_map,
        searchranges_p = searchranges_p, searchranges_u0 = searchranges_u0,
        predout = predout)
end

estimate_ksm(ks, λ = 0.9) = ks * (1/λ -1) # ksm = 0.11 ks
