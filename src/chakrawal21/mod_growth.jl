# simplified model from Chakrawal21

function chak21_growth_system()
    #using ModelingToolkit, DifferentialEquations, Plots
    @parameters t ks ksm kd Y HS HB 
    # t in hours
    # masses per soilmass in mol/g
    @variables s(t) b(t) cr_tot(t) r(t) r_tot(t) q(t) dec_s(t) tvr_b(t)
    @variables mm_s(t) u1(t) u2(t) r_gr(t) r_m(t) db(t)
    D = Differential(t)
    HG = HS - Y * HB # can take out of system because does not involve t
    eqs = [
        D(s) ~ -dec_s + tvr_b,
        db ~ Y * u1 - tvr_b, # store output
        D(b) ~ db,
        D(cr_tot) ~ r_tot,
        D(r) ~ Y * u1/b * (1 - r),
        #D(r) ~ Y * u1/b * (s/(s+kmr) - r),
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
        ksm => 0.15, # uptake for P (maint comp) equal that for U (growth)
        #km => 2.62, 
        #kd => 1.35e-2, 
        kd => 0.0, 
        Y => 0.72, 
        HB => -492, HS => -469,
        #kmr => 1e-9, # half-saturation low so that optimum r -> 1
    ]
    searchranges_p = OrderedDict(
        ks => (0.005, 1.0),   #ks
        ksm => (0.005, 1.0),   #ksm
        #kd => (1/200, 1/10),  #kd
        Y => (0.1, 0.8),  #Y
    )
    x0 = [s => 4.17e-5*1e6, b => 1.14e-5*1e6, cr_tot => 0, r => 0.5]
    searchranges_u0 = OrderedDict(
        #s => (0.1, 1/10),  #s0 is fixed
        b => (1.0, 100.0),  #b0
        r => (1e-8, 1-1e-8), 
    )
    parmspos = OrderedDict(Symbol(x) => indexof(x,parameters(syss)) for x in 
    (ks,ksm,kd,Y,HB,HS))
    statespos = OrderedDict(Symbol(x.val.f.name) => indexof(x,states(syss)) for x in 
    (s,b,cr_tot,r))
    states_map = OrderedDict(Symbol(x.f.name) => x for x in states(syss))
    observed_map = OrderedDict(Symbol(x.lhs.f.name) => x.lhs for x in observed(syss))
    prob = chak21_problem(syss, x0, parms)
    function predout(u,p,system)
        b,cr = u
        ks,ksm,kd,Y,HB,HS = (p[pos] for pos in values(system.parmspos))
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
    end
    (prob = prob, syss = syss, sys = chak21_growth, x0 = x0, parms = parms,
    parmspos = parmspos, statespos = statespos, 
    states = states_map, observed = observed_map,
    predout = predout,
        searchranges_p = searchranges_p, searchranges_u0 = searchranges_u0)
end

