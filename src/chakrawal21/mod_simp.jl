# simplified model from Chakrawal21
# assuming r=1 (not declining with substrate)

function chak21_simp_system()
    #using ModelingToolkit, DifferentialEquations, Plots
    @parameters t ks km kd Y HS HB HG
    # t in hours
    # masses per soilmass in mol/g
    @variables s(t) b(t) cr_tot(t) r_tot(t) q(t) dec_s(t) tvr_b(t)
    D = Differential(t)
    HG = HS - Y * HB # can take out of system because does not involve t
    eqs = [
        dec_s ~ (ks*s*b)/(s + km),
        tvr_b ~ kd*b,
        D(s) ~ -dec_s + tvr_b,
        D(b) ~ Y * dec_s - tvr_b,
        r_tot ~ (1-Y) * dec_s,
        D(cr_tot) ~ r_tot,
        q ~ -HG * dec_s,
        ]
    @named chak21_simp = ODESystem(eqs)
    syss = structural_simplify(chak21_simp) 
    # see S5
    # parms = [ks => 0.15, km => 2.62e-6, kd => 1.35e-2, Y => 0.72, 
    #         HB => -492, HS => -469]
    # #:x0 => [s => 4.17e-6, b => 1.14e-5, cr => 0] # twutz: wrong in paper
    # x0 = [s => 4.17e-5, b => 1.14e-5, cr_tot => 0] 
    # express in μmol instead of mol to avoid small numbers
    parms = [
        ks => 0.15, km => 2.62, kd => 1.35e-2, Y => 0.72, 
        HB => -492, HS => -469
    ]
    searchranges_p = OrderedDict(
        #ks => (0.005, 1.0, :lognormal),   #ks
        #ks => (0.005, 1.0, :uniform),   #ks
        #ks => (0.05, 1.0, :lognormal),   #ks
        ks => (0.05, 1.0, :uniform),   #ks
        km => (0.1, 10.0, :lognormal),    #km
        kd => (1/200, 1/10, :lognormal),  #kd
        Y => (0.1, 0.8, :normal),  #Y
    )
    x0 = [s => 4.17e-5*1e6, b => 1.14e-5*1e6, cr_tot => 0] 
    searchranges_u0 = OrderedDict(
        #s => (0.1, 1/10, :lognormal),  #s0 is fixed
        b => (1.0, 100.0, :lognormal),  #b0
    )
    parmspos = OrderedDict(Symbol(x) => indexof(x,parameters(syss)) for x in 
        (ks,km,kd,Y,HB,HS))
    #states_map = OrderedDict(x.f.name => x for x in states(syss))
    # observed_map = OrderedDict(Symbol(x.lhs.f.name) => x.lhs for x in observed(syss))
    prob = chak21_problem(syss, x0, parms)
    # uses closure parmspos
    function predout(u,p)
        local s,b,cr, ks, km, kd, Y, HB, HS,
            dec_s, r_tot, HG, q
        s,b,cr = u
        ks, km, kd, Y, HB, HS = (p[pos] for pos in values(parmspos))
        dec_s = (ks*s*b)/(s + km)
        r_tot = (1-Y) * dec_s
        HG = HS - Y * HB # can take out of system because does not involve t
        q = -HG * dec_s
        (dec_s = dec_s, r_tot = r_tot, q = q)
    end
    (prob = prob, syss = syss, 
        # sys = chak21_simp, x0 = x0, parms = parms, 
        # parmspos = parmspos, states = states_map, observed = observed_map,
        searchranges_p = searchranges_p, searchranges_u0 = searchranges_u0,
        predout = predout)
end

