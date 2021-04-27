# simplified model from Chakrawal21

function chak21_phys_system()
    #using ModelingToolkit, DifferentialEquations, Plots
    @parameters t ks km kd Y HS HB kmr ksm
    # t in hours
    # masses per soilmass in mol/g
    @variables s(t) b(t) cr_tot(t) r_tot(t) q(t) dec_s(t) tvr_b(t)
    @variables mm_s(t) u1(t) u2(t) r(t) r_gr(t) r_m(t) db(t)
    D = Differential(t)
    HG = HS - Y * HB # can take out of system because does not involve t
    eqs = [
        D(s) ~ -dec_s + tvr_b,
        db ~ Y * u1 - tvr_b, # store output
        D(b) ~ db,
        D(cr_tot) ~ r_tot,
        D(r) ~ Y * u1/b * (s/(s+kmr) - r),
        u1 ~ b * ks * r * s/(s + km),
        u2 ~ b * ksm * (1-r) * s/(s + km),
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
    parms = [ks => 0.15, km => 2.62, kd => 1.35e-2, Y => 0.72, 
            HB => -492, HS => -469,
            kmr => 1e-9, # half-saturation low so that optimum r -> 1
            ksm => 0.15, # uptake for P (maint comp) equal that for U (growth)
            ]
    x0 = [s => 4.17e-5*1e6, b => 1.14e-5*1e6, cr_tot => 0, r => 1]
    prob = chak21_problem(syss, x0, parms)
    (prob = prob, syss = syss, sys = chak21_phys, x0 = x0, parms = parms)
end

