# simplified model from Chakrawal21

function chak21_simp_system()
    #using ModelingToolkit, DifferentialEquations, Plots
    @parameters t ks km kd Y HS HB
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
    chak21_simp_s = structural_simplify(chak21_simp) # omit r,q
    p_straw = Dict(
        :parms => [ks => 0.15, km => 2.62e-6, kd => 1.35e-2, Y => 0.72, 
            HB => -492, HS => -469],
        #:x0 => [s => 4.17e-6, b => 1.14e-5, cr => 0] # twutz: wrong in paper
        :x0 => [s => 4.17e-5, b => 1.14e-5, cr_tot => 0]
    ) # see S5
    chak21_simp_s, chak21_simp, p_straw
end

