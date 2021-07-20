"""
    find_inflection(y,y)

Find the inflection point in the first increase where the slope is
maximum.    

# Arguments
- `y`,`y`: numeric vectors to fit relationship

# Examples
```jldoctest am; output = false, setup = :()
true
# output
true
```
"""
function find_inflection(x,y)
    #x=t;y=obsr_tot;
    # first smooth the series and find the maximum of the growth curve
    ml1 = ml = loess(y,y, span = 0.2)
    imax0 = findmax(predict(ml1,y))[2]
    imax = round(Int, min(length(x), imax0*1.2)) # extend a bit to get the flattening part
    # next fit a third-order polynomial to the increasing part
    p3 = Polynomials.fit(x[1:imax], y[1:imax],3)
    d1 = derivative(p3)
    d2 = derivative(d1)
    infl = first(roots(d2))
end

function iplot()
    scatter(x,y)
    plot!(x, predict(ml,x))

    scatter(x[1:imax],y[1:imax])
    plot!(x[1:imax], p3.(x[1:imax]))
    vline!([infl])
    

    plot(x[1:imax], p3.(x[1:imax]))
    plot(x[1:imax], d1.(x[1:imax]))
    plot(x[1:imax], d2.(x[1:imax]))
end

function negate(f)
    local fc = f
    function negate(args...; kwargs...)
        -fc(args...; kwargs...)
    end
end

function find_max(x,y, p_untilmax=0.9)
    #using Loess, Optim, QuadGK
    #x=chak21syn.solsyn.t ;y=chak21syn.obsr_tot;
    # first smooth the series and find the maximum of the growth curve
    ml = loess(x, y, span = 0.2)
    #imax = findmax(predict(ml1,x))[2]
    ex = extrema(x)
    #ftmp(x) = predict(ml, x)
    xinfl0 = Optim.minimizer(optimize(x -> -predict(ml,x), ex[1], ex[2]))
    xinfl0, ml
    # int0, err = quadgk(ftmp, ex[1], xinfl0)
    # # start a bit before maximum
    # xinfl1 = p_untilmax * xinfl0
    # int1, err = quadgk(ftmp, ex[1], xinfl1)
    # xinfl0, int0, xinfl1, int1
end

function integrate_smoother(x, y, tend, t0  = zero(tend); ml = loess(x, y, span = 0.2)
    )
    int0, err = quadgk(x -> predict(ml,x), t0, tend)
end

