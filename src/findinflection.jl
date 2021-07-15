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
find_inflection = function(x,y)
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