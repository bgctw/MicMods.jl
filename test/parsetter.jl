using Test
using MicMods
using ModelingToolkit
@parameters t ks km kd Y HS HB HG
@variables s(t) b(t) r(t) cr_tot(t) r_tot(t) q(t) dec_s(t) tvr_b(t)

using StaticArrays
ps = ParSetter( SA[ks,km,kd], SA[s,b,r], SA[ks,km], SA[b])

using DualNumbers

@testset "setting parameters SA" begin
    p = SA[1.1,2.2,3.3]
    u0 = SA[11.1,22.2,33.3]
    popt = getpopt(ps,p,u0)
    @test popt == SVector{3,Num}(p[1:2]...,u0[2])
    #
    po,u0o = setpu(ps, popt, p, u0)
    @test typeof(po) == typeof(p)
    @test po[3] == p[3]
    @test po[1:2] == popt[1:2]
    @test typeof(u0o) == typeof(u0)
    @test u0o[[1,3]] == u0[[1,3]]
    @test u0o[2] == popt[3]
    #
    popt2 = popt .* 200
    po,u0o = setpu(ps, popt2, p, u0)
    @test typeof(po) == typeof(p)
    @test po[3] == p[3]
    @test po[1:2] == popt2[1:2]
    @test typeof(u0o) == typeof(u0)
    @test u0o[[1,3]] == u0[[1,3]]
    @test u0o[2] == popt2[3]
end;

@testset "setting parameters Vector" begin
    p = [1.1,2.2,3.3]
    u0 = [11.1,22.2,33.3]
    popt = getpopt(ps,p,u0)
    @test popt == SVector{3,Num}(p[1:2]...,u0[2])
    #
    po,u0o = setpu(ps, popt, p, u0)
    @test typeof(po) == typeof(p)
    @test po[3] == p[3]
    @test po[1:2] == popt[1:2]
    @test typeof(u0o) == typeof(u0)
    @test u0o[[1,3]] == u0[[1,3]]
    @test u0o[2] == popt[3]
    #
    popt2 = popt .* 200
    po,u0o = setpu(ps, popt2, p, u0)
    @test typeof(po) == typeof(p)
    @test po[3] == p[3]
    @test po[1:2] == popt2[1:2]
    @test typeof(u0o) == typeof(u0)
    @test u0o[[1,3]] == u0[[1,3]]
    @test u0o[2] == popt2[3]
end;

@testset "setting parameters DualNumbers" begin
    p = [1.1,2.2,3.3]
    u0 = [11.1,22.2,33.3]
    popt = getpopt(ps,p,u0)
    poptd = map(Dual, popt)
    typeof(poptd)
    eltype(poptd)
    #
    po,u0o = setpu(ps, poptd, p, u0)
    @test eltype(po) == eltype(poptd)
    @test convert(eltype(p), po[3]) == p[3]
    @test po[1:2] == poptd[1:2]
    @test eltype(u0o) == eltype(poptd)
    @test u0o[[1,3]] == u0[[1,3]]
    @test convert(eltype(u0),u0o[2]) == poptd[3]
    #
    poptd2 = poptd .* 200
    po,u0o = setpu(ps, poptd2, p, u0)
    @test eltype(po) == eltype(poptd2)
    @test convert(eltype(p), po[3]) == p[3]
    @test po[1:2] == poptd2[1:2]
    @test eltype(u0o) == eltype(poptd2)
    @test u0o[[1,3]] == u0[[1,3]]
    @test convert(eltype(u0),u0o[2]) == poptd2[3]
end;

@testset "setting parameters DualNumbers and SVector" begin
    p = SA[1.1,2.2,3.3]
    u0 = SA[11.1,22.2,33.3]
    popt = getpopt(ps,p,u0)
    poptd = map(Dual, popt)
    typeof(poptd)
    eltype(poptd)
    #
    po,u0o = setpu(ps, poptd, p, u0)
    @test eltype(po) == eltype(poptd)
    @test convert(eltype(p), po[3]) == p[3]
    @test po[1:2] == poptd[1:2]
    @test eltype(u0o) == eltype(poptd)
    @test u0o[[1,3]] == u0[[1,3]]
    @test convert(eltype(u0),u0o[2]) == poptd[3]
    #
    poptd2 = poptd .* 200
    po,u0o = setpu(ps, poptd2, p, u0)
    @test eltype(po) == eltype(poptd2)
    @test typeof(po) == SVector{3, Dual128}
    @test convert(eltype(p), po[3]) == p[3]
    @test po[1:2] == poptd2[1:2]
    @test eltype(u0o) == eltype(poptd2)
    @test u0o[[1,3]] == u0[[1,3]]
    @test convert(eltype(u0),u0o[2]) == poptd2[3]
end;

@testset "getting parameters by name" begin
    #ps = ParSetter( SA[ks,km,kd], SA[s,b,r], SA[ks,km], SA[b])
    p = SA[1.1,2.2,3.3]
    @test p[parindex(ps, :km)] == 2.2
    @test isnothing(parindex(ps, :nonexisting))
end;

@testset "getting states name" begin
    #ps = ParSetter( SA[ks,km,kd], SA[s,b,r], SA[ks,km], SA[b])
    u0 = SA[11.1,22.2,33.3]
    @test u0[stateindex(ps, :r)] == 33.3
    @test isnothing(stateindex(ps, :nonexisting))
end;
