using Test
using MicMods
using ModelingToolkit

using StaticArrays
#ps = LabeledParSetter( SA[:ks,:km,:kd], SA[:s,:b,:r], SA[:ks,:km], SA[:b])
ps = LabeledParSetter( (:ks,:km,:kd), (:s,:b,:r), (:ks,:km), (:b,))

function tmpf()
    systemt = chak21_fixedr_system();
    ps = LabeledParSetter( systemt.syss, (:ks,:km), (:b,))
    #ps = LabeledParSetter( systemt.syss, SA[:ks,:km], SA[:b])
    ps = LabeledParSetter(systemt)
    # not type stable because paropt is not know at compile time
    # but using setpu and label_parsys when ps is passed is then type-stable
    # @code_warntype LabeledParSetter(systemt) 
    # @code_warntype LabeledParSetter( systemt.syss, (:ks,:km), (:b,))
end

using DualNumbers
using LabelledArrays

@testset "label SA" begin
    p = SA[1.1,2.2,3.3]
    u0 = SA[11.1,22.2,33.3]
    pl = @inferred label_parsys(ps, p)
    @test keys(pl) == getparsys(ps)
    u0l = @inferred label_statesys(ps, u0)
    @test keys(u0l) == getstatesys(ps)
    #@code_warntype label_parsys(ps, p)
end;

@testset "label MVector" begin
    p = @MVector [1.1,2.2,3.3]
    u0 = @MVector [11.1,22.2,33.3]
    pl = @inferred label_parsys(ps, p)
    #@code_warntype label_parsys(ps, p)
    @test pl isa LArray
    @test keys(pl) == getparsys(ps)
    u0l = @inferred label_statesys(ps, u0)
    @test u0l isa LArray
    @test keys(u0l) == getstatesys(ps)
    #@code_warntype label_parsys(ps, p)
    pl.ks = 4.4
    @test p[1] == 4.4 # changed the underlying vector
end;


@testset "setting parameters SA" begin
    p = SA[1.1,2.2,3.3]
    u0 = SA[11.1,22.2,33.3]
    popt = @inferred getpopt(ps,p,u0)
    @test popt == SLVector(ks=p[1],km=p[2],b=u0[2])
    popt = @inferred getpopt(ps,p,u0, Val(false))
    @test popt == SVector{3,Num}(p[1:2]...,u0[2])
    #
    po,u0o = @inferred setpu(ps, popt, p, u0, Val(false))
    po,u0o = @inferred setpu(ps, popt, p, u0)
    @test po isa LArray # but underlying storage is still Static
    @test_throws ErrorException po[3] = 4.4
    @test keys(po) == getparsys(ps)
    @test po[3] == p[3]
    @test po[1:2] == popt[1:2]
    @test u0o isa LArray
    @test keys(u0o) == getstatesys(ps)
    @test u0o[[1,3]] == u0[[1,3]]
    @test u0o[2] == popt[3]
    #
    popt2 = popt .* 200
    po,u0o = @inferred setpu(ps, popt2, p, u0)
    @test po isa LArray
    @test po[3] == p[3]
    @test po[1:2] == popt2[1:2]
    @test u0o isa LArray
    @test u0o[[1,3]] == u0[[1,3]]
    @test u0o[2] == popt2[3]
end;

@testset "setting parameters MVector" begin
    p = @MVector [1.1,2.2,3.3]
    u0 = @MVector [11.1,22.2,33.3]
    popt = @inferred getpopt(ps,p,u0)
    @test popt == LVector(ks=p[1],km=p[2],b=u0[2])
    popt = @inferred getpopt(ps,p,u0, Val(false))
    @test popt == MVector{3,Num}(p[1:2]...,u0[2])
    #
    po,u0o = @inferred setpu(ps, popt, p, u0, Val(false))
    #@code_warntype setpu(ps, popt, p, u0, Val(false))
    po,u0o = @inferred setpu(ps, popt, p, u0)
    #@code_warntype setpu(ps, popt, p, u0)
    @test po isa LArray
    @test keys(po) == getparsys(ps)
    @test po[3] == p[3]
    @test po[1:2] == popt[1:2]
    po.ks = 4.4
    @test p[1] == 1.1 # setpu created a copy -> not changed
    @test u0o isa LArray
    @test keys(u0o) == getstatesys(ps)
    @test u0o[[1,3]] == u0[[1,3]]
    @test u0o[2] == popt[3]
    #
    popt2 = popt .* 200
    po,u0o = @inferred setpu(ps, popt2, p, u0)
    @test po isa LArray
    @test po[3] == p[3]
    @test po[1:2] == popt2[1:2]
    @test u0o isa LArray
    @test u0o[[1,3]] == u0[[1,3]]
    @test u0o[2] == popt2[3]
end;


# with Vector labelling is not supported, because not type stable
@testset "setting parameters Vector" begin
    p = [1.1,2.2,3.3]
    u0 = [11.1,22.2,33.3]
    popt = @inferred getpopt(ps,p,u0)
    #
    #@code_warntype label_parsys(ps,p) 
    po,u0o = @inferred setpu(ps, popt, p, u0, Val(false)) 
    @test po isa typeof(p)
    @test u0o isa typeof(u0)
    po,u0o = @inferred setpu(ps, popt, p, u0) 
    @test po isa LArray
    @test po[3] == p[3]
    @test po[1:2] == popt[1:2]
    @test u0o isa LArray
    @test u0o[[1,3]] == u0[[1,3]]
    @test u0o[2] == popt[3]
    #
    popt2 = popt .* 200
    po,u0o = setpu(ps, popt2, p, u0)
    @test po isa LArray
    @test po[3] == p[3]
    @test po[1:2] == popt2[1:2]
    @test u0o isa LArray
    @test u0o[[1,3]] == u0[[1,3]]
    @test u0o[2] == popt2[3]
end;

@testset "setting parameters DualNumbers" begin
    p = @MVector [1.1,2.2,3.3]
    u0 = @MVector [11.1,22.2,33.3]
    popt = @inferred getpopt(ps,p,u0)
    poptd = map(Dual, popt)
    typeof(poptd)
    eltype(poptd)
    #pd = map(eltype(poptd), p)
    #u0d = map(eltype(poptd), u0)
    #
    po,u0o = @inferred setpu(ps, poptd, p, u0) # not type stable
    @test po isa LArray
    @test u0o isa LArray
    #po,u0o = @inferred setpu(ps, poptd, pd, u0d) # not type stable
    #@code_warntype setpu(ps, poptd, p, u0, Val(false))
    @test eltype(po) == eltype(poptd)
    @test convert(eltype(p), po[3]) == p[3]
    @test po[1:2] == poptd[1:2]
    @test eltype(u0o) == eltype(poptd)
    @test u0o[[1,3]] == u0[[1,3]]
    @test convert(eltype(u0),u0o[2]) == poptd[3]
    #
    poptd2 = poptd .* 200
    po,u0o = @inferred setpu(ps, poptd2, p, u0)
    #po,u0o = @inferred setpu(ps, poptd2, pd, u0d)
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
    # pd = map(eltype(poptd), p)
    # u0d = map(eltype(poptd), u0)
    #
    #po,u0o = setpu(ps, poptd, pd, u0d)
    po,u0o = setpu(ps, poptd, p, u0)
    @test po isa LArray
    @test u0o isa LArray
    @test eltype(po) == eltype(poptd)
    @test convert(eltype(p), po[3]) == p[3]
    @test po[1:2] == poptd[1:2]
    @test eltype(u0o) == eltype(poptd)
    @test u0o[[1,3]] == u0[[1,3]]
    @test convert(eltype(u0),u0o[2]) == poptd[3]
    #
    poptd2 = poptd .* 200
    #po,u0o = setpu(ps, poptd2, pd, u0d)
    po,u0o = setpu(ps, poptd2, p, u0)
    @test eltype(po) == eltype(poptd2)
    @test po isa LArray
    @test convert(eltype(p), po[3]) == p[3]
    @test po[1:2] == poptd2[1:2]
    @test u0o isa LArray
    @test u0o[[1,3]] == u0[[1,3]]
    @test convert(eltype(u0),u0o[2]) == poptd2[3]
end;

@testset "LabelledArrays" begin
    p = LVector(ks=1.1, km=2.2, kd=3.3)
    u0 = LVector(s=11.1, b=22.2, r= 33.3)
    label_parsys(ps, p)
    popt = getpopt(ps,p,u0)
    @test popt == SA[p[1:2]...,u0[2]]
    typeof(popt)
    #
    po,u0o = @inferred setpu(ps, popt, p, u0, Val(false))
    #@code_warntype setpu(ps, popt, p, u0, Val(false))
    #po,u0o = @inferred setpu(ps, popt, p, u0)
    #@code_warntype setpu(ps, popt, p, u0, Val(true))
    @test typeof(po) == typeof(p)
    @test typeof(u0o) == typeof(u0)
    @test keys(po) == keys(p)
    @test keys(u0o) == keys(u0)
    @test po[3] == p[3]
    @test po[1:2] == popt[1:2]
    @test u0o[[1,3]] == u0[[1,3]]
    @test u0o[2] == popt[3]
    #
    popt2 = popt .* 200
    po,u0o = setpu(ps, popt2, p, u0, Val(false))
    @test typeof(po) == typeof(p)
    @test po[3] == p[3]
    @test po[1:2] == popt2[1:2]
    @test typeof(u0o) == typeof(u0)
    @test u0o[[1,3]] == u0[[1,3]]
    @test u0o[2] == popt2[3]
    #
    po,u0o = setpu(ps, popt2, p, u0) # names stored twice, here
    @test po isa LArray
    @test keys(po) == getparsys(ps)
    @test po[3] == p[3]
    @test po[1:2] == popt2[1:2]
end;

