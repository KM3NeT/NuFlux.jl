using Test
using NuFlux

const DATA_DIR = joinpath(@__DIR__, "../data")

@testset "load fluxtable" begin
    @test !isa(try NuFlux.readfluxfile(joinpath(DATA_DIR, "frj-ally-20-12-solmin.d")) catch e e end, Exception)
end

@testset "power law flux" begin
    f = NuFlux.PowerLawFlux(-2)
    @test 1e-12 ≈ NuFlux.flux(f, 1e2)
end

@testset "power law flux w. cutoff" begin
    f = NuFlux.PowerLawFlux(-2, 1, 3e6)
    @test 1e-12 ≈ NuFlux.flux(f, 1e2) atol=1e-15
    @test 7.165e-21 ≈ NuFlux.flux(f, 1e6) atol=1e-24
end
