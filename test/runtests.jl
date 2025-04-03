using Test
using NuFlux

const DATA_DIR = joinpath(@__DIR__, "../data")

@testset "load fluxtable" begin
    @test !isa(try NuFlux.readfluxfile(joinpath(DATA_DIR, "frj-ally-20-12-solmin.d")) catch e e end, Exception)
    flux_table = NuFlux.readfluxfile(joinpath(DATA_DIR, "frj-ally-20-12-solmin.d"))
    @test !isa(try NuFlux.flux(flux_table[14], 10, 0; interpol=true) catch e e end, Exception)
end

