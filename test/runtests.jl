using Test
using NuFlux

const DATA_DIR = joinpath(@__DIR__, "../data")

@testset "load fluxtable" begin
    @test !isa(try NuFlux.readfluxfile(joinpath(DATA_DIR, "frj-ally-20-12-solmin.d")) catch e e end, Exception)
end

