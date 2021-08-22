module NuFlux

using DocStringExtensions
using DelimitedFiles
using StatsBase
using Statistics
using Corpuscles
using Interpolations

abstract type Flux end

struct FluxTable <: Flux
    coszentihbinedges::Vector{Float64}
    azimuthbinedges::Vector{Float64}
    energies::Vector{Float64}
    particle::Particle
    flux::Array{Float64, 3}
    sitp::Interpolations.Extrapolation
end

struct FluxBinInfo
    CosZenithBin::Tuple{Rational, Rational}
    AzimutBin::Tuple{Rational, Rational}
end

FluxBinInfo(cosZmin, cosZmax, azimutmin, azimutmax) = FluxBinInfo((cosZmin, cosZmax), (azimutmin, azimutmax)) 

function parsefluxbininfo(line)
    regex = eachmatch(r"-?\d+\.?\d*", line)
    numbers = collect(rationalize(parse(Float64, m.match)) for m in regex)
    FluxBinInfo(numbers...)
end

function _getbinmids(binedges)
    Float64.(mean.(binedges[2:end] .+ binedges[1:end-1]) ./ 2)
end

function _makerange(x)
    diffs = x[2:end] - x[1:end-1]
    std(diffs) > 1e-4 && error("value spacing not constant")
    minimum(x):mean(diffs):maximum(x)
end

function readfluxfile(io)
    lines = readlines(io)
    idx = findall(l->occursin("average",l), lines)
    fluxes = Dict()
    coszedges = Set(Vector{Rational}())
    azimuthedges = Set(Vector{Rational}())
    energies = nothing
    for (i,j) in collect(zip(idx[1:end], [idx[2:end]..., lastindex(lines)+1]))
        bininfo = parsefluxbininfo(lines[i])
        coszedges = union(coszedges, Set(bininfo.CosZenithBin))
        azimuthedges = union(azimuthedges, Set(bininfo.AzimutBin))
        bfr = IOBuffer()
        for l in lines[i+2:j-1]
            write(bfr, strip(l))
            write(bfr, "\n")
        end
        seekstart(bfr)
        tmp = readdlm(bfr, ' ','\n')
        fluxes[bininfo] = tmp[:,2:end]
        if isnothing(energies)
            energies = tmp[:,1]
        end
    end
    coszedges = sort(collect(coszedges))
    azimuthedges = sort(collect(azimuthedges))
    data = zeros(Float64, length(coszedges)-1, length(azimuthedges)-1, length(energies), 4)
    for (k,v) in fluxes
        coszbin = searchsortedfirst(coszedges, mean(k.CosZenithBin))-1
        azimutbin = searchsortedfirst(azimuthedges, mean(k.AzimutBin))-1
        data[coszbin, azimutbin, :, :] = v
    end
    retval = Vector{FluxTable}()
    for (i,p) in enumerate([Particle(14), Particle(-14), Particle(12), Particle(-12)])
        tmp = data[:,:,:,i]
        itp = interpolate(tmp, BSpline(Cubic(Line(OnGrid()))))
        sitp = scale(itp, _makerange(_getbinmids(coszedges)), _makerange(_getbinmids(azimuthedges)), _makerange(log10.(energies)))
        etp = extrapolate(sitp, Flat())
        push!(retval, FluxTable(coszedges, azimuthedges, energies, p, tmp, etp))
    end
    retval
end

function readfluxfile(filepath::String)
    io = open(filepath, "r")
    ret = readfluxfile(io)
    close(io)
    ret
end

function _getbinindex(binedges, x)
    if x ≈ first(binedges)
        return 1
    else
        return searchsortedfirst(binedges, x)-1
    end
end

"""
$(SIGNATURES)

# Arguments
- `flux`:           Flux data 
- `energy`:         Energy in GeV
- `cosθ`:           Cosine of the zenith angle
- `ϕ`:              Azimuth angle
- `interpolate`:    Interpolate the data
"""
function flux(f::FluxTable, energy::S, cosθ::T, ϕ::U; interpolate::Bool=false) where {S,T,U <: Real}
    if interpolate
        return f.sitp(cosθ, ϕ, log10(energy))
    else
        idx_energy = findmin(abs.(f.energies .- energy))[2]
        idx_azimuth = _getbinindex(f.azimuthbinedges, ϕ)
        idx_coszenith = _getbinindex(f.coszentihbinedges, cosθ)
        return f.flux[idx_coszenith, idx_azimuth, idx_energy]
    end
end

end # module
