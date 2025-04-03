module NuFlux

using DocStringExtensions
using DelimitedFiles
using StatsBase
using Statistics
using Corpuscles
using Interpolations
using LRUCache

const NUE_PDGID = Particle("nu(e)0").pdgid.value
const ANUE_PDGID = Particle("~nu(e)0").pdgid.value
const NUMU_PDGID = Particle("nu(mu)0").pdgid.value
const ANUMU_PDGID = Particle("~nu(mu)0").pdgid.value

"""
Abstract type representing a neutrino flux model.
"""
abstract type Flux end

"""
    FluxTable

A struct representing a tabulated neutrino flux model.

# Fields
- `coszentihbinedges::Vector{Float64}`: Bin edges for the cosine of the zenith angle.
- `azimuthbinedges::Vector{Float64}`: Bin edges for the azimuth angle.
- `energies::Vector{Float64}`: Energy values in GeV.
- `particle::Particle`: The particle type (e.g., neutrino or antineutrino).
- `flux::Array{Float64, 3}`: A 3D array of flux values, indexed by (cosθ, ϕ, energy).
"""
struct FluxTable <: Flux
    coszentihbinedges::Vector{Float64}
    azimuthbinedges::Vector{Float64}
    energies::Vector{Float64}
    particle::Particle
    flux::Array{Float64, 3}
end

"""
    FluxBinInfo

A struct representing bin information for a flux table.

# Fields
- `CosZenithBin::Tuple{Rational, Rational}`: Tuple of (min, max) cosine of the zenith angle for the bin.
- `AzimutBin::Tuple{Rational, Rational}`: Tuple of (min, max) azimuth angle for the bin.
"""
struct FluxBinInfo
    CosZenithBin::Tuple{Rational, Rational}
    AzimutBin::Tuple{Rational, Rational}
end

"""
    FluxBinInfo(cosZmin, cosZmax, azimutmin, azimutmax)

Construct a `FluxBinInfo` instance from the given bin edges.

# Arguments
- `cosZmin`: Minimum cosine of the zenith angle.
- `cosZmax`: Maximum cosine of the zenith angle.
- `azimutmin`: Minimum azimuth angle.
- `azimutmax`: Maximum azimuth angle.
"""
FluxBinInfo(cosZmin, cosZmax, azimutmin, azimutmax) = FluxBinInfo((cosZmin, cosZmax), (azimutmin, azimutmax)) 

"""
    parsefluxbininfo(line)

Parse a line of text into a `FluxBinInfo` instance.

# Arguments
- `line`: A string containing bin edge information.

# Returns
- `FluxBinInfo`: A struct containing the parsed bin edges.
"""
function parsefluxbininfo(line)
    regex = eachmatch(r"-?\d+\.?\d*", line)
    numbers = collect(rationalize(parse(Float64, m.match)) for m in regex)
    FluxBinInfo(numbers...)
end

"""
    _getbinmids(binedges)

Calculate the midpoints of bins given their edges.

# Arguments
- `binedges`: A vector of bin edges.

# Returns
- `Vector{Float64}`: A vector of bin midpoints.
"""
function _getbinmids(binedges)
    Float64.(mean.(binedges[2:end] .+ binedges[1:end-1]) ./ 2)
end

"""
    _makerange(x)

Create a range from a vector of values, ensuring constant spacing.

# Arguments
- `x`: A vector of values.

# Returns
- `StepRangeLen`: A range with constant spacing.

# Throws
- `ErrorException`: If the spacing between values is not constant.
"""
function _makerange(x)
    diffs = x[2:end] - x[1:end-1]
    std(diffs) > 1e-4 && error("value spacing not constant")
    range(minimum(x), maximum(x), length=length(x))
end

"""
    readfluxfile(io)

Read a flux file and construct a vector of `FluxTable` instances.

# Arguments
- `io`: An IO stream or file path to the flux data.

# Returns
- `Dict{FluxTable}`: A dictionary of `FluxTable` instances representing the flux data indexed by PDG ID.
"""
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
    for (i,p) in enumerate([Particle(NUMU_PDGID), Particle(ANUMU_PDGID), Particle(NUE_PDGID), Particle(ANUE_PDGID)])
        tmp = data[:,:,:,i]
        push!(retval, FluxTable(coszedges, azimuthedges, energies, p, tmp))
    end
    retdict = Dict(NUE_PDGID => retval[3],
        NUMU_PDGID => retval[1],
        ANUE_PDGID => retval[4],
        ANUMU_PDGID => retval[2],)
    retdict
end

"""
    readfluxfile(filepath::String)

Read a flux file from a file path and construct a vector of `FluxTable` instances.

# Arguments
- `filepath`: Path to the flux data file.

# Returns
- `Vector{FluxTable}`: A vector of `FluxTable` instances representing the flux data.
"""
function readfluxfile(filepath::String)
    io = open(filepath, "r")
    ret = readfluxfile(io)
    close(io)
    ret
end

"""
    _getbinindex(binedges, x)

Find the index of the bin that contains the value `x`.

# Arguments
- `binedges`: A vector of bin edges.
- `x`: The value to find the bin for.

# Returns
- `Int`: The index of the bin containing `x`.
"""
function _getbinindex(binedges, x)
    if x ≈ first(binedges)
        return 1
    else
        return searchsortedfirst(binedges, x)-1
    end
end

const lru_energy_etp = LRU{NuFlux.FluxTable, Interpolations.Extrapolation}(maxsize=20000)

"""
$(SIGNATURES)

# Arguments
- `flux`:       Flux data 
- `energy`:     Energy in GeV
- `interpol`:   Interpolate the data

# Returns
- `Float64`: The flux value at the given energy.
"""
function flux(f::FluxTable, energy::S; interpol::Bool=false) where {S <: Real}
    if interpol
        etp = get!(lru_energy_etp, f) do
            tmp = mean(f.flux, dims=(1,2))[1,1,:]
            itp = interpolate(tmp, BSpline(Cubic(Line(OnGrid()))))
            sitp = scale(itp, _makerange(log10.(f.energies)))
            extrapolate(sitp, Flat())
        end
        return etp(log10(energy))
    else
        idx_energy = findmin(abs.(f.energies .- energy))[2]
        return mean(f.flux[:, :, idx_energy])
    end
end

const lru_zenith_energy_etp = LRU{NuFlux.FluxTable, Interpolations.Extrapolation}(maxsize=20000)

"""
$(SIGNATURES)

# Arguments
- `flux`:       Flux data 
- `energy`:     Energy in GeV
- `cosθ`:       Cosine of the zenith angle
- `interpol`:   Interpolate the data

# Returns
- `Float64`: The flux value at the given energy and zenith angle.
"""
function flux(f::FluxTable, energy::S, cosθ::T; interpol::Bool=false) where {S,T <: Real}
    if interpol
        etp = get!(lru_zenith_energy_etp, f) do
            tmp = mean(f.flux, dims=2)[:,1,:]
            itp = interpolate(tmp, BSpline(Cubic(Line(OnGrid()))))
            sitp = scale(itp, _makerange(_getbinmids(f.coszentihbinedges)), _makerange(log10.(f.energies)))
            extrapolate(sitp, Flat())
        end
        return etp(cosθ, log10(energy))
    else
        idx_energy = findmin(abs.(f.energies .- energy))[2]
        idx_coszenith = _getbinindex(f.coszentihbinedges, cosθ)
        return mean(f.flux[idx_coszenith, :, idx_energy])
    end
end

const lru_zenith_azimuth_energy_etp = LRU{NuFlux.FluxTable, Interpolations.Extrapolation}(maxsize=20000)

"""
$(SIGNATURES)

# Arguments
- `flux`:       Flux data 
- `energy`:     Energy in GeV
- `cosθ`:       Cosine of the zenith angle
- `ϕ`:          Azimuth angle
- `interpol`:   Interpolate the data

# Returns
- `Float64`: The flux value at the given energy, zenith angle, and azimuth angle.
"""
function flux(f::FluxTable, energy::S, cosθ::T, ϕ::U; interpol::Bool=false) where {S,T,U <: Real}
    if interpol
        etp = get!(lru_zenith_azimuth_energy_etp, f) do
            itp = interpolate(f.flux, BSpline(Cubic(Line(OnGrid()))))
            sitp = scale(itp, _makerange(_getbinmids(f.coszentihbinedges)), _makerange(_getbinmids(f.azimuthbinedges)), _makerange(log10.(f.energies)))
            extrapolate(sitp, Flat())
        end
        return etp(cosθ, ϕ, log10(energy))
    else
        idx_energy = findmin(abs.(f.energies .- energy))[2]
        idx_azimuth = _getbinindex(f.azimuthbinedges, ϕ)
        idx_coszenith = _getbinindex(f.coszentihbinedges, cosθ)
        return f.flux[idx_coszenith, idx_azimuth, idx_energy]
    end
end

end # module
