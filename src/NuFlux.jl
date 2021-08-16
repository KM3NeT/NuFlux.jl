module NuFlux

using DelimitedFiles
using StatsBase
using Statistics

@enum NeutrinoType begin
    Muon = 1
    AntiMuon =2
    Electron = 3
    AntiElectron = 4
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

function _generateenergybinedges(energies::Vector{T}) where T <: Real
    vcat(energies, Inf)
end

function readfluxfile(io)
    lines = readlines(io)
    idx = findall(l->occursin("average",l), lines)
    fluxes = Dict()
    coszedges = Set(Vector{Rational}())
    azimuthedges = Set(Vector{Rational}())
    energyedges = nothing
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
        if isnothing(energyedges)
            energyedges = _generateenergybinedges(tmp[:,1])
        end
    end
    @show energyedges
    coszedges = sort(collect(coszedges))
    azimuthedges = sort(collect(azimuthedges))
    data = zeros(Float64, length(coszedges)-1, length(azimuthedges)-1, length(energyedges)-1, 4)
    for (k,v) in fluxes
        coszbin = searchsortedfirst(coszedges, mean(k.CosZenithBin))-1
        azimutbin = searchsortedfirst(azimuthedges, mean(k.AzimutBin))-1
        data[coszbin, azimutbin, :, :] = v
    end
    retval = Dict{NeutrinoType, Histogram}()
    for nutype in 1:NeutrinoType.size
        retval[NeutrinoType(nutype)] = Histogram((coszedges, azimuthedges, energyedges), data[:,:,:,nutype])
    end
    retval
end

function readfluxfile(filepath::String)
    io = open(filepath, "r")
    ret = readfluxfile(io)
    close(io)
    ret
end

function (h::Histogram)(x...)
    dims = length(h.edges)
    (length(x) != dims) && throw(DimensionMismatch("Size of arguments not equal to the histogram dimensions"))
    idx = zeros(Int64, dims)
    for i in 1:dims
        idx[i] = searchsortedfirst(h.edges[i], x[i])
    end
    h.weights[idx...]
end


end # module
