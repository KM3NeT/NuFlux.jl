![](https://github.com/KM3NeT/NuFlux.jl/raw/main/docs/src/assets/nuflux.svg)

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://km3net.github.io/NuFlux.jl/dev)
[![Build Status](https://github.com/KM3NeT/NuFlux.jl/workflows/CI/badge.svg)](https://github.com/KM3NeT/NuFlux.jl/actions)

# NuFlux.jl

**NuFlux.jl** is a package for providing neutrino flux related functionalities 
with focus on parsing the flux tables used in ["Atmospheric neutrino flux calculation using the NRLMSISE-00 atmospheric model" (10.1103/PhysRevD.92.023004)](https://link.aps.org/doi/10.1103/PhysRevD.92.023004) 
which are provided under [http://www.icrr.u-tokyo.ac.jp/~mhonda/](http://www.icrr.u-tokyo.ac.jp/~mhonda/).

## Installation

To install `NuFlux.jl`, clone this repo. Then open the Julia REPL as
```shell
julia --project=NuFlux.jl
```
Next step is to run the following command in the Julia REPL:

```julia
using Pkg
Pkg.instantiate()
```

Then you are good to go to use `NuFlux.jl`

## Example

Here’s a simple example to get started with `NuFlux.jl`. This example demonstrates how to load a flux table and calculate the flux for a given energy, zenith angle, and azimuth angle.

```julia
using NuFlux

# Load a flux table from a file, here we will use one of the default fluxes in the package, but you can use your own
NUFLUX_PATH = split(Base.pathof(NuFlux), "src")[1]
FLUX_DATA_DIR = joinpath(NUFLUX_PATH, "data")
flux_dict = NuFlux.readfluxfile(joinpath(FLUX_DATA_DIR, "frj-ally-20-12-solmin.d"))

# Select a specific flux table (e.g., for muon neutrinos)
muon_neutrino_flux = flux_dict[NuFlux.NUMU_PDGID]  # We can use NuFlux defined variables which are just instances of particles from Corpuscles.jl

# Calculate the flux for a given energy, zenith angle, and azimuth angle
energy = 10.0  # Energy in GeV
cosθ = 0.5     # Cosine of the zenith angle
ϕ = 0.0        # Azimuth angle in radians

flux_value = NuFlux.flux(muon_neutrino_flux, energy, cosθ, ϕ; interpol=true)
println("Flux value: ", flux_value)
```