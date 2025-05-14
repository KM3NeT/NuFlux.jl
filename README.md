![](https://git.km3net.de/simulation/NuFlux.jl/-/raw/main/docs/src/assets/nuflux.svg)


[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://simulation.pages.km3net.de/NuFlux.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://simulation.pages.km3net.de/NuFlux.jl/dev)
[![pipeline status](https://git.km3net.de/simulation/NuFlux.jl/badges/main/pipeline.svg)](https://git.km3net.de/simulation/NuFlux.jl/-/commits/main)
[![coverage report](https://git.km3net.de/simulation/NuFlux.jl/badges/main/coverage.svg)](https://git.km3net.de/simulation/NuFlux.jl/-/commits/main)


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


## Usage
The function `NuFlux.readfluxfile(my/flux/file)` will read the Honda flux file of any site and return a dictionary where each key is the the table for a different neutrino type.

The main function which extracts the flux value is `NuFlux.flux(flux_of_specific_neutrino_type, energy, cosθ, ϕ; interpol=true)`

This function extracts the flux at a given energy, zenith direction cosθ, and azimuth direction ϕ. If the argument of the azimuth direction is missing it will take the mean of the flux values, if the value of the zenith direction is missing another mean will be computed.
The function allows to compute an interpolation from the table using the boolean of `interpol`.

The `flux` function allows to define which kind of interpolation to be used with the argument `interpol_method` which can be `"linear", "quadratic", "cubic"`, the default method is `"cubic"`. 
Additionally, the argument `interp_logflux` is a boolean where if `true` it will interpolate using the $\log_10$ values of the flux to allow for a smoother interpolation, the default value is `false`.

It is recommended to use the combination of `interpol_method="linear",interp_logflux=true` as this is shown to better approximate the actual behavior of the Honda fluxes.


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
