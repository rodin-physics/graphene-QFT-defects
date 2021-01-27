using Cuba
using DelimitedFiles
using LinearAlgebra
using Plots
using QuadGK
using Statistics

## Parameters
const ν = 1e-2;     # Small number for relative tolerance
const α = 1e-8;     # Small number for absolute tolerance
const η = 1e-2;     # Small number for moving the contour off the real axis
const nevals = 1e6; # Maximum number of evaluations in integrals

# Graphene hopping integral in eV and lattice vectors in Angstroms
const t = 2.8;
const d = 2.46;
const d1 = d .* [+1, √(3)] ./ 2;
const d2 = d .* [-1, √(3)] ./ 2;

# Fermi-Dirac Distribution
function nF(x, T)
    if T == 0
        return (convert(Float64, x < 0))
    else
        return 1 / (1 + exp(x / T))
    end
end
