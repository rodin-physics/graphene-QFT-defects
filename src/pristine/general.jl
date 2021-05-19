using DelimitedFiles
using LinearAlgebra
using ProgressMeter
using QuadGK
using Statistics

## Parameters
const ν = 1e-2;     # Small number for relative tolerance
const α = 1e-10;     # Small number for absolute tolerance
const η = 1e-3;     # Small number for moving the contour off the real axis
const nevals = 1e8; # Maximum number of evaluations in integrals

# Graphene hopping integral in eV, overlap term, and lattice vectors in Angstroms
const t0 = 2.74;
const P = 0.065;
const d = 2.46;
const d1 = d .* [+1, √(3)] ./ 2;
const d2 = d .* [-1, √(3)] ./ 2;

const Γ = [0, 0];
const K = [2 * pi / 3 / d, 2 * pi / sqrt(3) / d];
const M = [0, 2 * pi / sqrt(3) / d];

# Pristine Hamiltonian
function Hπ(q)
    f1 = 1 + exp(1im * sum(d1 .* q)) + exp(1im * sum(d2 .* q))
    return ([
        0 (-t0*f1)
        conj(-t0 * f1) 0
    ])
end

# Overlap matrix
function Overlap(q)
    # Phase terms
    f1 = 1 + exp(1im * sum(d1 .* q)) + exp(1im * sum(d2 .* q))

    return ([
        1 (P*f1)
        conj(P * f1) 1
    ])
end

# Fermi-Dirac Distribution
function nF(x, T)
    if T == 0
        return (convert(Float64, x < 0))
    else
        return 1 / (1 + exp(x / T))
    end
end
