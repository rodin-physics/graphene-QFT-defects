using ProgressMeter
include("../src/computed_quantities/density.jl")

sys = new_graphene_system()
sys = set_μ(sys, 0.4)

# N_dopant = PerturbedAtom(-2.0, 0.0, graphene_A(0,0))
# sys = add_mod_atom(sys, N_dopant)

imp1 = new_impurity(0.5)
imp1 = add_coupling(imp1, -7.0, graphene_A(0,0))

sys = add_imp(sys, imp1)

nPts = 65;     # Number of grid points away from zero

# Arrays
d1s = -nPts:1:nPts; # d1 vectors
d2s = -nPts:1:nPts; # d2 vectors

# Lattice vector matrices
D1S = repeat(d1s, 1, 2 * nPts + 1);
D2S = repeat(d2s', 2 * nPts + 1, 1);
# scattering_atoms(Vector{ImpurityState}([]), [N_dopant])
# scattering_atoms(
# [add_coupling(new_impurity(0.5), -7.0, graphene_A(0,0))], [N_dopant])
#
