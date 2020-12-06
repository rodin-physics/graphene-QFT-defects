# Defects in Graphene: QFT approach

## Theory

The derivation of the formalism used in this package can be found [here](https://github.com/rodinalex/graphene-QFT-defects/blob/main/Theory/Theory.pdf).

For a more detailed discussion, see https://doi.org/10.1103/PhysRevB.102.195416

## Use

This package can calculate the electronic density at impurity states and graphene atoms, see [here](https://github.com/rodinalex/graphene-QFT-defects/blob/main/example/example.jl).

The general flow goes as follows:
1. Create a new system by calling `new_graphene_system()` and assigning it to a variable, say `my_system`. This produces a pristine graphene monolayer with μ = 0, T = 0, no defects, and no modified atoms.
2. If desired, one can create impurity states by calling `new_impurity(E::Float64)`, where E is the state energy, and assigning the output to a variable.
3. Once the impurities are created, it is possible to couple them to graphene sites using `add_coupling(imp::ImpurityState, V::Float64, coord::GrapheneCoord)`, where `V` is the coupling strength. For example, we want to connect `imp1` to a graphene site of the A sublattice at the unit cell (3,4). This is accomplished as follows: `imp1 = add_coupling(imp1, V, graphene_A(3,4))`. The procedure can be repeated to connect the impurity to more graphene sites. Couplings can be removed by calling `imp1 = remove_coupling(imp1, n)`, where `n` is the index of the coupling in the array.
4. Next, the impurities need to be included in `my_system` by calling `my_system = add_imp(my_system, imp1)`.
5. In a similar fashion, one can create perturbed atoms by calling `perturbed_atom = perturbed_atom(E::Float64, dt::Float64, coord::GrapheneCoord)`, where `E` is the new on-site energy, `dt` is the modification of the hopping parameter with the atom's neighbors, and coord is the atom's location. Just as with the impurities, the atoms must be added to the system by calling `my_system = add_mod_atom(my_system, perturbed_atom)`
6. It is possible to remove impurities or perturbed atoms by calling `my_system = remove_imp(my_system, n)` or `my_system = remove_mod_atom(my_system, n)`, where `n` is the index of the impurity/atom in the array.
7. The chemical potential and the temperature of the system can be changed by calling `my_system = set_μ(my_system, new_mu)` or `my_system = set_T(my_system, new_T)`.
8. Finally, one can calculate the charge densities. To get the defect-induced charge density in graphene, call `δρ_Graphene(coord, my_system)`, where coord is the coordinate of the atom. For the impurities, call `ρ_All(my_system)` if all the occupation numbers are desired, or `ρ_Impurity(n, my_system)` to only get the occupation number for the `n`'th impurity.
9. The spectral function at energy ω can be obtained by calling `spectral_graphene(ω, coord, my_system)` or `spectral_impurity(ω, n, my_system)`

## ToDo

1. Implement the energy calculation
2. Fix up the Fourier Transform part
3. Write up the plotting explanation
