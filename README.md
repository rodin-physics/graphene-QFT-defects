# Defects in Graphene: QFT approach

## Theory

The derivation of the formalism used in this package can be found [here](https://github.com/rodinalex/graphene-QFT-defects/blob/main/Theory/Theory.pdf).

For a more detailed discussion, see https://doi.org/10.1103/PhysRevB.102.195416

## Use

### Basics

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

### Scripts

To facilitate the use of the package, a few useful files can be found in the `calculation` folder.

#### Density
To take advantage of parallel computation when computing the induced charge density over many sites
1. Start by defining the system in the `computations_settings.jl` file. This includes creating a `new_graphene_system()` and populating it with defects, as outlined above
2. In the same file, define the size of the system by specifying how many unit vectors it spans by setting `nPts` to a desired value. This means that the lattice will go from -`nPts` to +`nPts` for each of the unit vectors so that the total number of points on the lattice to be computed will then be 2 x (2`nPts`+1) x (2`nPts`+1). The "2" in the front originates from the two sublattices.
3. In `rho_calculation.jl`, define the number of workers to be used in the parallel calculation by setting `procs` to a desired value.
4. Make sure that the targed directory to save the output exists and run the script! The induced density will be calculated at every lattice point and saved to two files, one for each sublattice.

#### Spectral function

### Analysis and Plotting

#### Plotting

When one performs the lattice-wide calculations, as described above, the output is saved to two files. These files only contain the values at the lattice sites without the atom coordinates. To plot these results, further processing is necessary. A function called `Data_Process(lattice_A_data, lattice_B_data)` in `src/analysis/dataprocess.jl` takes the data for the two sublattices and combines it by matching the data to the appropriate spatial coordinate, making it possible to plot the results. For more detail, see `calculation/rho_Plotter.jl`.

#### Fourier Transform



## ToDo

1. Implement the energy calculation
2. Fix up the Fourier Transform part
