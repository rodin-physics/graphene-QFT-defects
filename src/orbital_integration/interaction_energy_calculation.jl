include("orbital_data_import.jl")

# Coordinates of the host atom and its neighbors
host_atom = graphene_A(0, 0)
NN_1 = graphene_B(0, 0)
NN_2 = graphene_A(1, 0)
NN_3 = graphene_B(1, 1)
NN_4 = graphene_B(2, 0)
NN_5 = graphene_A(1, 1)
NN_6 = graphene_A(2, 0)
NN_7 = graphene_B(2, 1)
NN_8 = graphene_B(3, -1)
NN_9 = graphene_B(3, 0)
NN_10 = graphene_A(3, -1)
NN_11 = graphene_B(2, 2)
NN_12 = graphene_A(3, 0)
NN_13 = graphene_B(3, 1)
NN_14 = graphene_B(-2, -1)
NN_15 = graphene_A(2, 2)
NN_16 = graphene_B(4, 0)
NN_17 = graphene_A(3, 1)

# Array of coordinates
NN_Coords = [
    host_atom,
    NN_1,
    NN_2,
    NN_3,
    NN_4,
    NN_5,
    NN_6,
    NN_7,
    NN_8,
    NN_9,
    NN_10,
    NN_11,
    NN_12,
    NN_13,
    NN_14,
    NN_15,
    NN_16,
    NN_17,
]

# Convert the lattice coordinates of each neighbor to 3D Cartesian coordinates
# NOTE: We convert the distances into Bohr radii
NN_Coords_Cartesian_3D = map(NN_Coords) do c
    loc = crystal_to_cartesian(c)
    return [loc.x / a0, loc.y / a0, 0.0]
end

proton_energies =
    @showprogress map(coulomb_energy_proton, NN_Coords_Cartesian_3D)
ee_energies =
    @showprogress map(coulomb_energy_pz, NN_Coords_Cartesian_3D[2:end])

proton_energies_approx = 1 ./ norm.(NN_Coords_Cartesian_3D)
proton_energies_approx - map(x -> x[1], proton_energies)
