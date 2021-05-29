using DelimitedFiles

# Load the precomputed Coulomb potential due to a pz orbital
coulomb_precompute = readdlm(
    "src/orbital_integration/coulomb_potential_pz_precompute_Rmax_30.dat",
)

# Obtain the number of points and the maximum value of R (from the file name)
nPts = size(coulomb_precompute)[1]
Rmax = 30;
# Get the appropriate arrays
Rs = range(0, Rmax, length = nPts)
τs = range(0, π / 2, length = nPts)

# Perform the interpolation
interp_call = LinearInterpolation((Rs, τs), coulomb_precompute)

# Create an interpolated Coulomb potential function. If the radial distance is
# greater than Rmax, approximate as 1 / R. Because of the pz shape, the potential
# is symmetric across the xy plane. Hence, if the polar angle τ > π / 2,
# we use π - τ (since the precalculation was performed for 0 < τ < π / 2)
function coulomb_potential_pz_interp(R, τ)
    if τ > π / 2
        τ = π - τ
    end
    if R > Rmax
        return 1 / R
    else
        return interp_call(R, τ)
    end
end

# Use the interpolated coulomb potential to calculate the interaction energy
# between two pz orbitals.
# NOTE: r is a spherical vector with components (r, θ, ϕ)
@inline function coulomb_energy_pz_pz_Integrand(
    r::Vector{Float64},
    ρ::Vector{Float64},
)
    # ρ is a CARTESIAN vector between the two pz orbitals
    # r is a SPHERICAL vector pointing from the center of one of the orbitals to a
    # point in space
    # d is a CARTESIAN vector pointing from the center of the second pz to the same
    # point
    dist =
        ρ - [
            r[1] * sin(r[2]) * cos(r[3]),
            r[1] * sin(r[2]) * sin(r[3]),
            r[1] * cos(r[2]),
        ]
    res =
        Ψ_pz(r[1], r[2]) .^ 2 *
        r[1]^2 *
        sin(r[2]) *
        coulomb_potential_pz_interp(norm(dist), acos((dist[3]) / norm(dist)))[1]
end

function coulomb_energy_pz_pz(ρ::Vector{Float64})
    hcubature(
        r -> coulomb_energy_pz_pz_Integrand(r, ρ),
        [0, 0, 0],
        [40, π, 2 * π],
        reltol = 1e-5,
    )
end

function coulomb_energy_proton_pz(ρ::Vector{Float64})
    R = norm(ρ)
    if R > 0
        τ = acos(ρ[3] / R)
    else
        τ = 0
    end
    res = coulomb_potential_pz(R, τ)
    return res
end

# Coordinates of the host atom and its neighbors. NOTE: Only one neighbor at
# each distance is given
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

proton_energies = @showprogress map(
    x -> coulomb_energy_proton_pz(x)[1],
    NN_Coords_Cartesian_3D,
)
ee_energies =
    @showprogress map(x -> coulomb_energy_pz_pz(x)[1], NN_Coords_Cartesian_3D)
