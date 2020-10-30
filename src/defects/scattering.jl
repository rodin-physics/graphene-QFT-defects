include("propagator.jl")


# A function that takes a list of impurities and another one of directly
# modified atoms and returns a list of all the scattering atoms
function scattering_atoms(
    imps::Vector{ImpurityState},
    mod_atoms::Vector{PerturbedAtom},
)
    # Start by obtaining the atoms that are coupled to the impurities
    # Get all the couplings from the impurity states
    all_couplings = collect(Iterators.flatten(map(x -> x.coupling, imps)))
    # Turn the couplings into PerturbedAtoms
    hosts = map(x -> PerturbedAtom(0, x.δt, x.coord), all_couplings)
    # Combine all the atoms into one array
    all_perturbed_atoms = vcat(hosts, mod_atoms)
    # Some of the atoms could be coupled to multiple impurities and also be
    # modified by an external potential/substitution, so they will show up more
    # than once in all_perturbed_atoms, possibly with different values of t and ϵ
    # since different impurities could modify the hopping differently.
    # The way we resolve the conflicts is as follows:
    # For a given coordinate, if there is an entry with ϵ ≠ 0, it takes
    # precendence. If there are multiple different values of t ≠ 0, we take the
    # average. This is done using the conf_resolve() function defined below

    # Get the unique coordinates of the atoms
    u_perturbed_c = unique(map(x -> x.coord, all_perturbed_atoms))
    # Obtain the unique PerturbedAtom's through conf_resolve()
    u_perturbed_atoms = map(
        y -> conf_resolve(filter(x -> x.coord == y, all_perturbed_atoms)),
        u_perturbed_c,
    )

    # Finally, we add unique neighbors to the list of the perturbed atoms
    # First, obtain the coordinates of all the perturbed atoms' neighbors
    all_neighbors_c = collect(Iterators.flatten(map(neighbors, u_perturbed_c)))
    # Turn the neighbors into PerturbedAtoms with no energy or bond modification
    all_neighbors_atoms = map(x -> PerturbedAtom(0, 0, x), all_neighbors_c)
    # Combine all the unique perturbed atoms and their neighbors
    all_atoms = vcat(u_perturbed_atoms, all_neighbors_atoms)
    # It is possible that a given graphene atom will show up as a perturbed
    # atom and as a neighbor. We use the conf_resolve() function here to remove
    # duplicate coordinates.

    # Get the unique coordinates of the atoms
    u_c = unique(map(x -> x.coord, all_atoms))
    u_atoms = map(y -> conf_resolve(filter(x -> x.coord == y, all_atoms)), u_c)
    return u_atoms
end

# This function takes a list of PerturbedAtom's with the same coordinate and
# returns a single instance of a PerturbedAtom
function conf_resolve(atoms::Vector{PerturbedAtom})
    ϵs = filter(y -> y > 0, map(x -> x.ϵ, atoms))
    δts = filter(y -> y > 0, map(x -> x.δt, atoms))

    if length(ϵs) > 1
        error("Multiple definitions of on-site energy")
    elseif length(ϵs) == 0
        ϵ = 0
    else
        ϵ = ϵs[1]
    end

    if length(δts) == 0
        δt = 0
    else
        δt = mean(δts)
    end

    return PerturbedAtom(ϵ, δt, (atoms[1]).coord)
end

function Δ(scatterers::Vector{PerturbedAtom})
    n_atoms = length(scatterers)
    # Arrange the atoms into a matrix
    atoms_M = repeat(scatterers, 1, n_atoms)
    atoms_M_T = permutedims(atoms_M)
    # Check if there are modified bonds between pairs of atoms. If both atoms
    # have modified bonds with their neighbors, take the larger of the two
    Δ_ = map(
        (x, y) -> (x.coord in neighbors(y.coord)) * max(x.δt, y.δt),
        atoms_M,
        atoms_M_T,
    )
    # Add the on-site potential for the dopants
    Δ_ = Δ_ .+ Diagonal(map(x -> x.ϵ, scatterers))
    return Δ_
end

function Vj(scatterers::Vector{PerturbedAtom}, imp::ImpurityState)
    # Get the coordinates of the scattering atoms
    atom_coords = map(x -> x.coord, scatterers)
    # Get all the couplings from the impurity
    couplings = imp.coupling
    # Obtain the coupling strengths for each individual coupling
    Vs = map(x -> x.V, couplings)
    # ... and the coordinates of the graphene atoms that correspond to each V
    coords = map(x -> x.coord, couplings)

    return map(y -> sum((map(x -> y == x, coords)) .* Vs)::Float64, atom_coords)
end

function V(scatterers::Vector{PerturbedAtom}, imps::Vector{ImpurityState})
    Vj_vecs = map(ii -> Vj(scatterers, ii), imps)
    return reduce(hcat, Vj_vecs)::Array{Float64,2}
end
