include("defect_types.jl")
include("../pristine/propagator.jl")

# A function that takes a list of impurities and another one of directly
# modified atoms and returns a list of all the scattering atoms
function scattering_atoms(imps::Vector{ImpurityState}, mod_atoms::Vector{PerturbedAtom})
    # Start by obtaining the atoms that are coupled to the impurities
    # Get all the couplings from the impurity states
    all_couplings = collect(Iterators.flatten(map(x -> x.coupling, imps)))
    # all_couplings = reduce(vcat, map(x -> x.coupling, imps))
    # Turn the couplings into PerturbedAtoms and get the unique hosts
    hosts = map(x -> PerturbedAtom(0.0, 0.0, x.coord), all_couplings) |> unique
    # Next, we get the neighbors of mod_atoms
    neighbor_coords =
        collect(Iterators.flatten(map(neighbors, map(x -> x.coord, mod_atoms))))
    all_neighbors = map(x -> PerturbedAtom(0.0, 0.0, x), neighbor_coords)
    # We combine the hosts with the neighbors and get a unique list
    hosts_neighbors = vcat(hosts, all_neighbors) |> unique
    # If the same coordinate appears in hosts_neighbors and in mod_atoms, take
    # the PerturbedAtom from mod_atoms. First, identify the mod_atoms coords
    mod_coords = map(x -> x.coord, mod_atoms)
    # Next, keep only the host coordinates that do NOT appear in the mod_coords
    hosts_neighbors = filter(x -> x.coord ∉ mod_coords, hosts_neighbors)
    # Now, we can combine the perturbed and host atoms
    all_perturbed_atoms = vcat(hosts_neighbors, mod_atoms)
    return all_perturbed_atoms
end

function Δ(scatterers::Vector{PerturbedAtom})
    n_atoms = length(scatterers)
    # Arrange the atoms into a matrix
    atoms_M = repeat(scatterers, 1, n_atoms)
    atoms_M_T = permutedims(atoms_M)
    # Check if there are modified bonds between pairs of atoms. If both atoms
    # have modified bonds with their neighbors, take the larger of the two
    Δ_ =
        map((x, y) -> (x.coord in neighbors(y.coord)) * max(x.δt, y.δt), atoms_M, atoms_M_T)
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
