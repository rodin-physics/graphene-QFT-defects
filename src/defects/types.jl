include("general.jl")


# Atom that is perturbed either by an impurity state or via substitution. This
# changes the atom's on-site energy and its coupling to its 3 neighbors.
struct PerturbedAtom
    ϵ::Float64
    δt::Float64
    coord::GrapheneCoord
end

# Coupling for impurity states
struct Coupling
    V::Float64              # Coupling to a graphene atom
    δt::Float64             # The modification to the hopping parameter for the
    # graphene atom coupled to the impurity
    coord::GrapheneCoord    # Location of the graphene atom
end

# Impurity state defined by its on-site energy and its coupling to graphene
struct ImpurityState
    ϵ::Float64                  # Impurity state energy
    U::Float64                  # Hubbard term

    coupling::Vector{Coupling}  # Coupling array for the impurity
end

# The structure describing the system
struct GrapheneSystem
    μ::Float64                          # Chemical potential
    T::Float64                          # Temperature
    imps::Vector{ImpurityState}         # Impurity states in the system
    occ_num::Vector{Float64}            # Occupation numbers of impurity states
    # by the opposite spin
    mod_atoms::Vector{PerturbedAtom}    # Atoms modified by substitution or
    # external potential
end

## Manipulation functions

# Setting up the graphene system
function new_graphene_system()
    return GrapheneSystem(0.0, 0.0, ImpurityState[], Float64[], PerturbedAtom[])
end

function set_T(s::GrapheneSystem, new_T::Float64)
    return GrapheneSystem(s.μ, new_T, s.imps, s.occ_num, s.mod_atoms)
end

function set_μ(s::GrapheneSystem, new_μ::Float64)
    return GrapheneSystem(new_μ, s.T, s.imps, s.occ_num, s.mod_atoms)

end

function add_mod_atom(s::GrapheneSystem, new_mod_atom::PerturbedAtom)
    return GrapheneSystem(s.μ, s.T, s.imps, s.occ_num, push!(s.mod_atoms, new_mod_atom))
end

function add_imp(s::GrapheneSystem, new_imp::ImpurityState)
    return GrapheneSystem(
        s.μ,
        s.T,
        push!(s.imps, new_imp),
        push!(s.occ_num, 0.0),
        s.mod_atoms,
    )
end

function remove_mod_atom(s::GrapheneSystem, ind::Int)
    return GrapheneSystem(s.μ, s.T, s.imps, s.occ_num, deleteat!(s.mod_atoms, ind))
end

function remove_imp(s::GrapheneSystem, ind::Int)
    return GrapheneSystem(
        s.μ,
        s.T,
        deleteat!(s.imps, ind),
        deleteat!(s.occ_num, ind),
        s.mod_atoms,
    )
end

function new_impurity(ϵ, U)
    return ImpurityState(ϵ, U, Coupling[])
end

function add_coupling(imp::ImpurityState, V, δt, coord)
    new_coup = Coupling(V, δt, coord)
    if (coord in map(x -> x.coord, imp.coupling))
        error("Impurity is already coupled to this atom")
    else
        return ImpurityState(imp.ϵ, imp.U, push!(imp.coupling, new_coup))
    end
end

function remove_coupling(imp::ImpurityState, ind::Int)
    return ImpurityState(imp.ϵ, imp.U, deleteat!(imp.coupling, ind))
end
