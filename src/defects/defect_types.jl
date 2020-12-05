include("../pristine/propagator.jl")

## Structures describing defects
# Coupling for impurity states
struct Coupling
    V::Float64              # Coupling to a graphene atom
    coord::GrapheneCoord    # Location of the graphene atom
end

# Impurity state defined by its on-site energy and its coupling to graphene
struct ImpurityState
    ϵ::Float64                  # Impurity state energy
    coupling::Vector{Coupling}  # Coupling array for the impurity
end

# Atom that is perturbed either by an impurity state or via substitution. This
# changes the atom's on-site energy and its coupling to its 3 neighbors.
struct PerturbedAtom
    ϵ::Float64
    δt::Float64
    coord::GrapheneCoord
end

# The structure describing the system
struct GrapheneSystem
    μ::Float64                          # Chemical potential
    T::Float64                          # Temperature
    imps::Vector{ImpurityState}         # Impurity states in the system
    mod_atoms::Vector{PerturbedAtom}    # Atoms modified by substitution or
    # external potential
end

## Manipulation functions
# Setting up the graphene system
function new_graphene_system()
    return GrapheneSystem(0.0, 0.0, ImpurityState[], PerturbedAtom[])
end

# Changing the temperature of an initialized system
function set_T(s::GrapheneSystem, new_T::Float64)
    return GrapheneSystem(s.μ, new_T, s.imps, s.mod_atoms)
end

# Changing the chemical potential of an initialized system
function set_μ(s::GrapheneSystem, new_μ::Float64)
    return GrapheneSystem(new_μ, s.T, s.imps, s.mod_atoms)
end

# Adding a new modified atom to an initialized system
function add_mod_atom(s::GrapheneSystem, new_mod_atom::PerturbedAtom)
    return GrapheneSystem(s.μ, s.T, s.imps, push!(s.mod_atoms, new_mod_atom))
end

# Adding a new impurity to an initialized system
function add_imp(s::GrapheneSystem, new_imp::ImpurityState)
    return GrapheneSystem(s.μ, s.T, push!(s.imps, new_imp), s.mod_atoms)
end

# Removing a modified atom from a system
function remove_mod_atom(s::GrapheneSystem, ind::Int)
    return GrapheneSystem(s.μ, s.T, s.imps, deleteat!(s.mod_atoms, ind))
end

# Removing an impurity from a system
function remove_imp(s::GrapheneSystem, ind::Int)
    return GrapheneSystem(s.μ, s.T, deleteat!(s.imps, ind), s.mod_atoms)
end

# Initialize an impurity
function new_impurity(ϵ)
    return ImpurityState(ϵ, Coupling[])
end

# Add coupling to an existing impurity
function add_coupling(imp::ImpurityState, V, coord)
    new_coup = Coupling(V, coord)
    if (coord in map(x -> x.coord, imp.coupling))
        error("Impurity is already coupled to this atom")
    else
        return ImpurityState(imp.ϵ, push!(imp.coupling, new_coup))
    end
end

# Remove coupling from an impurity
function remove_coupling(imp::ImpurityState, ind::Int)
    return ImpurityState(imp.ϵ, deleteat!(imp.coupling, ind))
end
