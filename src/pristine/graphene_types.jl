include("general.jl")

## Structures and functions for defect-free graphene
# Graphene Coordinate with the position R = u * d1 + v * d2 and the sublattice
struct GrapheneCoord
    u::Int
    v::Int
    sublattice::String
end

# Location structure in Angstroms
struct Location
    x::Float64
    y::Float64
end

# Functions to initialize graphene coordinates
function graphene_A(u, v)
    return GrapheneCoord(u, v, "●")
end

function graphene_B(u, v)
    return GrapheneCoord(u, v, "○")
end

# Function for determining the neighboring atoms of a lattice site
function neighbors(atom::GrapheneCoord)
    u = atom.u
    v = atom.v
    if atom.sublattice == "●"
        return [
            graphene_B(u, v)
            graphene_B(u + 1, v)
            graphene_B(u, v + 1)
        ]
    elseif atom.sublattice == "○"
        return [
            graphene_A(u, v)
            graphene_A(u - 1, v)
            graphene_A(u, v - 1)
        ]
    else
        error("Illegal sublattice parameter")
    end
end

# Switching sublattice
function lattice_swap(atom::GrapheneCoord)
    if atom.sublattice == "●"
        return graphene_B(atom.u, atom.v)
    elseif atom.sublattice == "○"
        return graphene_A(atom.u, atom.v)
    else
        error("Illegal sublattice parameter")
    end
end

# Moving the graphene coordinate
function atom_move(atom::GrapheneCoord, u, v)
    return GrapheneCoord(atom.u + u, atom.v + v, atom.sublattice)
end

# Converting Crystal coordinates to Cartesian in Angstroms
function crystal_to_cartesian(coord::GrapheneCoord)
    u = coord.u
    v = coord.v
    x = d1[1] * u + d2[1] * v
    y = d1[2] * u + d2[2] * v
    if coord.sublattice == "●"
        shft = 0.0
    elseif coord.sublattice == "○"
        shft = lattice_shift
    else
        error("Illegal sublattice parameter")
    end

    return (Location(x, y + shft))
end
