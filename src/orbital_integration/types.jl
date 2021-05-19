include("../pristine/graphene_types.jl")
## Structures

#=
Data structure that contains the positions of ions, their unit cell indices,
and and their sublattice indices within the unit cell
=#
struct IonData
        numberatoms::Any
        sca1::Tuple
        sca2::Tuple
        unitcellist::Array{Tuple,1}
        sublatticelist::Array{Integer,1}
        ionpositionlist::Array{Tuple,1}
end
