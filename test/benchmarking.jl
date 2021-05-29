using BenchmarkTools
include("../src/computed_quantities/spectral.jl")

## Propagator Integrals
@benchmark Ω(0.2 + 1im, rand(-20:20), rand(-20:20))
@benchmark Ωn(0.2 + 1im, rand(-20:20), rand(-20:20))
@benchmark Ωp(0.2 + 1im, rand(-20:20), rand(-20:20))

@benchmark propagator(
    graphene_A(rand(-20:20), rand(-20:20)),
    graphene_A(rand(-20:20), rand(-20:20)),
    1 + 1im,
)

@benchmark propagator(
    graphene_A(rand(-20:20), rand(-20:20)),
    graphene_B(rand(-20:20), rand(-20:20)),
    1 + 1im,
)


# Create a new system
my_system = new_graphene_system()

# Change the doping of the system
my_system = set_μ(my_system, -0.4)

pert_atom_1 = PerturbedAtom(-8, 0.0, graphene_A(0, 0))

my_system = add_mod_atom(my_system, pert_atom_1)
@benchmark spectral_graphene(2, graphene_A(0, 1), my_system)
@time δρ_Graphene_Integrand(
    graphene_A(0, 3),
    1im + 0,
    my_system.imps,
    my_system.mod_atoms,
)


@time δρ_Graphene(graphene_A(4, 1), my_system)
quadgk(
    x -> spectral_graphene(x, graphene_A(3, 0), my_system),
    -Inf,
    Inf,
    atol = 1e-4,
    rtol = 1e-3,
)
# Create a couple of impurities
imp1 = new_impurity(1.0)
imp2 = new_impurity(1.0)

# Add couplings to these impurities
imp1 = add_coupling(imp1, 2, graphene_A(0, 0))
imp2 = add_coupling(imp2, 2, graphene_A(2, 2))
my_system
# Add the impurities to the system
my_system = add_imp(my_system, imp1)
my_system = add_imp(my_system, imp2)


@benchmark spectral_graphene(2, graphene_A(0, 1), my_system)

@time δρ_Graphene_Integrand(
    graphene_A(0, 1),
    1im + 1,
    my_system.imps,
    my_system.mod_atoms,
)

@benchmark δρ_Graphene(graphene_A(0, 1), my_system)
δρ_Graphene(graphene_A(0, 0), my_system)
# spectral_graphene(ω, loc, s::GrapheneSystem)


@time Ω(2 + 1im, 10, 0)

# Create a new system
my_system = new_graphene_system()

# Change the doping of the system
my_system = set_μ(my_system, -0.4)

# Create a couple of impurities
imp1 = new_impurity(1.0)
imp2 = new_impurity(1.0)

# Add couplings to these impurities
imp1 = add_coupling(imp1, 2, graphene_A(0, 0))
imp2 = add_coupling(imp2, 2, graphene_A(2, 2))
my_system
# Add the impurities to the system
my_system = add_imp(my_system, imp1)
my_system = add_imp(my_system, imp2)

# Compute the local density in graphene
@time δρ_Graphene(graphene_A(1, 1), my_system)

# Remove the first impurity
my_system = remove_imp(my_system, 1)

# Compute the local density in graphene
δρ_Graphene(graphene_A(1, 1), my_system)

# Change the doping of the system and then recompute the local density
my_system = set_μ(my_system, -0.2)
δρ_Graphene(graphene_A(1, 1), my_system)

# Change the doping and the temperature
my_system = set_μ(my_system, 0.9)
my_system = set_T(my_system, 1e-2)

# Compute the density  at  the impurity
ρ_Impurity(1, my_system)
