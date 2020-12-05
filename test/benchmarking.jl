using BenchmarkTools

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

# Add the impurities to the system
my_system = add_imp(my_system, imp1)
my_system = add_imp(my_system, imp2)

# Compute the local density in graphene
δρ_Graphene(graphene_A(1, 1), my_system)

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
