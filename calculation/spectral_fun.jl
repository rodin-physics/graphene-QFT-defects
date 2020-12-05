include("../src/spectral.jl")
#
# sys = new_graphene_system()
# sys = set_μ(sys, 0.2)
# imp1 = new_impurity(2.5, 0.0)
# imp1 = add_coupling(imp1, -8.5, 0.5 * t, graphene_A(0,0))
# imp1 = add_coupling(imp1, -0.4, 0.0 * t, graphene_B(0,0))
# imp1 = add_coupling(imp1, -0.4, 0.0 * t, graphene_B(1,0))
# imp1 = add_coupling(imp1, -0.4, 0.0 * t, graphene_B(0,1))
# sys = add_imp(sys, imp1)
#
# @time δρ_Graphene(graphene_A(21,2), sys)
#
#
# @time ρ_Impurity(1, sys)
#
# spectral_graphene(0.1, graphene_B(2,0), sys)
#
#
# spectral_impurity(3.1, 1, sys)
#
#
#
#
# # ωs = range(-3 * t, 3 * t, length = 400)
# # # sp_fun = map(x -> spectral_impurity(x, 1, sys), ωs)
# # sp_fun = map(x -> spectral_graphene(x,atom_move(graphene_B(), 1, 1), sys), ωs)
# # # spectral_graphene(1.1, graphene_B(), sys)
# #
# # pyplot()
# # plot(ωs, sp_fun)
# # # savefig("T2.pdf")
# # savefig("T6.pdf")
# @time propagator(graphene_B(0,0), graphene_B(0,0), 0.1 + 1im * η)
