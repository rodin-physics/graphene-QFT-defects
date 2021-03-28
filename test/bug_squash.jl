res = map(n -> spectral_graphene(0.4, graphene_B(n, n), sys), 1:70)
res = map(
    n ->
        -δρ_Graphene_Integrand(
            graphene_B(n, n),
            0.6 + 1im * η,
            sys.imps,
            sys.mod_atoms,
        ) / π |> imag,
    1:70,
)

res = map(n -> Ω(0.4 + 1im * η, n, n), 0:70)

res
scatter(0:70, real.(res))

@time spectral_graphene(0.6, graphene_B(27, 27), sys)
@time spectral_graphene(0.6, graphene_B(28, 28), sys)
@time spectral_graphene(0.6, graphene_B(29, 29), sys)
@time spectral_graphene(0.6, graphene_B(30, 30), sys)

spectral_graphene(0.6, graphene_B(-27, -27), sys)
spectral_graphene(0.6, graphene_B(-28, -28), sys)
spectral_graphene(0.6, graphene_B(-29, -29), sys)
spectral_graphene(0.6, graphene_B(-30, -30), sys)

propagator(graphene_A(0, 0), graphene_A(-29, -29), 0.6 + 1im * η)
Ω(0.4 + 1im * η, 41, 41)
Ω(0.4 + 1im * η, 42, 42)
Ω(0.4 + 1im * η, 43, 43)
