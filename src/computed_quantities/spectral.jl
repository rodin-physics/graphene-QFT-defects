include("density.jl")

function spectral_graphene(ω, loc, s::GrapheneSystem)
    pristine_spectral = -ω * Ω(ω + 1im * η, 0, 0) / π |> imag
    correction_spectral =
        -δρ_Graphene_Integrand(loc, ω + 1im * η, s.imps, s.mod_atoms) / π |>
        imag
    return (pristine_spectral + correction_spectral)
end

function spectral_graphene_17(ω, loc, s::GrapheneSystem)
    pristine_spectral =
        -propagator_17(graphene_A(0, 0), graphene_A(0, 0), ω + 1im * η) / π |>
        imag
    correction_spectral =
        -δρ_Graphene_Integrand_17(loc, ω + 1im * η, s.imps, s.mod_atoms) / π |>
        imag
    return (pristine_spectral + correction_spectral)
end

function spectral_impurity(ω, n, s::GrapheneSystem)
    Γ0_inv = map(x -> ω + 1im * η - x.ϵ, s.imps)
    pristine_spectral = -1 / (Γ0_inv[n]) / π |> imag
    correction_spectral =
        -δρ_Impurity_Integrand(n, ω + 1im * η, s.imps, s.mod_atoms) / π |> imag
    return (pristine_spectral + correction_spectral)
end
