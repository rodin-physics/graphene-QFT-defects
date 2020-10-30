include("density.jl")

function spectral_graphene(ω, loc, s::GrapheneSystem)
    pristine_spectral = -ω * Ω(ω + 1im * η, 0, 0) / π |> imag
    correction_spectral =
        -δρ_Graphene_Integrand(
            loc,
            ω + 1im * η,
            s.imps,
            s.occ_num,
            s.mod_atoms,
        ) / π |> imag
    return (pristine_spectral + correction_spectral)
end

function spectral_impurity(ω, n, s::GrapheneSystem)
    Γ0_inv = map((x, y) -> ω + 1im * η - x.ϵ - x.U * y, s.imps, s.occ_num)
    pristine_spectral = -1 / (Γ0_inv[n]) / π |> imag
    correction_spectral =
        -δρ_Impurity_Integrand(n, ω + 1im * η, s.imps, s.occ_num, s.mod_atoms) /
        π |> imag
    return (pristine_spectral + correction_spectral)
end
