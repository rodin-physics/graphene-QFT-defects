include("scattering.jl")

## Graphene Energy
function F_G_Integrand(
    z,
    imps::Vector{ImpurityState},
    occ_num::Vector{Float64},
    mod_atoms::Vector{PerturbedAtom},
)
    atoms = scattering_atoms(imps, mod_atoms)
    n_atoms = length(atoms)

    Δ_ = Δ(atoms)

    prop_mat = propagator_matrix(z, map(x -> x.coord, atoms))
    res = Matrix{Int}(I, n_atoms, n_atoms) .- prop_mat * Δ_ |> det |> log
    # The factor of 2 is due to the two spins
    return (-2 * res)
end

function F_G(s::GrapheneSystem)
    μ = s.μ
    T = s.T
    imps = s.imps
    occ_num = s.occ_num
    mod_atoms = s.mod_atoms
    if T == 0
        res = quadgk(
            x -> real(F_G_Integrand(μ + 1im * x, imps, occ_num, mod_atoms)),
            0,
            Inf,
            rtol = 1e-5,
        )
        return (res[1] / π)
    else
        # NOTE: THE ISSUE IN FINITE T SEEMS TO COME FROM THE FACT THAT LOG
        # IS MULTIVALUED!!!

        error("Finite T not implemented")
        # res = quadgk(
        #     x -> -imag(F_G_Integrand(
        #         x + μ + 1im * η,
        #         imps,
        #         occ_num,
        #         mod_atoms,
        #     )) * nF(x, T),
        #     -Inf,
        #     Inf,
        #     rtol = 1e-5,
        # )
        # return (res[1] / π)
    end
end
