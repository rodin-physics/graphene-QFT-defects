include("scattering.jl")

# Integrand used to calculate the local density in graphene
function δρ_Graphene_Integrand(
    loc::GrapheneCoord,
    z,
    imps::Vector{ImpurityState},
    occ_num::Vector{Float64},
    mod_atoms::Vector{PerturbedAtom},
)
    atoms = scattering_atoms(imps, mod_atoms)
    n_atoms = length(atoms)
    prop_mat = propagator_matrix(z, map(x -> x.coord, atoms))
    PropVectorR = map(x -> propagator(x.coord, loc, z), atoms)
    Δ_ = Δ(atoms)
    if length(imps) != 0
        V_ = V(atoms, imps)
        Γ0 = map((x, y) -> z - x.ϵ - x.U * y, imps, occ_num) |> Diagonal |> inv
        D =
            (Δ_ .+ V_ * Γ0 * transpose(V_)) * inv(
                Matrix{Int}(I, n_atoms, n_atoms) .-
                prop_mat * (Δ_ .+ V_ * Γ0 * transpose(V_)),
            )
    else
        D = Δ_ * inv(Matrix{Int}(I, n_atoms, n_atoms) .- prop_mat * Δ_)
    end
    return (transpose(PropVectorR)*D*PropVectorR)[1]

end

# Local density function in graphene. For T = 0, we integrate along the
# imaginary axis. Otherwise, we perform a contour integration to enclose the
# Matsubara frequencies.
function δρ_Graphene(loc, s::GrapheneSystem)
    μ = s.μ
    T = s.T
    imps = s.imps
    occ_num = s.occ_num
    mod_atoms = s.mod_atoms
    if T == 0
        res = quadgk(
            x ->
                real(δρ_Graphene_Integrand(loc, μ + 1im * x, imps, occ_num, mod_atoms)),
            0,
            Inf,
            rtol = ν,
        )
        return (res[1] / π)
    else
        res = quadgk(
            x ->
                -imag(δρ_Graphene_Integrand(
                    loc,
                    x + μ + 1im * η,
                    imps,
                    occ_num,
                    mod_atoms,
                )) * nF(x, T),
            -Inf,
            Inf,
            rtol = ν,
        )
        return (res[1] / π)
    end
end

# Integrand used to calculate the graphene-induced density in impurity states
function δρ_Impurity_Integrand(
    n,
    z,
    imps::Vector{ImpurityState},
    occ_num::Vector{Float64},
    mod_atoms::Vector{PerturbedAtom},
)
    atoms = scattering_atoms(imps, mod_atoms)
    n_atoms = length(atoms)
    Δ_ = Δ(atoms)
    V_ = V(atoms, imps)

    prop_mat = propagator_matrix(z, map(x -> x.coord, atoms))

    Λ =
        prop_mat .+
        prop_mat * Δ_ * inv(Matrix{Int}(I, n_atoms, n_atoms) .- prop_mat * Δ_) * prop_mat

    Γ0 = map((x, y) -> z - x.ϵ - x.U * y, imps, occ_num) |> Diagonal |> inv
    res =
        Γ0 *
        transpose(V_) *
        Λ *
        inv(Matrix{Int}(I, n_atoms, n_atoms) .- V_ * Γ0 * transpose(V_) * Λ) *
        V_ *
        Γ0

    return (res[n]::ComplexF64)
end

# Graphene-induced density in impurity states. For T = 0, we integrate along the
# imaginary axis. Otherwise, we perform a contour integration to enclose the
# Matsubara frequencies.
function δρ_Impurity(n, s::GrapheneSystem)
    μ = s.μ
    T = s.T
    imps = s.imps
    occ_num = s.occ_num
    mod_atoms = s.mod_atoms
    if T == 0
        res = quadgk(
            x -> real(δρ_Impurity_Integrand(n, μ + 1im * x, imps, occ_num, mod_atoms)),
            0,
            Inf,
            rtol = ν,
        )
        return (res[1] / π)
    else
        res = quadgk(
            x ->
                -imag(δρ_Impurity_Integrand(n, x + μ + 1im * η, imps, occ_num, mod_atoms)) * nF(x, T),
            -Inf,
            Inf,
            rtol = ν,
        )
        return (res[1] / π)
    end
end

# Total impurity density is the sum of the isolated impurity filling factor
# plus a correction δρ
function ρ_Impurity(n, s::GrapheneSystem)
    en = map((x, y) -> x.ϵ + x.U * y, s.imps, s.occ_num)
    return δρ_Impurity(n, s) + nF(en[n] - s.μ, s.T)
end

# Total impurity density for all the impurities
function ρ_All(s::GrapheneSystem)
    res = map(n -> ρ_Impurity(n, s), 1:length(s.imps))
end
