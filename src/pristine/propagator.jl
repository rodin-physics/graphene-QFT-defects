include("graphene_types.jl")
include("tight_binding.jl")

@inline function Ω_Integrand(z, u, v, x::Float64)
    W = ((z / t)^2 - 1.0) / (4.0 * cos(x)) - cos(x)
    return (
        exp(1.0im * (u - v) * x) / cos(x) *
        ((W - √(W - 1) * √(W + 1))^abs.(u + v)) / (√(W - 1) * √(W + 1))
    )
end

@inline function Ω(z, u, v)
    return ((quadgk(
        x -> Ω_Integrand(z, u, v, x) / (8.0 * π * t^2),
        0.0,
        2.0 * π,
    ))[1])
end

@inline function Ωp_Integrand(z, u, v, x::Float64)
    W = ((z / t)^2 - 1.0) / (4.0 * cos(x)) - cos(x)
    return (
        2 *
        exp(1.0im * (u - v) * x) *
        ((W - √(W - 1) * √(W + 1))^abs.(u + v + 1)) / (√(W - 1) * √(W + 1))
    )
end

@inline function Ωp(z, u, v)
    return ((quadgk(
        x -> Ωp_Integrand(z, u, v, x) / (8.0 * π * t^2),
        0.0,
        2.0 * π,
    ))[1])
end

@inline function Ωn_Integrand(z, u, v, x::Float64)
    W = ((z / t)^2 - 1.0) / (4.0 * cos(x)) - cos(x)
    return (
        2 *
        exp(1.0im * (u - v) * x) *
        ((W - √(W - 1) * √(W + 1))^abs.(u + v - 1)) / (√(W - 1) * √(W + 1))
    )
end

@inline function Ωn(z, u, v)
    return ((quadgk(
        x -> Ωn_Integrand(z, u, v, x) / (8.0 * π * t^2),
        0.0,
        2.0 * π,
    ))[1])
end

# The propagator function picks out the correct element of the Ξ matrix based
# on the sublattices of the graphene coordinates
function propagator(a_l::GrapheneCoord, a_m::GrapheneCoord, z)
    u = a_l.u - a_m.u
    v = a_l.v - a_m.v
    if a_l.sublattice == a_m.sublattice
        return (z * Ω(z, u, v))
    elseif ([a_l.sublattice, a_m.sublattice] == ["●", "○"])
        return (-t * (Ω(z, u, v) + Ωp(z, u, v)))
    elseif ([a_l.sublattice, a_m.sublattice] == ["○", "●"])
        return (-t * (Ω(z, u, v) + Ωn(z, u, v)))
    else
        error("Illegal sublattice parameter")
    end
end

# The (I^T Ξ I) Matrix. We use the fact that the matrix is symmetric to speed
# up the calculation
function propagator_matrix(z, Coords::Vector{GrapheneCoord})
    len_coords = length(Coords)
    out = zeros(ComplexF64, len_coords, len_coords)
    for ii = 1:len_coords
        @inbounds for jj = ii:len_coords
            out[ii, jj] = propagator(Coords[ii], Coords[jj], z)
            out[jj, ii] = out[ii, jj]
        end
    end
    return out
end

# The propagator function for the 17-parameter TB
function propagator_17_same_sublattice(u, v, z)
    function integrand(x, f)
        tmp =
            exp(-1im * π * (u - v)) *
            Gπ(z, [4 * π * x[1] / d - 2 * π / d, 4 * π * x[2] / (√(3) * d)])[
                1,
                1,
            ] *
            exp(2im * π * ((u - v) * x[1] + (u + v) * x[2]))
        f[1] = real(tmp)
        f[2] = imag(tmp)
    end
    res = cuhre(integrand, 2, 2, rtol = ν, maxevals = nevals)
end

function propagator_17_AB(u, v, z)
    function integrand(x, f)
        tmp =
            exp(-1im * π * (u - v)) *
            Gπ(z, [4 * π * x[1] / d - 2 * π / d, 4 * π * x[2] / (√(3) * d)])[
                1,
                2,
            ] *
            exp(2im * π * ((u - v) * x[1] + (u + v) * x[2]))
        f[1] = real(tmp)
        f[2] = imag(tmp)
    end
    res = cuhre(integrand, 2, 2, rtol = ν, maxevals = nevals)
end

function propagator_17_BA(u, v, z)
    function integrand(x, f)
        tmp =
            exp(-1im * π * (u - v)) *
            Gπ(z, [4 * π * x[1] / d - 2 * π / d, 4 * π * x[2] / (√(3) * d)])[
                2,
                1,
            ] *
            exp(2im * π * ((u - v) * x[1] + (u + v) * x[2]))
        f[1] = real(tmp)
        f[2] = imag(tmp)
    end
    res = cuhre(integrand, 2, 2, rtol = ν, maxevals = nevals)
end

function propagator_17(a_l::GrapheneCoord, a_m::GrapheneCoord, z)
    u = a_l.u - a_m.u
    v = a_l.v - a_m.v
    if a_l.sublattice == a_m.sublattice
        res = propagator_17_same_sublattice(u, v, z)
    elseif ([a_l.sublattice, a_m.sublattice] == ["●", "○"])
        res = propagator_17_AB(u, v, z)
    elseif ([a_l.sublattice, a_m.sublattice] == ["○", "●"])
        res = propagator_17_BA(u, v, z)
    else
        error("Illegal sublattice parameter")
    end
    return (res[1][1] + 1im * res[1][2])
end

function propagator_matrix_17(z, Coords::Vector{GrapheneCoord})
    len_coords = length(Coords)
    out = zeros(ComplexF64, len_coords, len_coords)
    for ii = 1:len_coords
        @inbounds for jj = ii:len_coords
            out[ii, jj] = propagator_17(Coords[ii], Coords[jj], z)
            out[jj, ii] = out[ii, jj]
        end
    end
    return out
end
