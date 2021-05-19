include("../pristine/graphene_types.jl")
using Cubature
using Interpolations
# Import the data
radial_orb = readdlm("src/orbital_integration/radial_orbital.dat")
rs_DFT = pushfirst!(float.(radial_orb[2:end, 1]), 0)
orb_2p_radial_DFT = pushfirst!(float.(radial_orb[2:end, 2]), 0)
orb_2s_radial_DFT = pushfirst!(float.(radial_orb[2:end, 3]), 0)
orb_1s_radial_DFT = pushfirst!(float.(radial_orb[2:end, 4]), 0)

# Interpolate the data
orb_2p_radial = LinearInterpolation(rs_DFT, orb_2p_radial_DFT)
orb_2s_radial = LinearInterpolation(rs_DFT, orb_2s_radial_DFT)
orb_1s_radial = LinearInterpolation(rs_DFT, orb_1s_radial_DFT)

# 2pz orbital function. r is measured from the center of the orbital
function Ψ_pz(r::Vector{Float64})
    radius = norm(r)
    if radius == 0.0 || radius > 40
        return 0.0
    else
        # Spherical harmonic for 2pz (n=2, l=1, m=0): cos(theta) = z/r
        # Note that the DFT values for the radial portion are r R(r)
        # so we divide by an extra factor of r when computing the WF
        return (orb_2p_radial(radius) * r[3] / radius^2 * √(3 / (4 * π)))
    end
end

# Coulomb energy due to a proton at distance R (in Bohr) from the 2pz orbital
function coulomb_energy_proton(R::Vector{Float64})
    res = hcubature(
        x -> Ψ_pz(x) .^ 2 ./ (norm(x - R) + 1e-9),
        -10 * ones(3),
        10 * ones(3),
        reltol = 1e-4,
    )
end

# Coulomb energy due to the interaction between two 2pz orbitals at distance R
function coulomb_energy_pz(R::Vector{Float64})
    res = hcubature(
        x ->
            Ψ_pz(x[1:3]) .^ 2 .* Ψ_pz(x[4:6]) .^ 2 ./
            (norm(x[1:3] - x[4:6] + R) + 1e-6),
        -10 * ones(6),
        10 * ones(6),
        reltol = 1e-2,
    )
end
