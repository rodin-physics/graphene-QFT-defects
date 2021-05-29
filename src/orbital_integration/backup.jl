



#
#
#
#
#
# function coulomb_energy_pz_new(ρ::Vector{Float64})
#     hcubature(
#         r ->
#             Ψ_pz(r[1], r[2]) .^ 2 *
#             r[1]^2 *
#             sin(r[2]) *
#             coulomb_potential_pz(norm(ρ - r), acos(((ρ-r)[3]) / norm(ρ - r)))[1],
#         [0, 0, 0],
#         [20, π, 2 * π],
#         reltol = 1e-2,
#     )
# end
#
# coulomb_energy_pz_new([5e0, 0, 0])
# coulomb_potential_pz(5e0, 0)
# # acos(((ρ-r)[3]) / norm(ρ - r)
# r = 1e-5
# coulomb_energy_pz([1e-3, 0, 0])
#
#
# ρ = [10, 2, 1.0]
# r = [0.1, 2, 11]
# Ψ_pz_new(r[1], r[2]) .^ 2 *
# r[1]^2 *
# sin(r[2]) *
# coulomb_potential_pz(norm(ρ - r), 0)
#
#
# @time coulomb_potential_pz(r, pi / 2)[1]
#
# @benchmark coulomb_potential_pz(norm(rr), acos(rr[3] / norm(rr)))
# @benchmark coulomb_potential_pz_new(norm(rr), acos(rr[3] / norm(rr)))
#
#
#
#
# function Ψ_pz(r::Vector{Float64})
#     radius = norm(r)
#     if radius == 0.0 || radius > 40
#         return 0.0
#     else
#         # Spherical harmonic for 2pz (n=2, l=1, m=0): cos(theta) = z/r
#         # Note that the DFT values for the radial portion are r R(r)
#         # so we divide by an extra factor of r when computing the WF
#         return (orb_2p_radial(radius) * r[3] / radius^2 * √(3 / (4 * π)))
#     end
# end
#
#
#
#
# r = [1.2, 4, 1]
# nr = norm(r)
# θ = acos(r[3] / nr)
# @time Ψ_pz(r)
# @time Ψ_pz_new(nr, 1.0)
#
# acos(1 / norm([1.2, 4, 1]))
# # Coulomb energy due to a proton at distance R (in Bohr) from the 2pz orbital
# function coulomb_energy_proton(R::Vector{Float64})
#     res = hcubature(
#         x -> Ψ_pz(x) .^ 2 ./ (norm(x - R) + 1e-9),
#         -22 * ones(3),
#         22 * ones(3),
#         reltol = 1e-4,
#     )
# end
#
# # Coulomb energy due to a proton at distance R (in Bohr) from the 2pz orbital
# function coulomb_potential_pz(R, τ)
#     res = hcubature(
#         x ->
#             Ψ_pz_new(x[1], x[2]) .^ 2 * x[1]^2 * sin(x[2]) ./ (
#                 √(
#                     x[1]^2 + R^2 -
#                     2 *
#                     x[1] *
#                     R *
#                     (cos(τ) * cos(x[2]) + cos(x[3]) * sin(τ) * sin(x[2])),
#                 ) + 1e-9
#             ),
#         [0, 0, 0],
#         [22, π, 2 * π],
#         reltol = 1e-4,
#     )
# end
#
# function coulomb_potential_pz_new(R, τ)
#     pt = R .* [sin(τ), 0, cos(τ)]
#     res = hcubature(
#         x -> Ψ_pz(x) .^ 2 ./ (norm(x - pt) + 1e-9),
#         -22 * ones(3),
#         22 * ones(3),
#         reltol = 1e-4,
#     )
# end
#
# rr = [1.0, 1, 1]
# @time coulomb_energy_proton(rr + [0, 0, 1e-4])
# @time coulomb_potential_pz(norm(rr), acos(rr[3] / norm(rr)))
#
# @time coulomb_energy_proton(rr + [0, 0, 1e-4])
# @time coulomb_potential_pz(0, 0)
#
#
# # Coulomb energy due to the interaction between two 2pz orbitals at distance R
# function coulomb_energy_pz(R::Vector{Float64})
#     res = hcubature(
#         x ->
#             Ψ_pz(x[1:3]) .^ 2 .* Ψ_pz(x[4:6]) .^ 2 ./
#             (norm(x[1:3] - x[4:6] + R) + 1e-6),
#         -7 * ones(6),
#         7 * ones(6),
#         reltol = 0.5e-2,
#     )
# end
# hcubature(
#     x -> Ψ_pz(x) .^ 2 ./ (norm(x - [0.0, 0, 0]) + 1e-9),
#     -30 * ones(3),
#     30 * ones(3),
#     reltol = 1e-4,
# )
# hcubature(
#     x -> Ψ_pz(x) .^ 2 ./ (norm(x - [0.0, 0, 0]) + 1e-9),
#     -20 * ones(3),
#     20 * ones(3),
#     reltol = 1e-4,
# )
# hcubature(
#     x -> Ψ_pz(x) .^ 2 ./ (norm(x - [0.0, 0, 0]) + 1e-9),
#     -10 * ones(3),
#     10 * ones(3),
#     reltol = 1e-4,
# )
#
# Ψ_pz([7.0, 7, 7])
# coulomb_energy_pz(NN_Coords_Cartesian_3D[1] + [0, 0, 1e-2])
# coulomb_energy_pz(NN_Coords_Cartesian_3D[1] + [0, 0, 1e-2])
# @time coulomb_energy_proton(NN_Coords_Cartesian_3D[17] + [0, 0, 1e-4])
#
# 1 / norm([0, 0, 10.0])



#
# proton_energies =
#     @showprogress map(coulomb_energy_proton, NN_Coords_Cartesian_3D)
# ee_energies =
#     @showprogress map(coulomb_energy_pz, NN_Coords_Cartesian_3D[2:end])
#
# proton_energies_approx = 1 ./ norm.(NN_Coords_Cartesian_3D)
# proton_energies_approx - map(x -> x[1], proton_energies)
# @time coulomb_energy_proton(NN_Coords_Cartesian_3D[1])
# coulomb_energy_pz(NN_Coords_Cartesian_3D[1] + [0, 0, 1e-2])
# coulomb_energy_pz(NN_Coords_Cartesian_3D[1] + [0, 0, 1e-2])
#
# .65*27
#
# rr =NN_Coords_Cartesian_3D[2]
# @time coulomb_potential_pz(norm(rr), pi/2)
# NN_Coords_Cartesian_3D
