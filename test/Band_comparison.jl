using CairoMakie
using DelimitedFiles

include("../src/pristine/general.jl")

## Data
DFT_Bands_G = readdlm("test/graphene.pbe.bands.dat")

# Distances between different point of symmetry
ΓK_L = √((Γ[1] - K[1])^2 + (Γ[2] - K[2])^2)
KM_L = √((K[1] - M[1])^2 + (K[2] - M[2])^2)
MΓ_L = √((M[1] - Γ[1])^2 + (M[2] - Γ[2])^2)

# BZ path: Γ-K-M-Γ
nPts = 101;

# x coordinates used for plotting the bands
xCoord = vcat(
        (1:nPts) / nPts * ΓK_L,
        ΓK_L .+ (1:nPts) / nPts * KM_L,
        ΓK_L .+ KM_L .+ (1:nPts) / nPts * MΓ_L,
)

# Momenta
QX = vcat(
        range(Γ[1], stop = K[1], length = nPts),
        range(K[1], stop = M[1], length = nPts),
        range(M[1], stop = Γ[1], length = nPts),
)

QY = vcat(
        range(Γ[2], stop = K[2], length = nPts),
        range(K[2], stop = M[2], length = nPts),
        range(M[2], stop = Γ[2], length = nPts),
)

# Calculated energies along the predefined momentum path
energy_states =
        map((qx, qy) -> inv(Overlap([qx, qy])) * Hπ([qx, qy]) |> eigen, QX, QY)
Energies = reduce(hcat, map(x -> x.values, energy_states))

fig = Figure(resolution = (1800, 1800))
ax =
        fig[1, 1] = Axis(
                fig,
                # xlabel = "E",
                # ylabel = "DOS",
                xlabelpadding = 0,
                ylabelpadding = 0,
                xlabelsize = 12,
                ylabelsize = 12,
                xticklabelsize = 12,
                yticklabelsize = 12,
                aspect = AxisAspect(1),
                xticklabelfont = "serif-roman",
                yticklabelfont = "serif-roman",
                xlabelfont = "serif-italic",
                ylabelfont = "serif-roman",
        )

sc = scatter!(
        DFT_Bands_G[:, 1] ./ maximum(DFT_Bands_G[:, 1]) .* (ΓK_L + MΓ_L + KM_L),
        DFT_Bands_G[:, 2] .+ 0.7668,
        strokewidth = 0,
        marker = :circle,
        markersize = 2.0,
)

l1 = lines!(xCoord, real.(Energies[1, :]))
l2 = lines!(xCoord, real.(Energies[2, :]))

ylims!(ax, (-10, 15))
# ylims!(ax, (-200, 200))
fig
