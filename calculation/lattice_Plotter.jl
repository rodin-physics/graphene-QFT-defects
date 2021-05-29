using CairoMakie

include("../src/analysis/dataprocess.jl")

A_Lattice = readdlm("Data/Spectral/spectral_A.dat")
B_Lattice = readdlm("Data/Spectral/spectral_B.dat")

data = Data_Process(A_Lattice, B_Lattice)

XS = data[:, 1]
YS = data[:, 2]
ρ = data[:, 3]

fig = Figure(resolution = (1800, 1800))
ax =
        fig[1, 1] = Axis(
                fig,
                xlabel = "x/Å",
                ylabel = "y/Å",
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
                ylabelfont = "serif-italic",
        )

sc = scatter!(
        ax,
        XS,
        YS,
        color = ρ .- ρ[1],
        strokewidth = 0,
        marker = :circle,
        markersize = 1.0,
        colormap = :coolwarm,
        colorrange = (-.001, .001),
)
xlims!(ax, (-200, 200))
ylims!(ax, (-200, 200))
fig

save("Spectral_04eV_1N.png", fig)
