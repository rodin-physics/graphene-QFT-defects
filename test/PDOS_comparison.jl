using DelimitedFiles
data = readdlm("test/uc.P-Cp.pdos")

# DFT data
ens = data[:, 1] .+ 1.7668
pdos = data[:, 3]

# QFT data
ωs = range(-10, 15, length = 300)
sp_fun = map(
        x -> spectral_graphene(x, graphene_A(0, 0), new_graphene_system()),
        ωs,
)

fig = Figure(resolution = (1800, 1800))
ax =
        fig[1, 1] = Axis(
                fig,
                xlabel = "E",
                ylabel = "DOS",
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

l = lines!(ens, pdos)
l2 = lines!(ωs, 5 .*sp_fun, color = :red)
xlims!(ax, (-10, 15))
fig
