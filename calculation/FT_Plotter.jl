using CairoMakie
using DelimitedFiles

QX = readdlm("Data/FT/QX.dat")
QY = readdlm("Data/FT/QY.dat")
FT_real = readdlm("Data/FT/FT_real.dat")
FT_imag = readdlm("Data/FT/FT_imag.dat")

FT_res = FT_real + 1im .* FT_imag

fig = Figure(resolution = (1800, 1800))
ax =
        fig[1, 1] = Axis(
                fig,
                xlabel = "qx",
                ylabel = "qy",
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

hm = heatmap!(qx, qy, angle.(FT_res), colormap = :oslo)
fig
save("FT_phase_04eV_1N.png", fig)
