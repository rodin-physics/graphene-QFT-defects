using DelimitedFiles
using LaTeXStrings
using Plots

include("dataprocess.jl")

A_Lattice = readdlm("Data/Density/rho_A.dat")
B_Lattice = readdlm("Data/Density/rho_B.dat")

data = Data_Process(A_Lattice, B_Lattice)

XS = data[:, 1]
YS = data[:, 2]
ρ = data[:, 3]

function FT_component(qx, qy)
    res = sum(map((x, y, z) -> exp(1im * (x * qx + y * qy)) * z, XS, YS, ρ))
end

qx_min = -1 * π;
qx_max = 1 * π;

n_pts = 100;

qx = range(qx_min, qx_max, length = n_pts)
qy = range(qx_min, qx_max, length = n_pts)

QX = repeat(qx, 1, n_pts)
QY = repeat(qy, 1, n_pts) |> permutedims

FT_res = map((qx, qy) -> FT_component(qx, qy), QX, QY)

gr()
heatmap(qx, qy, abs.(FT_res), aspect_ratio = 1)
# savefig("Phase.pdf")
