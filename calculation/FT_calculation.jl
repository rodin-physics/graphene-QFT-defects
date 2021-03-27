@everywhere include("../src/analysis/dataprocess.jl")
## Data preparation
A_Lattice = readdlm("Data/Spectral/spectral_A.dat")
B_Lattice = readdlm("Data/Spectral/spectral_B.dat")
data = Data_Process(A_Lattice, B_Lattice)

XS = data[:, 1]
YS = data[:, 2]
ρ = data[:, 3]

ρ = ρ .- mean(ρ)

## Momenta
qx_min = -1 * π;
qx_max = 1 * π;

n_pts = 500;

qx = range(qx_min, qx_max, length = n_pts)
qy = range(qx_min, qx_max, length = n_pts)

QX = repeat(qx, 1, n_pts)
QY = repeat(qy, 1, n_pts) |> permutedims

## Calculation
FT_res = @showprogress pmap((qx, qy) -> FT_component(qx, qy, ρ, XS, YS), QX, QY)

writedlm("Data/FT/QX.dat", QX)
writedlm("Data/FT/QY.dat", QY)
writedlm("Data/FT/FT_real.dat", real.(FT_res))
writedlm("Data/FT/FT_imag.dat", imag.(FT_res))
