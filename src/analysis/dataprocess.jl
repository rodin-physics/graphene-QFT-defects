include("../pristine/general.jl")

## Data Processing Functions

# In order to plot the calculated results as a lattice, one needs to create
# the correct coordinate arrays
function Data_Process(A_Lattice, B_Lattice)
    sz = size(A_Lattice)
    nPts = floor(Int, (sz[1] - 1) / 2)
    # Lattice shift between the two sublattices
    latticeShift = -1 / √(3) * 2.46
    # Arrays
    d1s = -nPts:1:nPts # d1 vectors
    d2s = -nPts:1:nPts # d2 vectors

    D1S = repeat(d1s, 1, 2 * nPts + 1)
    D2S = repeat(d2s', 2 * nPts + 1, 1)

    # Coordinates of the carbon atoms
    XS = d1[1] .* D1S .+ d2[1] .* D2S
    YS = d1[2] .* D1S .+ d2[2] .* D2S

    # Flatten the coordinates and the data
    XS_A = reshape(XS, 1, sz[1]^2)
    YS_A = reshape(YS, 1, sz[1]^2)
    A_Lattice = reshape(A_Lattice, 1, sz[1]^2)

    XS_B = reshape(XS, 1, sz[1]^2)
    YS_B = reshape(YS, 1, sz[1]^2) .+ latticeShift
    B_Lattice = reshape(B_Lattice, 1, sz[1]^2)

    # Combine all the coordinates and data for both sublattices
    XS = vcat(XS_A', XS_B')
    YS = vcat(YS_A', YS_B')
    dta = vcat(A_Lattice', B_Lattice')

    return [XS YS dta]
end

# Fourier Transform

function FT_component(qx, qy, ρ, XS, YS)
    res = sum(map((x, y, z) -> exp(1im * (x * qx + y * qy)) * z, XS, YS, ρ))
end
