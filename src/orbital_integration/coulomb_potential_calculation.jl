using Distributed
numProc = 4
if nprocs() < numProc
    addprocs(numProc - nprocs())
end

# Calculate a fine grid of coulomb_potential_pz(...). This will be interpolated
# and used to calculate the interaction between pz orbital
@everywhere include("orbital_data_import.jl")
@everywhere nPts = 400;
@everywhere Rmax = 30;
@everywhere Rs = range(0, Rmax, length = nPts) .* ones(nPts)'
@everywhere τs = range(0, π / 2, length = nPts)' .* ones(nPts)

res = @showprogress pmap((R, τ) -> coulomb_potential_pz(R, τ)[1], Rs, τs)
writedlm("coulomb_potential_pz_precompute_Rmax_$Rmax.dat", res)
