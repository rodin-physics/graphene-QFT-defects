using Distributed

procs = 4;
if nprocs() < procs
    addprocs(procs - nprocs())
end

@everywhere include("computation_settings.jl")

# Calculation
@everywhere function f_A_Sublattice(x1, x2)
    Loc = graphene_A(x1, x2)
    res = spectral_graphene(ω_slice, Loc, sys)
    return res[1]
end

@everywhere function f_B_Sublattice(x1, x2)
    Loc = graphene_B(x1, x2)
    res = spectral_graphene(ω_slice, Loc, sys)
    return res[1]
end

resA = @showprogress pmap(f_A_Sublattice, D1S, D2S)
resB = @showprogress pmap(f_B_Sublattice, D1S, D2S)

writedlm("Data/Spectral/spectral_A.dat", resA)
writedlm("Data/Spectral/spectral_B.dat", resB)
