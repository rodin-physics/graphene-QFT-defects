using Distributed

@everywhere include("calculation/computation_settings.jl")

# Calculation
@everywhere function f_A_Sublattice(x1, x2)
    Loc = graphene_A(x1, x2)
    res = δρ_Graphene(Loc, sys)
    return res[1]
end

@everywhere function f_B_Sublattice(x1, x2)
    Loc = graphene_B(x1, x2)
    res = δρ_Graphene(Loc, sys)
    return res[1]
end

resA = @showprogress pmap(f_A_Sublattice, D1S, D2S)
resB = @showprogress pmap(f_B_Sublattice, D1S, D2S)

writedlm("Data/Density/rho_A.dat", resA)
writedlm("Data/Density/rho_B.dat", resB)
