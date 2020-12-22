# Let F(β) be the NV-dimensional cumulative distribution of coefficient vector β in the population.
# We assume F(β) can be arbitrarily approximated by a finite support set S, with β_r ∈ S.
# We let Pr(β = β_r) = exp(z(β_r)'ω) / Σ_r exp(z(β_r)'ω); then Pr(choose i) = Σ_r Pr(choose i | β_r) × Pr(β = β_r).
# BETAS holds the draws b_r from the support of β_r for each individual and variable (size NP x NV x NDRAWS).
# PROBS holds Π_t Pr(choose i at t| b_r) ∀ i, b_r (size NP x NDRAWS)
# Written in Matlab by Kenneth Train, first version Oct 21, 2015 (correction to file CreateZ.m made on Feb 10, 2020).
# Ported to Julia by Giuseppe Forte Dec, 22, 2020.

function createdrawswtpprobs()

    if subsetS == 1                                                                   # S_n ⊂ S is a subset of R randomly selected values of β_r, 
        BETAS = rand(1:NGridPts, NP, NV, NDRAWS);                                     # with all elements of S having the same probability of being selected. (Train 2016, p. 42)
    else
        global NDRAWS = NGridPts * floor(Int64, NDRAWS / NGridPts);                   # Ensure NDRAWS is an exact mulyiple of NGridPts
        BETAS = zeros(NP, NV, NDRAWS);                                                # Fix S_n = S, varying support points across i.
        for j = 1:NP, i = 1:NV, k = 0:(floor(Int64, NDRAWS / NGridPts) - 1)
            BETAS[j, i, (k * NGridPts + 1):(k * NGridPts + NGridPts)] = sample(1:NGridPts, NGridPts, replace = false);
        end
    end
    
    BETAS = (BETAS .- 1) ./ (NGridPts .- 1);                                          # Now BETAS go from zero to 1 inclusive
    for r = 1:NV
        BETAS[:, r, :] = COEF[r, 1] .+ (COEF[r, 2] .- COEF[r, 1]) .* BETAS[:, r, :];  # Now BETAS go from lower limit to upper limit for each coefficient
    end;

    # Probabilities conditional on β_r
    v = zeros(Float64, NROWS, NDRAWS);
    for j in unique(Int.(XMAT[:, 1]))
        a = findfirst(x -> x == j, Int.(XMAT[:, 1]))
        b = findlast(x -> x == j,  Int.(XMAT[:, 1]))
        if wtpspace == 1                                                        # WTP space 
            v[a:b, :] = BETAS[j, 1, :]' .* (XMAT[a:b, IDV] * BETAS[j, 2:end, :] .- XMAT[a:b, IDPRICE]); # (x'γ - p) * α to let γ = β / α be the coefficients scaled to WTP. NROWS x NDRAWS.
        else                                                                    # Preference space
            v[a:b, :] = XMAT[a:b, IDV] * BETAS[j, 2:end, :] .- XMAT[a:b, IDPRICE:IDPRICE] * BETAS[j, 1:1, :]; # - αp + x'β, standard discrete choice approach. NROWS x NDRAWS.
        end
    end
    v = exp.(v);
    sparsematrix = sparse(collect(1:NCS) .== XMAT[:, 2]');                      # allowing for unbalanced histories. NCS x NROWS
    denom = sparsematrix * v;                                                   # NCS x NDRAWS
    p = v[XMAT[:, 3] .== 1, :] ./ denom;                                        # NCS x NDRAWS
    sparsematrix= sparse(collect(1:NP) .== XMAT[XMAT[:, 3] .== 1, 1]');         # NP x NCS
    PROBS = sparsematrix * log.(p);                                             # NP x NDRAWS
    
    return exp.(PROBS), BETAS                                                   # exp(Σ_j log p_j) = Π_j p_j, without as much numerical error. NP x NDRAWS
    
end
