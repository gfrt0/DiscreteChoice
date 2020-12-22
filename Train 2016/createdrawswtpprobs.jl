# Take draws from the parameter space. 
# These draws are held in BETAS which is size NP x NV x NDRAWS
# Then calculate the probability of each person's sequence of choices at each draw for that person. 
# These probabilities are held in PROBS which is size NP x NDRAWS
# Written in Matlab by Kenneth Train, first version Oct 21, 2015 (correction to file CreateZ.m made on Feb 10, 2020).
# Ported to Julia by Giuseppe Forte Dec, 22, 2020.

function createdrawswtpprobs()

    BETAS = rand(1:NGridPts, NP, NV, NDRAWS); # BETAS go from 1 to NGridPts - uniformly distributed.
    BETAS = (BETAS .- 1) ./ (NGridPts .- 1);    # Now BETAS go from zero to 1 inclusive
    for r = 1:NV
    BETAS[:, r, :] = COEF[r, 1] .+ (COEF[r, 2] .- COEF[r, 1]) .* BETAS[:, r, :];  # Now BETAS go from lower limit to upper limit for each coefficient
    end;

    # Add up the non-price variables times their WTP
    v = zeros(Float64, NROWS, NDRAWS);
    for j in unique(Int.(XMAT[:, 1]))
        a = findfirst(x -> x == j, Int.(XMAT[:, 1]))
        b = findlast(x -> x == j,  Int.(XMAT[:, 1]))
        v[a:b, :] = BETAS[j, 1, :]' .* (XMAT[a:b, IDV] * BETAS[j, 2:end, :] .- XMAT[a:b, IDPRICE]); # (x'γ - p) * α to let γ = β / α be the coefficients scaled to WTP. NROWS x NDRAWS.
    end
    v = exp.(v);
    sparsematrix = sparse(collect(1:NCS) .== XMAT[:, 2]');                      # allowing for unbalanced histories. NCS x NROWS
    denom = sparsematrix * v;                                                   # NCS x NDRAWS
    p = v[XMAT[:, 3] .== 1, :] ./ denom;                                        # NCS x NDRAWS
    sparsematrix= sparse(collect(1:NP) .== XMAT[XMAT[:, 3] .== 1, 1]');         # NP x NCS
    PROBS = sparsematrix * log.(p);                                             # NP x NDRAWS
    
    return exp.(PROBS), BETAS                                                   # exp(Σ_j log p_j) = Π_j p_j, without as much numerical error. NP x NDRAWS
    
end
