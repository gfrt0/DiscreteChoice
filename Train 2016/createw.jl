function createw(θ)

    # Z: NPxNDRAWSxNZ array of variables that explain distribution of coefficients
    # NZ: scalar number of variables in Z
    # NDRAWS: scalar number of draws in PROBS

    # Input θ is NZ x 1 vector of coefficients

    w = dropdims(sum(permutedims(Z, (3, 1, 2)) .* θ, dims = 1), dims = 1); # NP x NDRAWS
    w[w .> 500] .= 500;                                                    # As a precaution to prevent machine ∞ when exponentiating
    w = exp.(w);
    w = w ./ sum(w, dims = 2);                                             # NP x NDRAWS

    return w
end