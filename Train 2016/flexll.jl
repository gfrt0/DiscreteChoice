# Likelihood function
# Written in Matlab by Kenneth Train, first version Oct 21, 2015 (correction to file CreateZ.m made on Feb 10, 2020).
# Ported to Julia by Giuseppe Forte Dec, 22, 2020.

function flexll!(F, G, θ)

    # PROBS: NP x NDRAWS matrix of probabilities for each person at each sampled coefficients
    # Z:     NP x NDRAWS x NZ array of variables that explain distribution of coefficients
    # NZ: scalar number of variables in Z
    # NDRAWS: scalar number of draws in PROBS

    # Input θ is NZ x 1 vector of coefficients

    w = dropdims(sum(permutedims(Z, (3, 1, 2)) .* θ, dims = 1), dims = 1); # NP x NDRAWS
    w[w .< -500] .= -500;                                                  # As precaution against extreme parameters
    w[w .> 500]  .=  500;
    w = exp.(w);
    w = w ./ sum(w, dims = 2);                                             # NP x NDRAWS
    logit_w = w .* PROBS;                                                  # NP x NDRAWS
    mix_probs = sum(logit_w, dims = 2);                                    # NP x 1
    
    if G != nothing 
        h = logit_w ./ mix_probs;                                          # NP x NDRAWS
        g  = (h - w) .* Z;                                                 # Train 2016, eq. 10: [h(β_r | α) - w(β_r | α)] × z(β_r). NP x NDRAWS x NZ
        G .= vec(- sum(sum(g, dims = 1), dims = 2));                       # 1 x NZ, to minimise
    end;

    if F != nothing
        ll = - sum(log.(mix_probs), dims = 1);                             # 1 x 1, to minimize
        return ll[1]
    end;
    
end
