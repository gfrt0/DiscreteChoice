# Calculate logit probability for chosen alternatives for each person
# with multiple choice situations for each person and multiple people,
# using globals for all inputs except coefficients.
# Written by Kenneth Train, first version on July 14, 2006; latest edits on Sept 24, 2006.
#
# Simulated Mixed Logit Probability and gradient
# Logit probability is Π_t {exp(V_*t) / Σ_j[exp(V_jt)]} = Π_t { 1 / (1 + Σ_{j ≂̸ *} [exp(V_jt - V-*t)]}
# where * denotes chosen alternative, j is alternative, t is choice situation.
#
# Using differences from chosen alternative reduces computation time (with
# one less alternative), eliminates need to retain and use the dependent
# variable, and avoids numerical problems when exp(V) is very large or
# small, since prob is 1/(1+k) which can be evaluated for any k, including
# infinite k and infinitesimally small k. In contrast, e(V)/sum e(V) can
# result in numerical "divide by zero" if denominator is sufficiently small
# and NaN if numerator and denominator are both numerically zero.
# 
# Input f contains the fixed coefficients, and has dimension NF x 1.
# Input c contains the random coefficients for each person, and has dimensions NV x NP.
# Either input can be an empty matrix. 
# Output p contains the logit probabilities, which is a row vector of dimension 1 x NP.
# Output g contains the gradients of log(p), which is a matrix (NF+NV+NV) x NP
# Code assumes that all GLOBALS already exist.

function llgrad2(f, b, w) 

    p = zeros(np);                                                          
    g = zeros(nf + nv + nv, np);

    c = trans(b, w, dr);                                                       # Gets simulated β = β0 + σ * ν. c is NV x NP x ndraws
    v = zeros(ndraws, naltmax - 1, ncsmax, np);
    
    if nf > 0
        ff = reshape(f, (1, 1, nf, 1));
        ff = repeat(ff, inner = [naltmax - 1, ncsmax, 1, np]);                 # fixed coefficients ϕ for all alternatives, choices, people
        vf = reshape(sum(xf .* ff, dims = 3), (naltmax - 1, ncsmax, np));      # Σ x_j ϕ_j for variables x_j with fixed coefficients. 
    else
        vf = zeros(naltmax - 1, ncsmax, np);
    end
    vf = repeat(vf, inner = [1, 1, 1, ndraws]);                                # Σ x_j ϕ_j is the same irrespective of random draw.

    if nv > 0
        cc = reshape(c, (1, 1, nv, np, ndraws));
        cc = repeat(cc, inner = [naltmax - 1, ncsmax, 1, 1, 1]);               # β = β0 + σ * ν created for each alternative and choice
        v  = repeat(x,  inner = [1, 1, 1, 1, ndraws]) .* cc;                   # x_jβ_j created for each alternative, choice, variable, person, draw.
        v  = reshape(sum(v, dims = 3), naltmax - 1, ncsmax, np, ndraws);       # Σ_j x_jβ_j created for each alternative, choice, person, draw.
        v  = v .+ vf;                                                          # v is (NALTMAX - 1) x NCSMAX x NP x NDRAWS
    else
        v  = vf;
    end

    v = exp.(v);                                                               # exp(V) = exp(x'β)
    v[isinf.(v)] .= 10^20;                                                     # As precaution when exp(v) is too large for machine
    v = v .* repeat(s, inner = [1, 1, 1, ndraws]);                             # Only consider choices that were actually available to the consumer.
    pp = 1 ./ (1.0 .+ sum(v, dims = 1));                                       # Likelihood of individual history: Π_t {exp(V_*t) / Σ_j[exp(V_jt)]} = Π_t { 1 / (1 + Σ_{j ≂̸ *} [exp(V_jt - V-*t)]}

    # Calculate gradient

    gg = v .* repeat(pp, inner = [naltmax - 1, 1, 1, 1]);                      # Probs for all nonchosen alts NALTMAX-1 x NCSMAX x NP x ndraws 
    gg = reshape(gg, (naltmax - 1, ncsmax, 1, np, ndraws));                    

    # ∂ϕ
    if nf > 0
        grf = - repeat(gg, inner = [1, 1, nf, 1, 1]) .* repeat(xf, inner = [1, 1, 1, 1, ndraws]); # ∂xβ / ∂ϕ = x
        grf = reshape(sum(sum(grf, dims = 1), dims = 2), (nf, np, ndraws));    # NF x NP x ndraws. Revelt Train (1998, eq.3): Σ_t Σ_j
    else
        grf = [];
    end
    
    # ∂β, ∂σ, etc. 
    if nv > 0
        gg = - repeat(gg, inner = [1, 1, nv, 1, 1]) .* repeat(x, inner = [1, 1, 1, 1, ndraws]);
        grb, grw = der(b, w, dr);                                              # ∂xβ / ∂β0, ∂xβ / ∂σ ... 
        grb = reshape(grb, (1, 1, nv, np, ndraws));
        grw = reshape(grw, (1, 1, nv, np, ndraws));
        grb = gg .* repeat(grb, inner = [naltmax - 1, ncsmax, 1, 1, 1]);
        grw = gg .* repeat(grw, inner = [naltmax - 1, ncsmax, 1, 1, 1]);
        grb = reshape(sum(sum(grb, dims = 1), dims = 2), nv, np, ndraws);      # NV x NP x ndraws. Revelt Train (1998, eq.3): Σ_t Σ_j
        grw = reshape(sum(sum(grw, dims = 1), dims = 2), nv, np, ndraws);      # NV x NP x ndraws. Revelt Train (1998, eq.3): Σ_t Σ_j
    else
        grb = [];
        grw = [];
    end
    
    # Back to prob
    pp = reshape(pp, (ncsmax, np, ndraws));                                     # pp is now ncsmax x np x ndraws
    pp = prod(pp, dims = 1);                                                    # pp is 1 x np x ndraws

    # Gradient
    gr = [grf; grb; grw];
    gr = gr .* repeat(pp, inner = [nf + nv + nv, 1, 1]);                        # Revelt Train (1998, eq.3): S_n(β) * Σ_t Σ_j (- L (∂βx / ∂θ)).
    g = g + sum(gr, dims = 3);                                                  

    # Back to prob (1 x NP)
    pp = sum(pp, dims = 3);                                                     # Revelt Train (1998, eq.3): Σ_r S_n(β) * Σ_t Σ_j (- L (∂βx / ∂θ)).
    p  = p + reshape(pp, (ndraws));

    # Probabilities
    p = p ./ ndraws;
    p[isnan.(p)] .= 1; # Change missing values to 1, as a precaution. 
    # Gradient
    g = reshape(g ./ ndraws, (nf + nv + nv, ndraws));                            # Revelt Train (1998, eq.3): (1/R) * Σ_r S_n(β) * Σ_t Σ_j (- L (∂βx / ∂θ)).
    g = g ./ repeat(p', inner = [nf + nv + nv, 1]);                              # Revelt Train (1998, eq.3): (1/SP) * (1/R) * Σ_r S_n(β) * Σ_t Σ_j (- L (∂βx / ∂θ)).

    return p, g
end