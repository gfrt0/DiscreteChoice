# Calculates log-likelihood function value for mixed logit model
# Written by Kenneth Train, July 27, 2006 (revised July 31, 2006).
# Ported to Julia by Giuseppe Forte on Dec, 10 2020.
#
# This code is input to Julia's optim() 
#
# Input param is a column vector of parameters, dimension (NF+NV+NV)x1
#     containing the fixed coefficients, the first parameters of the random
#     coefficients, and then the second parameters of the random coefficients
# Output ll is the scalar value of the negative of the simulated log-likelihood 
#     at the input parameters

function llg!(F, G, param)
    
    if nf > 0
        ϕ = param[1:nf];                                                           # fixed parameters ϕ
    else
        ϕ = [];
    end;

    if nv > 0
        if sum(idv[:, 2] .== 5) > 0
            β = zeros(nv);
            β[idv[:, 2] .!= 5] = param[nf + 1:nf + sum(idv[:, 2] .!= 5)];          # location parameters β
            σ = param[nf + sum(idv[:, 2] .!= 5) + 1:end];                          # scale parameters σ
        else
            β = param[nf + 1:nf + nv];
            σ = param[nf+ nv + 1:nf + nv + nv];
        end;
    else
        β = [];
        σ = [];
    end;

    # - start of llgrad2.m

    p = zeros(np);                                                                 # Likelihood of each individual history                                                  
    g = zeros(nf + nv + nv, np);                                                   # Score of each individual history wrt ϕ, β, σ.
    c = trans(β, σ, dr);                                                           # Simulated β = β0 + σ * ν. c is NV x NP x ndraws

    if nf > 0
        vf = sum((@view xf[:, :, j, :]) .* ϕ[j] for j in 1:nf);                    # Σ x_j ϕ_j for variables with fixed coefficients. 
    else
        vf = zeros(naltmax - 1, ncsmax, np);
    end;

    hh = zeros(naltmax - 1, ncsmax, np, ndraws);                                   # Σ x_j ϕ_j is the same irrespective of random draw.
    for j = 1:ndraws 
        hh[:, :, :, j] .= vf
    end;

    if nv > 0
        vvvv = zeros(ndraws, np, naltmax - 1, ncsmax);
        xx = permutedims(x, [3, 4, 1, 2]);                                          # x from [naltmax - 1, ncsmax, nv, np] to [nv, np, naltmax - 1, ncsmax] to match c
        for k = 1:ndraws
            vvvv[k, :, :, :] = sum((@view xx[:, :, :, :]) .* c[:, :, k], dims = 1); # Σ_j x_jβ_j, β = β0 + σ * ν, created for each alternative, choice, person, draw.
        end;
        v = permutedims(vvvv, [3, 4, 2, 1]);
        v = v .+ vf;                                                                # v is (NALTMAX - 1) x NCSMAX x NP x NDRAWS
    else
        v = vf;
    end;

    v = exp.(v);                                                                   # exp(V) = exp(x'β)
    v[isinf.(v)] .= 10^20;                                                         # As precaution when exp(v) is too large for machine
    v = v .* s;                                                                    # Only consider choices that were actually available to the consumer.
    pp = 1 ./ (1.0 .+ sum(v, dims = 1));                                           # Likelihood of individual history: Π_t {exp(V_*t) / Σ_j[exp(V_jt)]} = Π_t { 1 / (1 + Σ_{j ≂̸ *} [exp(V_jt - V_*t)]}

    ppp = dropdims(prod(pp, dims = 2), dims = (1, 2))                              # Π_t ∀ i. 
    p = reshape(sum(ppp, dims = 2) ./ ndraws, (np))                                # Revelt Train (1998, eq.3): (1/R) * Σ_r S_n(β) * Σ_t Σ_j (- L (∂βx / ∂θ)).                                   
    p[isnan.(p)] .= 1;                                                             # Change missing values to 1, as a precaution. 

    if G != nothing   # Calculate gradient

        gg = permutedims(permutedims(v, (3, 4, 2, 1)) .* permutedims(pp, (3, 4, 2, 1)), (4, 3, 1, 2));        # Probs for all nonchosen alts. NALTMAX-1 x NCSMAX x NP x ndraws                
        gg = reshape(gg, (naltmax - 1, ncsmax, 1, np, ndraws))

        # ∂ϕ
        if nf > 0
            grϕ = copy(gg);
            jj = 1
            while jj < nf
                grϕ = cat(grϕ, grϕ, dims = 3)
                jj += 1
            end
            grϕ = - grϕ .* xf;                                                     # (- L (∂βx / ∂ϕ)), where ∂xβ / ∂ϕ = x
            gϕ  = reshape(sum(sum(grϕ, dims = 1), dims = 2), (nf, np, ndraws));    # NF x NP x ndraws. Revelt Train (1998, eq.3): Σ_t Σ_j (- L (∂βx / ∂ϕ))
        else
            gϕ  = [];
        end;
        
        # ∂β, ∂σ, etc. 
        if nv > 0
            ggg = zeros(naltmax - 1, ncsmax, nv, np, ndraws);
            for k = 1:nv
                ggg[:, :, k, :, :] = - gg[:, :, 1, :, :] .* x[:, :, k, :]          # Revelt Train (1998, eq.3): (- L (∂βx / ∂θ)) =  (- L x(∂β / ∂θ)).
            end 
            grβ, grσ = der(β, σ, dr);                                              # ∂β / ∂β0, ∂β / ∂σ ... 
            grβ = permutedims(ggg, (3, 4, 5, 1, 2)) .* grβ;
            grσ = permutedims(ggg, (3, 4, 5, 1, 2)) .* grσ;
            gβ  = reshape(sum(sum(grβ, dims = 4), dims = 5), (nv, np, ndraws));    # NV x NP x ndraws. Revelt Train (1998, eq.3): Σ_t Σ_j
            gσ  = reshape(sum(sum(grσ, dims = 4), dims = 5), (nv, np, ndraws));    # NV x NP x ndraws. Revelt Train (1998, eq.3): Σ_t Σ_j
        else
            gβ = [];
            gσ = [];
        end;
        
        # Gradient
        gr = [gϕ; gβ; gσ];
        gr = permutedims(gr, (2, 3, 1)) .* ppp;                                      # Revelt Train (1998, eq.3): S_n(β) * Σ_t Σ_j (- L (∂βx / ∂θ)).
        g = permutedims(dropdims(sum(gr, dims = 2), dims = 2), (2, 1)) ./ ndraws;  # Revelt Train (1998, eq.3): (1/R) * Σ_r S_n(β) * Σ_t Σ_j (- L (∂βx / ∂θ)).                         
        g = g ./ p';                                                               # Revelt Train (1998, eq.3): (1/SP) * (1/R) * Σ_r S_n(β) * Σ_t Σ_j (- L (∂βx / ∂θ)).
        
        if wantwgt == 0
            g  = - sum(g, dims = 2);
        else
            g  = - sum(repeat(wgt, inner = [size(g, 1), 1]) .* g, dims = 2);
        end

        if (nv > 0) & (sum(idv[:, 2] .== 5) > 0)                                   # μ = 0 error components
          z = [ones(nf); idv[:, 2] .!= 5; ones(nv)];
          g = g[z .== 1];
        end
        
        G .= vec(g)
    end

    if F != nothing # Calculate likelihood
        
        if wantwgt == 0
            ll = - sum(log.(p));
        else
            ll = - sum(wgt .* log(p));
        end

        return ll
    end
    # - end of llgrad2.m
end

#########################################################################################

function ll(param)

    if nf > 0
        ϕ = param[1:nf];                                                           # fixed parameters ϕ
    else
        ϕ = [];
    end;

    if nv > 0
        if sum(idv[:, 2] .== 5) > 0
            β = zeros(nv);
            β[idv[:, 2] .!= 5] = param[nf + 1:nf + sum(idv[:, 2] .!= 5)];          # location parameters β
            σ = param[nf + sum(idv[:, 2] .!= 5) + 1:end];                          # scale parameters σ
        else
            β = param[nf + 1:nf + nv];
            σ = param[nf+ nv + 1:nf + nv + nv];
        end;
    else
        β = [];
        σ = [];
    end;

    # - start of llgrad2.m

    p = zeros(np);                                                                 # Likelihood of each individual history                                                  
    c = trans(β, σ, dr);                                                           # Simulated β = β0 + σ * ν. c is NV x NP x ndraws

    if nf > 0
        vf = sum((@view xf[:, :, j, :]) .* ϕ[j] for j in 1:nf);                    # Σ x_j ϕ_j for variables with fixed coefficients. 
    else
        vf = zeros(naltmax - 1, ncsmax, np);
    end;

    hh = zeros(naltmax - 1, ncsmax, np, ndraws);                                   # Σ x_j ϕ_j is the same irrespective of random draw.
    for j = 1:ndraws 
        hh[:, :, :, j] .= vf
    end;

    if nv > 0
        vvvv = zeros(ndraws, np, naltmax - 1, ncsmax);
        xx = permutedims(x, [3, 4, 1, 2]);                                          # x from [naltmax - 1, ncsmax, nv, np] to [nv, np, naltmax - 1, ncsmax] to match c
        for k = 1:ndraws
            vvvv[k, :, :, :] = sum((@view xx[:, :, :, :]) .* c[:, :, k], dims = 1); # Σ_j x_jβ_j, β = β0 + σ * ν, created for each alternative, choice, person, draw.
        end;
        v = permutedims(vvvv, [3, 4, 2, 1]);
        v = v .+ vf;                                                                # v is (NALTMAX - 1) x NCSMAX x NP x NDRAWS
    else
        v = vf;
    end;

    v = exp.(v);                                                                   # exp(V) = exp(x'β)
    v[isinf.(v)] .= 10^20;                                                         # As precaution when exp(v) is too large for machine
    v = v .* s;                                                                    # Only consider choices that were actually available to the consumer.
    pp = 1 ./ (1.0 .+ sum(v, dims = 1));                                           # Likelihood of individual history: Π_t {exp(V_*t) / Σ_j[exp(V_jt)]} = Π_t { 1 / (1 + Σ_{j ≂̸ *} [exp(V_jt - V-*t)]}

    ppp = dropdims(prod(pp, dims = 2), dims = (1, 2))
    p = reshape(sum(ppp, dims = 2) ./ ndraws, (np))                                # Revelt Train (1998, eq.3): (1/R) * Σ_r S_n(β) * Σ_t Σ_j (- L (∂βx / ∂θ)).                                   
    p[isnan.(p)] .= 1;                                                             # Change missing values to 1, as a precaution. 
        
    if wantwgt == 0
        ll = - sum(log.(p));
    else
        ll = - sum(wgt .* log(p));
    end

    return ll
    # - end of llgrad2.m
end
