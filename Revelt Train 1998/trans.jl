# Transform normally distributed terms into coefficients
# Written by Kenneth Train, July 14, 2006 (revised on July 28, 2006).
# Ported to Julia by Giuseppe Forte on Dec, 10 2020.

# Input b has dimension NV x 1. 
# Input w has dimension NV x 1.
# Input dr are the draws and have dimension NV x NP x NDRAWS. 
# Output c has dimension NV x NP x NDRAWS.
# Uses IDV to determine transformations of draws in dr.

function trans(b, w, dr)
    if (nv > 0)
        c = repeat(b, inner = [1 np ndraws]) .+ repeat(w, inner = [1 np ndraws]) .* dr;
        c[idv[:, 2] .== 2, :, :] = exp.(c[idv[:, 2] .== 2, :, :]);                                              # lognormal
        c[idv[:, 2] .== 3, :, :] = c[idv[:, 2] .== 3, :, :] .* (c[idv[:, 2] .== 3, :, :] .> 0);                 # truncated normal 
        c[idv[:, 2] .== 4, :, :] = exp.(c[idv[:, 2] .== 4, :, :]) ./ (1.0 .+ exp.(c[idv[:, 2] .== 4, :, :]));   # Johnson S_B distribution 
    else
        c=[];
    end
    c
end 