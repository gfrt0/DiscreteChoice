# Take derivates of random coefficients wrt b and w
# Written by Kenneth Train, July 28, 2006 (latest edits on Aug 8, 2006).
# Ported to Julia by Giuseppe Forte on Dec, 10 2020.

# Input β has dimension NVx1. 
# Input σ has dimension NVx1.
# Input dr are the draws and have dimension NV x NP x ndraws  
# Output der has dimension NV x NP x ndraws.
# Uses IDV to determine transformations of draws in dr.

function der(β, σ, dr)
 
    dβ = ones(nv, np, ndraws);

    if sum((idv[:, 2] .== 2) + (idv[:, 2] .== 3) + (idv[:, 2] .== 4)) > 0

        c = β .+ σ .* dr;

        dβ[idv[:, 2] .== 2, :, :] = exp.(@view c[idv[:, 2] .== 2, :, :]);                                           # Lognormal derivative

        dβ[idv[:, 2] .== 3, :, :] = (@view c[idv[:, 2] .== 3, :, :]) .> 0;                                          # Truncated normal derivative

        dβ[idv[:, 2] .== 4, :, :] = exp.(@view c[idv[:, 2] .== 4, :, :])./ (1.0 .+ exp.(@view c[idv[:, 2] .== 4, :, :]));
        dβ[idv[:, 2] .== 4, :, :] = (@view dβ[idv[:, 2] .== 4, :, :]) .- ((@view dβ[idv[:, 2] .== 4, :, :]).^2);    # Johnson S_B derivative
        
    end

    dσ = dβ .* dr;
    
    return dβ, dσ
end