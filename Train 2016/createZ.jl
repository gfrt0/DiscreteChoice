#Create the z variables that describe density over coefficient space
#The variables are given in an array that is NPxNDRAWSxNZ
# where NZ is the number of these variables

#ZTYPE determines the kind of variables to create
# 1=polynomial, 2= step function, 3=spline

function createZ()

    # Polynomial
    if ZTYPE == 1;  
        NZ = NV * PolyOrder;
        if CrossCorr == 0
            Z = zeros(NP, NDRAWS, NZ);
        else 
            NewNZ = Int(NZ + ((NV - 1) * (NV - 2) ./ 2));                                                # NZ + upper triangle of correlation matrix
            Z = zeros(NP, NDRAWS, NewNZ);
        end 
        κ = 0;
        for r = 1:NV;
            for o = 1:PolyOrder;
                κ += 1;
                yy = legendre(LegendreSphereNorm(), o, o, - 1 .+ 2 .* (BETAS[:, r, :] .- COEF[r, 1]) ./ (COEF[r, 2] .- COEF[r, 1]));
                Z[:, :, κ] = yy;
            end;
        end;
        if CrossCorr == 1;
            κ = 0;
            for j = (PolyOrder + 1):PolyOrder:(NV * PolyOrder);                                          # First coef is for price/scale coef, which is not correlated
                for i = (PolyOrder + 1):PolyOrder:(j - 1);  
                    κ += 1;
                    Z[:, :, NZ + κ] = Z[:, :, j] .* Z[:, :, i];
                end;
            end
        end;
    # Step functions
    elseif ZTYPE == 2;  
        NZ = (NLevels - 1) * NV;                                                                         # One parameter for each level minus one by normalization
        if CrossCorr == 0
            Z = zeros(NP, NDRAWS, NZ);
        else 
            NewNZ = Int(NZ + 2 .* (NV - 1) + ((NV - 2) * (NV - 1) ./ 2));                                # NZ + 2(NV - 1) + upper triangle of correlation matrix
            Z = zeros(NP, NDRAWS, NewNZ);
        end;
        CUTS = zeros(NV, NLevels + 1);
        CUTS[:, 1] = COEF[:, 1];
        cutsize = (COEF[:, 2] - COEF[:, 1]) ./ NLevels;                                                  # Evenly spaced steps
        for k = 1:NLevels;
            CUTS[:, k + 1] = CUTS[:, 1] + k .* cutsize;
        end
        κ = 0;
        for r = 1:NV;
            for k = 1:NLevels - 1;
                κ += 1;
                Z[:, :, κ] = (BETAS[:, r, :] .>= CUTS[r, k]) .& (BETAS[:, r, :] .< CUTS[r, k + 1]);
            end;
        end;
        if CrossCorr == 1;
            for r = 1:NP
                Z[r, :, (NZ + 1):(NZ + 2 * (NV - 1))] = [BETAS[r, 2:end, :]; BETAS[r, 2:end, :].^2]';    # Do not include price/scale coef
            end
            NNZZ = NZ + 2 * (NV - 1);
            κ = 0;
            for j = 1:(NV - 1);                                                                          # First coef is for price/scale coef, which is not correlated
                for i = 1:(j - 1);
                    κ += 1;
                    Z[:, :, NNZZ + κ] = Z[:, :, NZ + j] .* Z[:, :, NZ + i];
                end;
            end;
        end;
    # Linear spline
    elseif ZTYPE == 3; 
        NZ = (NKnots + 1) * NV;  # Number of parameters per variable is [1 startpoint] + [1 endpoint] + NKnots - [1 for normalization] = NKnots + 1. Height at endpoint is normalized to zero.
        if CrossCorr == 0
            Z = zeros(NP, NDRAWS, NZ);
        else 
            NewNZ = Int(NZ + 2 .* (NV - 1) + ((NV - 2) * (NV - 1) ./ 2));                                # NZ + 2(NV - 1) + upper triangle of correlation matrix
            Z = zeros(NP, NDRAWS, NewNZ);
        end;
        CUTS = zeros(NV, NKnots + 2);                                                                    # Knots and endpoints (evenly spaced)
        CUTS[:, 1] = COEF[:, 1];
        cutsize = (COEF[:, 2] - COEF[:, 1]) ./ (NKnots + 1); 
        for k = 1:(NKnots + 1);
            CUTS[:, k + 1] = CUTS[:, 1] + k .* cutsize;
        end
        κ = 0;
        for r = 1:NV;
            for χ = 1:(NKnots + 1);
                κ += 1;
                inseg = (BETAS[:, r, :] .>= CUTS[r, χ]) .& (BETAS[:, r, :] .< CUTS[r, χ + 1]);           # NP x 1 x NDRAWS
                m = (BETAS[:, r, :] .- CUTS[r, χ]) ./ (CUTS[r, χ + 1] - CUTS[r, χ]);                     # NP x 1 x NDRAWS
                Z[:, :, κ] = Z[:, :, κ] .+ (1 .- m) .* inseg;
                if χ < (NKnots + 1)
                    Z[:, :, κ + 1] = Z[:, :, κ + 1] + m .* inseg;
                end;
            end;
        end;
        if CrossCorr == 1;
            for r = 1:NP
                Z[r, :, (NZ + 1):(NZ + 2 * (NV - 1))] = [BETAS[r, 2:end, :] ; BETAS[r, 2:end, :].^2]';   # Do not include price/scale coef
            end
            NNZZ = NZ + 2 .* (NV - 1);
            κ = 0;
            for j = 1:(NV - 1);                                                                          # First coef is for price/scale coef, which is not correlated
                for i = 1:(j - 1);
                    κ += 1;
                    Z[:, :, NNZZ + κ] = Z[:, :, NZ + j] .* Z[:, :, NZ + i];
                end;
            end;
        end;
    end;
    
    return Z

end

