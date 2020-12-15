# Calculate draws of standardized random terms for random coefficients and save if specified
# Written by Kenneth Train, August 6, 2006.
# Ported to Julia by Giuseppe Forte on Dec, 10 2020.

# Array of draws has dimension NDRAWS x NP x NV.

triang(draws) = (sqrt(2 .* draws) - 1) .* (draws .≤ .5) + (1 - sqrt(2 .* (1 - draws))) .* (draws .> .5) 

function makedraws()

    dr = zeros(Float64, ndraws, np, nv);

    if drawtype == 1                                # Random draws
        for j = 1:nv
            if idv[j, 2] != 6                       # Normal draws
                dr[:, :, j] = randn(ndraws, np);
            else                                    # Triangular draws
                draws = rand(ndraws, np);
                dr[:, :, j] = triang(draws);
            end
        end
    end

    if drawtype ∈ 2:3                               # Halton draws
        h = primes(100);
        k = 1;
        while length(h) < nv
            h = primes(k .* 100);
            k += 1;
        end
        h = h[1:nv];
        for j = 1:nv
            hh = h[j];
            draws = [0.0];
            test = 0;
            b = 1;
            while test == 0
                drawsold = draws;
                for m = 1:(hh - 1);
                    dd = m / (hh^b);
                    draws = [draws; drawsold .+ [dd]];
                    test = length(draws) .≥ ((np .* ndraws) + 10);
                    if test == 1
                        break
                    end
                end
                b += 1;    
            end
            draws = draws[11:(10 + np*ndraws), 1];
            if drawtype == 3
                draws = draws .+ rand();                      # Shift: one shift for entire sequence
                draws = draws - floor.(draws);
                draws = reshape(draws, ndraws, np);
                for n = 1:np                                  # Shuffle for each person separately
                    draws[:, n] = draws[sortperm(rand(ndraws)), n];
                end
                draws = reshape(draws, np * ndraws, 1);
            end
            if idv[j, 2] != 6
                draws = - sqrt(2) * erfcinv.(2 * draws);      # Take inverse cum normal
            else
                draws = triang(draws); 
            end
            dr[:, :, j] = reshape(draws, ndraws, np);
        end
    end

    if drawtype == 4                                          # MLHS
        h = 0:(ndraws - 1);
        h = h' / ndraws;
        for j = 1:nv
            for n = 1:np
                draws = h .+ rand() / ndraws;                 # Shift: Different shift for each person
                draws = draws[sortperm(rand(ndraws))];                          # Shuffle
                if idv[j, 2] != 6
                    draws = -sqrt(2) * erfcinv.(2 * draws);    # Take inverse cum normal
                else
                    draws = triang(draws); 
                end
                dr[:, n, j] = draws;
            end
        end
    end
    dr;
end

