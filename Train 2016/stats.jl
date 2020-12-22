# Function to create summary statistics of distribution approximations.
# Written in Matlab by Kenneth Train, first version Oct 21, 2015 (correction to file CreateZ.m made on Feb 10, 2020).
# Ported to Julia by Giuseppe Forte Dec, 22, 2020.

function stats(θ,NBins)

    w = createw(θ); # NP x NDRAWS
    mn = zeros(NV, 1);
    v  = zeros(NV, 1);
    cc = zeros(NV, NV);
    freqbins = zeros(NV, NBins);
    midbins  = zeros(NV, NBins);
    wtsum = sum(w);
    ww = reshape(w, NP * NDRAWS, 1);
    for r = 1:NV
        thiscoef = BETAS[:, r, :];                                          # NP x NDRAWS
        mn[r, 1] = sum(thiscoef .* w) ./ wtsum;
        v[r, 1]  = sum(((thiscoef .- mn[r, 1]) .^ 2) .* w) ./ wtsum;
        cc[r, r] = v[r, 1];
        if r > 1 & CrossCorr == 1;
            for thisc = 2:NV;
                thiscoef2 = BETAS[:, thisc, :];                             # NP x NDRAWS
                thismn2 = sum(thiscoef2 .* w) ./ wtsum;
                cc[r, thisc] = sum(((thiscoef .- mn[r, 1]) .* (thiscoef2 .- thismn2)) .* w) ./ wtsum;
                cc[thisc, r] = cc[r, thisc];
            end
        end
        histw, midpt = histwgt(reshape(thiscoef, NP * NDRAWS, 1), ww, COEF[r, 1], COEF[r, 2], NBins); 
        freqbins[r, :] = histw ./ sum(histw);
        midbins[r, :]  = midpt;
    end
    stdv = sqrt.(v);

    return mn, stdv, cc, freqbins, midbins
end