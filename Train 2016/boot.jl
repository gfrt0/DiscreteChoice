# Bootstraps NReps the model run in doestimation.jl and provides summary table.
# Written in Matlab by Kenneth Train, first version Oct 21, 2015 (correction to file CreateZ.m made on Feb 10, 2020).
# Ported to Julia by Giuseppe Forte Dec, 22, 2020.

function boot(XMAT)

    ParamSe = zeros(NZ, 1);                          # Standard error of estimated coefficients of Z variables
    MeanSE  = zeros(NV, 1);                          # Standard error of Mean
    StdSE   = zeros(NV, 1);                          #  " of Std Devs
    CorrSE  = zeros(NV, NV);                         # " of Correlation matrix
    FreqSE  = zeros(NV, NBins);                      # " of Freq in each bin
    MidSE   = zeros(NV, NBins);                      # " of Midpoint for each bin, (std err should be zero)

    paramhold = zeros(NZ, NReps);
    mnhold   = zeros(NV, NReps);
    stdhold  = zeros(NV, NReps);
    corrhold = zeros(NV - 1, NV - 1, NReps);
    freqhold = zeros(NV, NBins, NReps);
    midhold  = zeros(NV, NBins, NReps);

    XMAT_Original = XMAT;
    PID_Original  = XMAT[:, 1];
    CSID_Original = XMAT[:, 2];

    for ξ = 1:NReps
        btsample = rand(1:NP, NP, 1);
        XMAT = reshape([], 0, size(XMAT_Original, 2));
        for r = 1:NP;
            thisp = btsample[r, 1];
            thisx = XMAT_Original[PID_Original .== thisp, :];
            thisx[:, 1] = r .* ones(Int64, size(thisx, 1), 1);
            XMAT = cat(XMAT, thisx, dims = 1);   
        end
        println(@sprintf "Done with optimisation of bootstrap sample %d" ξ)
        XMAT[:, 2] = CSID_Original;
        PROBS, BETAS = createdrawswtpprobs()
        Z = createZ()
        res = optimize(Optim.only_fg!(flexll!), StartB, BFGS(),
        Optim.Options(g_tol = GTOL, f_calls_limit = 10000, x_tol = PARAMTOL, f_tol = LLTOL, iterations = MAXITERS, store_trace = true, extended_trace = true))
        paramhold[:, ξ] = Optim.minimizer(res);
        mnhold[:, ξ], stdhold[:, ξ], cc, freqhold[:, :, ξ], midhold[:, :, ξ] = stats(Optim.minimizer(res), NBins);
        corrhold[:, :, ξ] = Statistics.cov2cor!(cc[2:end, 2:end], sqrt.(diag(cc[2:end, 2:end])));
    end

    global ParamSE = std(paramhold; corrected = true, dims = 2);
    global MeanSE  = std(mnhold; corrected = true, dims = 2);
    global StdSE   = std(stdhold; corrected = true, dims = 2);
    global CorrSE  = dropdims(std(corrhold; corrected = true, dims = 3), dims = 3); 
    global FreqSE  = dropdims(std(freqhold; corrected = true, dims = 3), dims = 3); 
    global MidSE   = dropdims(std(midhold; corrected = true, dims = 3), dims = 3);
    XMAT = XMAT_Original;

    println("\nBootstrapped Utility Parameters:\n");
    println("                    Mean                   Std Dev");
    println("              ------------------   -----------------------");
    println("                 Est     SE            Est         SE");
    for r = 1:length(NAMES);
        println(@sprintf "%-10s %10.4f %10.4f %10.4f %10.4f\n" NAMES[r, 1] MeanEst[r, 1] MeanSE[r, 1] StdEst[r, 1] StdSE[r, 1]);
    end

    println("\nCorrelations of WTPs:");
    show(stdout, "text/plain", round.(Statistics.cov2cor!(CovMatEst[2:end, 2:end], sqrt.(diag(CovMatEst[2:end, 2:end]))), digits = 5))
    println("\nT-stats on Correlations (diagonal has Inf because diagonals of correlation matrix are always 1):");
    show(stdout, "text/plain", round.(Statistics.cov2cor!(CovMatEst[2:end, 2:end], sqrt.(diag(CovMatEst[2:end, 2:end]))) ./ CorrSE, digits = 5))
    println("\nTo obtain histogram for utility parameter k, type: histogram(MidEst[k, :], bins = length(MidEst[k, :]), weights = FreqEst[k, :])\n");

end