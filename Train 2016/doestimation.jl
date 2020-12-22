# Runs the flexible logit model specified in FlexibleWtp.jl
# Written in Matlab by Kenneth Train, first version Oct 21, 2015 (correction to file CreateZ.m made on Feb 10, 2020).
# Ported to Julia by Giuseppe Forte Dec, 22, 2020.

function doestimation()

    println("\nStart estimation.\n");
    println("\nThe negative of the log-likelihood is minimized, which is the same as maximizing the log-likelihood.\n");

    global res = optimize(Optim.only_fg!(flexll!), StartB, BFGS(),
                Optim.Options(g_tol = GTOL, f_calls_limit = 10000, x_tol = PARAMTOL, f_tol = LLTOL, iterations = MAXITERS, store_trace = true, extended_trace = true))

    global θhat = Optim.minimizer(res);
    global converged = Optim.converged(res);
    
    global trace = res.trace
    
    global invHessian = res.trace[end].metadata["~inv(H)"];
    global gradient   = res.trace[end].metadata["g(x)"];
    global θse = sqrt.(diag(invHessian));

    t_stats = θhat ./ θse;
    
    println(@sprintf "\n Estimation took approximately %2.2f minutes (%d seconds).\n" (res.trace[end].metadata["time"] / 60) res.trace[end].metadata["time"]);
    
    if Optim.converged(res)
        println("\n Convergence achieved.\n");
    else
        println("\n Convergence not achieved.\n");
        return
    end

    println("\nCalculating summary statistics for random utility parameters.\n");
    global MeanEst, StdEst, CovMatEst, FreqEst, MidEst = stats(θhat, NBins);

    println("\nEstimated coefficients of Z variables are held in θhat.\n");
    println("\nInverse Hessian at convergence is held in invHessian.\n");

    println("\nMeans, Standard Deviations, and Covariance Matrix of Utility Parameters are held as MeanEst, StdEst, and CovMatEst.");
    println("Share of density in each bin and midpoint for each bin are held in FreqEst and MidEst; each is size NV x NBins.\n");

    println("\nMeans and StdDevs of Random Utility Paramaters\n");
    println("                Mean     Std Dev");
    println("              -------------------");
    for r = 1:length(NAMES);
        println(@sprintf "%-10s %10.4f %10.4f\n" NAMES[r] MeanEst[r, 1] StdEst[r, 1]);
    end
    println("\nCorrelation Matrix for WTPs:");
    show(stdout, "text/plain", round.(Statistics.cov2cor!(CovMatEst[2:end, 2:end], sqrt.(diag(CovMatEst[2:end, 2:end]))), digits = 5))
    
    println("\nTo obtain histogram for utility parameter k, type: histogram(MidEst[k, :], bins = length(MidEst[k, :]), weights = FreqEst[k, :])\n");
    
end