   # This script calls check() to check the data and input specifications, transforms the data
   # into a more easily useable form, calls estimation routine, and prints out results.
   # Written by Kenneth Train on July 19, 2006 (latest edits on Sept 24, 2006). 
   # Ported to Julia by Giuseppe Forte on Dec, 10 2020.

function doit()

   # Check the input data and specifications

   println("\n Checking inputs.\n");
   ok = check();
   if ok == 1
      println("\n Inputs have been checked and look fine.\n");
   else
      return;
   end

   # Create Global variables to use in estimation

   println("\n Creating data arrays for run.\n");
   cp = xmat[:, 1]; # person

   nn = zeros(Int64, ncs, 1);
   for n = 1:ncs
      nn[n, 1] = sum(xmat[:, 2] .== n);
   end
   global naltmax = maximum(nn);     # Maximum number of alternatives in any choice situation

   nn = zeros(Int64, np, 1);
   for n = 1:np
      k = xmat[xmat[:, 1] .== n, 2];
      nn[n, 1] = 1 + k[end, 1] - k[1, 1];
   end
   global ncsmax = maximum(nn);      # Maximum number of choice situations faced by any person

   if wantwgt == 1;
      wgt = zeros(1, np);
      for r = 1:np
         wgt[1, r] = mean(xmat[cp .== r, idwgt]);
      end
      wgt = wgt .* (np / sum(wgt));
   else
      wgt = [];
   end

   # Data arrays 
   # All variables are differenced from the chosen alternative
   # Only nonchosen alternatives are included, since V for chosen alt = 0
   # This reduces number of calculations for logit prob and eliminates need to
   # retain dependent variable.

   global x  = zeros(naltmax - 1, ncsmax, nv, np);     # Explanatory variables with random coefficients for all choice situations, for each person 
   global xf = zeros(naltmax - 1, ncsmax, nf, np);     # Explanatory variables with fixed  coefficients for all choice situations, for each person 
   global s  = zeros(naltmax - 1, ncsmax, np);         # Identification of the alternatives in each choice situation, for each person

   for n = 1:np                                 # loop over people
      cs = Int.(xmat[cp .== n, 2]);
      yy = Int.(xmat[cp .== n, 3]);
      if nv > 0
         xx = xmat[cp .== n, idv[:, 1]];
      end
      if nf > 0
         xxf = xmat[cp .== n, idf[:, 1]];
      end
      t1 = cs[1, 1];
      t2 = cs[end, 1];
      for t = t1:t2;                            # loop over choice situations
         k = sum(cs .== t) - 1;                 # One less than number of alts = number of nonchosen alts
         s[1:k, 1 + t - t1, n] = ones(k, 1);
         if nv > 0
            x[1:k, 1 + t - t1, :, n] = xx[(cs .== t) .& (yy .== 0), :] - repeat(xx[(cs .== t) .& (yy .== 1), :], k);
         end
         if nf > 0
            xf[1:k, 1 + t - t1, :, n] = xxf[(cs .== t) .& (yy .== 0), :] - repeat(xxf[(cs .== t) .& (yy .== 1), :], k);
         end
      end
   end

   Random.seed!(seed1);

   # Create draws

   println("\n Creating draws. \n");
   global dr = permutedims(makedraws(), [3, 2, 1]);   # NMEM x NP x NV -> NV x NP x NDRAWS

   if (nv > 0) & (nf > 0)
      param = [F; B[idv[:, 2] .!= 5]; W];
   elseif (nv > 0) & (nf == 0)
      param = [B[idv[:, 2] .!= 5]; W];
   elseif (nv == 0) & (nf > 0);
      param = F;
   else
      println("\n Model has no explanatory variables: IDV and IDF are both empty. \n Program terminated.\n");
      return
   end

   println("\n Start estimation: The negative of the log-likelihood is minimized.\n");

   if usegradient != 1
      global res = optimize(ll, param, NelderMead(),
                   Optim.Options(f_calls_limit = 10000, x_tol = paramtol, f_tol = lltol, iterations = maxiters, store_trace = true, extended_trace = true))
   else 
      global res = optimize(Optim.only_fg!(llg!), param, BFGS(),
                   Optim.Options(g_tol = 1e-12, f_calls_limit = 10000, x_tol = paramtol, f_tol = lltol, iterations = maxiters, store_trace = true, extended_trace = true))
   end

   global θhat = Optim.minimizer(res);
   global converged = Optim.converged(res);

   global trace = res.trace

   if usegradient != 1
      println("\n Analytical gradient not provided, Hessian and Score not calculated.\n");
      global invHessian = NaN .* zeros(length(param), length(param));
      global gradient   = NaN .* zeros(length(param))
      global θse = NaN .* zeros(length(param));

      t_stats = NaN .* zeros(length(param));
   else 
      global invHessian = res.trace[end].metadata["~inv(H)"];
      global gradient   = res.trace[end].metadata["g(x)"];
      global θse = sqrt.(diag(invHessian));

      t_stats = θhat ./ θse;
   end

   println(@sprintf "\n Estimation took approximately %d second(s).\n" res.trace[end].metadata["time"]);
   if Optim.converged(res)
      println("\n Convergence achieved.\n");
   else
      println("\n Convergence not achieved.\n");
      return
   end

   println(@sprintf "\n Value of the log-likelihood function at convergence: %10.5f \n" (- Optim.minimum(res)));

   println(@sprintf "\n The value of grad*inv(hessian)*grad is: %2.10f \n" gradient' * invHessian * gradient);

   if nf > 0
      ϕhat = θhat[1:nf, 1];
      ϕse  = θse[1:nf, 1];
   end;

   if nv > 0
      if sum(idv[:, 2] .== 5) > 0;
         βhat = zeros(nv);
         βse  = zeros(nv)
         βhat[idv[:, 2] .!= 5] = θhat[nf + 1:nf + sum(idv[:, 2] .!= 5)];
         βse[idv[:, 2] .!= 5]  = θse[nf + 1:nf + sum(idv[:, 2] .!= 5)];
         σhat = θhat[nf + sum(idv[:, 2] .!= 5) + 1:end];
         σse  = θse[nf +  sum(idv[:, 2] .!= 5) + 1:end];
      else;
         βhat = θhat[nf + 1:nf + nv];
         βse  = θse[nf + 1:nf + nv];
         σhat = θhat[nf + nv + 1:nf + nv + nv];
         σse  = θse[nf + nv + 1:nf + nv + nv];
      end;
   end

   println("\n RESULTS\n");
   if nf > 0
      println("              FIXED COEFFICIENTS ϕ\n              --------------------\n                Est         SE\n");
      for r = 1:length(names_idf);
         println(@sprintf "%-10s %10.4f %10.4f\n" names_idf[r] ϕhat[r] ϕse[r]);
      end
   end

   if nv > 0;
      println("\n                        RANDOM COEFFICIENTS \n                       β                      σ");
      println("              ------------------   -----------------------\n                 Est        SE         Est        SE");
      for r = 1:length(names_idv);
         println(@sprintf "%-10s %10.4f %10.4f %10.4f %10.4f\n" names_idv[r] βhat[r] βse[r] σhat[r] σse[r]);
      end
   end

   # Create draws of coefficients from βhat and σhat

   c = reshape(trans(βhat, σhat, dr), nv, np*ndraws);

   println(@sprintf "\n Distribution of coefficients in population implied by βhat and σhat, for NP x NDRAWS draws (%d).\n" ndraws * np);
   jj = ["Normal"; "Lognormal"; "Truncated Normal"; "Johnson S_B"; "Normal (μ = 0)"; "Triangular"];
   println("\n                                     Mean      StdDev   Share < 0  Share = 0    Median\n");
   kk = [mean(c, dims = 2) std(c, dims = 2) mean(c .< 0, dims = 2) mean(c .== 0, dims = 2) median(c, dims = 2)];
   for r = 1:length(names_idv);
      println(@sprintf "%-10s %-20s %10.4f %10.4f %10.4f %10.4f %10.4f\n" names_idv[r] jj[idv[r, 2]] kk[r, 1] kk[r, 2] kk[r, 3] kk[r, 4] kk[r, 5]);
   end

   print("\n ESTIMATED PARAMETERS AND FULL COVARIANCE MATRIX\n");
   println("\n The estimated values of the parameters are:");
   show(stdout, "text/plain", round.(θhat', digits = 5))
   println("\n The covariance matrix for these parameters is:");
   show(stdout, "text/plain", round.(invHessian, digits = 5))

   println("\n\n You can access the estimated parameters as variable θhat, the gradient of the negative of the log-likelihood function as variable gradient,")
   println("the inverse of the hessian (calculated by the BFGS updating procedure) as variable invHessian, and the trace as variable trace.\n");

end

    

