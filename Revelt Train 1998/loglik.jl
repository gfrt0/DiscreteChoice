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

function loglik(param)

  if nf > 0
    f = param[1:nf];
  else
    f=[];
  end

  if nv > 0
    if sum(idv[:, 2] .== 5) > 0
      b = zeros(nv);
      b[idv[:, 2] .!= 5] = param[nf + 1:nf + sum(idv[:, 2] .!= 5)];
      w = param[nf + sum(idv[:, 2] .!= 5) + 1:end];
    else
      b = param[nf + 1:nf + nv];
      w = param[nf+ nv + 1:nf + nv + nv];
    end;
  else
    b=[];
    w=[];
  end

  p, g = llgrad2(f, b, w); 

  if wantwgt == 0
      ll = - sum(log.(p),dims = 2);
      g  = - sum(g, dims = 2);
  else
      ll = - sum(wgt .* log(p));
      g  = - sum(repeat(wgt, inner = [size(g, 1), 1]) .* g, dims = 2);
  end

  if (nv > 0) & (sum(idv[:, 2] .== 5) > 0)  # Zero mean error components
    z = [ones(nf); idv[:, 2] .!= 5; ones(nv)];
    g = g[z .== 1, 1];
  end

  return ll, g
end
