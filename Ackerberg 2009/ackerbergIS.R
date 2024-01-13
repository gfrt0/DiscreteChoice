### Data setup 

N <- 5000 # number of simulated individuals
J <- 3    # number of inside choices
R <- 100  # simulated shocks per individual

theta <- c(- 1, .5)

alpha_i <- rnorm(N, mean = theta[1], sd = theta[2])

p_ij <- matrix(rexp(J * N), ncol = J) # prices vary at the individual level.

v_ij <- matrix(NA, nrow = N, ncol = J + 1)
s_ij <- data.frame(matrix(NA, nrow = N, ncol = J + 1))

# choice probabilities under the model 
v_ij <- alpha_i * p_ij
s_ij[, 1:J] <- exp(v_ij)
s_ij[, J + 1] <- 1
s_ij <- s_ij / rowSums(s_ij[, 1:(J + 1)])

# observed choices
y_ij <- matrix(0, nrow = N, ncol = J + 1)
for (i in 1:N) {
    y_ij[i, sample(1:(J + 1), 1, prob = s_ij[i, 1:(J + 1)])] <- 1 
}

### maximum simulated likelihood 

id  <- rep(1:N, R)
o_ir <- matrix(rnorm(N * R, mean = 0, sd = 1), nrow = N * R)
p_irj <- rep(1, R) %x% p_ij

s_ij_sim <- function(the) {
  alpha_ir <- the[1] + o_ir * the[2]
  s_irj <- cbind( exp( alpha_ir %*% matrix(1, ncol = J) * p_irj ), 1 )
  s_irj <- s_irj / rowSums(s_irj) 
  s_ij  <- aggregate(s_irj, list(id), mean)
  return( s_ij[, 2:(J + 2)])
}

loglikelihood_sl <- function(the) {
  s_sim <- s_ij_sim(the)
  return( - sum( log( rowSums( y_ij * s_sim ) ) ) )
}

msl <- optim( c(0, 1), loglikelihood_sl )
msl 

### importance sampling with normal proposal density

id  <- rep(1:N, R)
p_irj <- rep(1, R) %x% p_ij

importancesampling <- function(proposal_theta) {
  
  o_ir <- matrix(rnorm(N * R, mean = proposal_theta[1], sd = proposal_theta[2]), nrow = N * R)  # draws from proposal density
  
  s_irj_is <- function(theta) {
    s_irj <- cbind( exp( o_ir %*% matrix(1, ncol = J) * p_irj ), 1 )
    s_irj <- s_irj / rowSums(s_irj) 
    g     <- dnorm(o_ir, mean = proposal_theta[1], sd = proposal_theta[2]) # density under proposal density
    return( list(s = s_irj, g = g, u = o_ir) )
  }
  
  s_is <- s_irj_is()
  
  loglikelihood_is <- function(theta) {
    pu <- dnorm(s_is$u, mean = theta[1], sd = theta[2])
    is_factor <- pu / s_is$g
    
    s_irj_is <- s_is$s * matrix(is_factor, length(is_factor), J + 1)
    s_ij_is <- aggregate(s_irj_is, list(id), mean)[, 2:(J+2)]
    
    return( - sum( log( rowSums( y_ij * s_ij_is ) ) ) )
  }
  
  is <- optim( c(0, 1), loglikelihood_is )
  
  return( list( proposal = proposal_theta, theta0 = theta, thetahat = is$par, 
                ll0 = - loglikelihood_is(theta), llhat = - is$value ) )
}

proposal_thetas <- rbind( cbind( seq(.5, 1.5, length.out = 10) * theta[1], theta[2]), 
                          cbind( theta[1], seq(.5, 1.5, length.out = 10) * theta[2]) )

is_pdensity <- lapply(1:nrow(proposal_thetas), function(x) importancesampling(proposal_thetas[x, ]))

thetahat <- as.data.frame(t(sapply(1:length(is_pdensity), function(x) is_pdensity[[x]]$thetahat)))
thetahat