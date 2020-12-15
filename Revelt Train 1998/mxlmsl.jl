# Giuseppe Forte, UCL - December, 10th 2020
# Julia port of Kenneth Train's mixed logit estimated with MSL.

#                                               * NOTES *
# Eliminated the possibility of keeping in memory only a subset of draws - no more NMEM
# Eliminated the possibility of using stored draws. 
# Added a likelihood function that does not compute the gradient (for gradient-free methods) and the possibiity to not pass gradient info.

using DelimitedFiles, Printf, Random, Statistics, Primes, SpecialFunctions, Optim, NLSolversBase, ForwardDiff, LinearAlgebra

include("check.jl")
include("makedraws.jl")
include("trans.jl")
include("der.jl")
include("llg.jl")
include("doit.jl")

# TITLE
print("Mixed logit with 100 MLHS draws: 3 lognormal, 2 normal, 2 fixed.")

# DATA

# Number of people (decision-makers) in dataset
np = 100

# Number of choice situations in dataset. This is the number faced by all the people combined.
ncs = 1484

# Total number of alternatives faced by all people in all choice situations combined.
# This is the number of rows of data in XMAT below.
nrows = 1484 * 3

# Load and/or create XMAT, a matrix that contains the data.
#
# XMAT must contain one row of data for each alternative in each choice situation for each person.
# The rows are grouped by person, and by choice situations faced by each person.
# The number of rows in XMAT must be NROWS, specified above.
# The columns in XMAT are variable that describe the alternative.
#
# The *first* column of XMAT identifies the person who faced this alternative.
# The people must be numbered sequentially from 1 to NP, in ascending order.
# All alternatives for a given person must be grouped together.
# The *second* column of XMAT identifies the choice situation. The choice
# situations must be numbered sequentially from 1 to NCS.
# All alternatives for a given choice situation must be grouped together.
# The *third* column of XMAT identifies the chosen alternatives (1 for
# chosen, 0 for not). One and only one alternative must be chosen for each
# choice situation.
# The remaining columns of XMAT can be any variables.

xmat = readdlm("data.txt");         # The variables are described below

xmat[:, 4:5] = - xmat[:, 4:5];      # To make price and opcost negative so coef can be positive.
xmat[:, 4] = xmat[:, 4] ./ 10000;   # To scale price to be in tens of thousands of dollars.

# To help you keep up with the variables, list the variables in XMAT here.
# Start each line with # so that matlab sees that it is a comment rather than a command.
# NOTES for XMAT for sample run:
# This dataset is for people's choice among vehicles in stated-preference
# experiments. Each person faced up to 15 experiments (some faced fewer
# than 15 because they did not complete all the experiments.) Each
# experiment contained 3 alternatives representing three different vehicles
# whose price and other attributes were described. The person stated which
# of the three vehicle he/she would buy if facing this choice in the real world.
# The variables in XMAT are:
# 1. Person number (1-NP)            MUST BE THIS. DO NOT CHANGE.
# 2. Choice situation number (1-NCS) MUST BE THIS. DO NOT CHANGE.
# 3. Chosen alternative (1/0)        MUST BE THIS. DO NOT CHANGE.
# 4. Negative of Price in tens of thousands of dollars
# 5. Negative of Operating cost in dollars per month
# 6. Range in hundreds of miles (0 if not electric)
# 7. Electric (1/0)
# 8. Gas (1/0)
# 9. Hybrid (1/0)
# 10. High performance (1/0)
# 11. Medium or high performance (1/0)

# MODEL SPECIFICATION

# RANDOM COEFFICIENTS
# List the variables in XMAT that enter the model with random coefficients and
# give the distribution for the coefficient of each variable.
# IDV contains one row for each random coefficient and two columns.
# The *first* column gives the number of a variable in XMAT that has a random coefficient,
# and the *second* column specifies the distribution of the coefficient for that variable.
# The distributions can be
# 1. normal: N(b,w^2) where mean b and standard deviation w are estimated.
# 2. lognormal: coefficient is exp(beta) where beta~N(b,w^2) with b and w estimated
# 3. truncated normal, with the share below zero massed at zero: max(0,beta) where
#                      beta~N(b,w^2) with b and w estimated.
# 4. S_B: exp(beta)/(1+exp(beta))  where beta~N(b,w^2) with b and w estimated.
# 5. normal with zero mean (for error components): N(0,w^2) where w is estimated.
# 6. triangular: b+w*t where t is triangular between -1 and 1 and mean b and spread w are estimated.
# If no random coefficients, put IDV=[];
# Notes:
# The lognormal, truncated normal, and S_B distributions give positive
# coefficients only. If you want a variable to have only negative coefficients,
# create the negative of the variable (in the specification of XMAT above).
# The S_B distribution gives coefficients between 0 and 1. If you want
# coefficients to be between 0 and k, then multiply the variable by k (in the specification
# of XMAT above), since b*k*x for b~(0-1) is the same as b*x for b~(0-k).
# If no random coefficients, put IDV=[];

idv =[4 2;
      5 2;
      6 2;
      7 1;
      9 1];

nv = size(idv, 1); #Number of random coefficients. Do not change this line.

# Give a name to each of the explanatory variables in IDV. They can
# have up to ten characters including spaces. Put the names in single quotes and separate
# the quotes with semicolons. If IDV=[], then set NAMES=[];
names_idv = ["price"; "opcost"; "range"; "ev"; "hybrid"];

# Starting values
# Specify the starting values for b and w for each random coeffient.
# B contains the first parameter, b, for each random coefficient.
# It is a column vector with the same length as IDV. For distribution 5 (normal with zero mean),
# put 0 for the starting value for the mean. The code will keep it at 0.
# W contains the second parameter, w, for each random coefficient.
# It is a column vector with the same length as IDV.
# Put semicolons between the elements of B and W (so they will be column vectors).

B = [.2; .1; .1; -1.44; .412];
W = [0.01; 0.01; 0.01; 0.01; 0.01];

# FIXED COEFFICIENTS
# List the variables in XMAT that enter with fixed coefficients.
# Put semicolons between the numbers.
# If no fixed coefficients, put IDF=[];

idf = [10; 11];

nf = size(idf, 1); #Number of fixed coefficients. Do not change this line.

# Give a name to each of the variables in IDF.
names_idf = ["hiperf"; "medhiperf"];

# Starting values.
# Specify the starting values for the fixed coefficients F.
# F must have the same length as IDF and have one column.
# Put semicolons between the elements (so F will be a column vector.)
F = [0; 0];

# Type of draws to use in simulation
# 1 = pseudo-random draws
# 2 = standard Halton draws
# 3 = shifted and shuffled Halton draws
# 4 = modified Latin hypercube sampling, shifted and shuffled
drawtype = 4;

# Number of draws from to use per person in simulation.
ndraws = 100;

# Set seed for the random number generator.
seed1 = 14239;

# WEIGHTS.
# Do you want to apply weights to the people?
# Set WANTWGT=1 if you want to apply weights; otherwise set WANTWGT=0.
wantwgt = 0;

# If WANTWGT=1, identify the variable in XMAT that contains the weights.
# This variable can vary over people but must be the same for all rows of
# data for each person. Weights cannot vary over choice situations for
# each person or over alternatives for each choice situation -- only over people.
# The code normalizes the weights such that the sum
# of weights over people is to equal NP (to assure that standard errors
# are correctly calculated.) If WANTWGT=0, set IDWGT=0.
idwgt = 0;

# OPTIMIZATION
# Maximum number of iterations for the optimization routine.
# The code will abort after ITERMAX iterations, even if convergence has
# not been achieved. The default is 1000.
maxiters = 10000;

# Convergence criterion based on the maximum change in parameters that is considered
# to represent convergence. If all the parameters change by less than PARAMTOL
# from one iteration to the next, then the code considers convergence to have been
# achieved. The default is 0.
paramtol = 0.0;

# Convergence criterion based on change in the log-likelihood that is
# considered to represent convergence. If the log-likelihood value changes
# less than LLTOL from one iteration to the next, then the optimization routine
# considers convergence to have been achieved. The default is 0.
lltol = 0.0;

# Can also try minimising the likelihood without providing gradient information. 
# Not sure why one would want to do so, however. I SUGGEST NOT ALTERING THIS.
# If a gradient is not provided, a Hessian is not returned, so no standard errors etc.
usegradient = 0; 

# Do not change the next line. It runs the model.
doit()

