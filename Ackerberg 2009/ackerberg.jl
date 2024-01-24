using BenchmarkTools, Distributions, LinearAlgebra, StatsBase, Optim

### Data setup 

N = 5000 # number of simulated individuals
J = 3    # number of inside choices
R = 100  # simulated shocks per individual


μ = [-1, .5]
Σ = [.5 0; 0 .9]

βᵢ = rand(MvNormal(μ, Σ), N)

pᵢⱼ = rand(Exponential(1), N, J) # characteristics vary at the individual level.
xᵢⱼ = [2.2 3.8 1.4] .+ rand(Exponential(2.5), N, J)

sᵢⱼ = zeros(Float64, N, J + 1)

# choice probabilities under the model 
vᵢⱼ = βᵢ[1, :] .* pᵢⱼ + βᵢ[2, :] .* xᵢⱼ
sᵢⱼ = [ exp.(vᵢⱼ) ones(Float64, N) ]
sᵢⱼ = sᵢⱼ ./ sum(sᵢⱼ, dims = 2)

σⱼ = mean(sᵢⱼ, dims = 1)

# observed choices
yᵢⱼ = zeros(Float64, N, J + 1)
for i ∈ 1:N
    yᵢⱼ[i, sample(collect(1:(J + 1)), Weights(sᵢⱼ[i, :]), 1)] .= 1.0 
end

# objects useful below 
pᵢᵣⱼ = repeat(pᵢⱼ, inner = (R, 1))
xᵢᵣⱼ = repeat(xᵢⱼ, inner = (R, 1))
βᵢᵣ = zeros(Float64, 2, N * R)

####################################
### maximum simulated likelihood ###
####################################

# T(a, b) = [exp(max(a, -4.0)) 0.0; 0.0 exp(max(b, -4.0))] * [exp(max(a, -4.0)) 0.0; 0.0 exp(max(b, -4.0))]'
T(a, b) = [exp(a) 0.0; 0.0 exp(b)] * [exp(a) 0.0; 0.0 exp(b)]'
T(θ) = T(θ[3], θ[4])

oᵢᵣ = rand(MvNormal([0.0, 0.0], [1.0 0.0; 0.0 1.0]), N * R)

function sᵢⱼ_sim(ϑ, βᵢᵣ) 
    βᵢᵣ  .= ϑ[1:2] .+ T(ϑ) * oᵢᵣ 
    @views sᵢᵣⱼ  = [exp.(βᵢᵣ[1, :] .* pᵢᵣⱼ + βᵢᵣ[2, :] .* xᵢᵣⱼ) ones(Float64, N * R)]
    sᵢᵣⱼ .= sᵢᵣⱼ ./ sum(sᵢᵣⱼ, dims = 2) 
    return dropdims(mean(reshape(sᵢᵣⱼ, R, N, J + 1), dims = 1), dims = 1)
end 

lnℒ_msl(ϑ) = - sum( log.( sum( yᵢⱼ .* sᵢⱼ_sim(ϑ, βᵢᵣ), dims = 2 ) ) )

MSL = optimize(lnℒ_msl, [ -0.25 0.25 -.7 -.7 ], NelderMead(), 
               Optim.Options(g_tol = 1e-12, time_limit = 100, show_trace = true, show_every = 10))
θMSL = Optim.minimizer( MSL )
[θMSL[1:2]... diag(T(θMSL))...]
[μ... diag(Σ)...] # diag(Σ) terms hard to get right  

########################################################
### importance sampling with normal proposal density ###
########################################################

function sᵢᵣⱼ_is(ϑ, oᵢᵣ)
    βᵢᵣ  .= oᵢᵣ
    @views sᵢᵣⱼ = [exp.(βᵢᵣ[1, :] .* pᵢᵣⱼ + βᵢᵣ[2, :] .* xᵢᵣⱼ) ones(Float64, N * R)]
    sᵢᵣⱼ .= sᵢᵣⱼ ./ sum(sᵢᵣⱼ, dims = 2) 
    𝒻ᵢᵣ   = pdf(MvNormal([ϑ[1], ϑ[2]], T(ϑ)), oᵢᵣ) # density under proposal density
    return sᵢᵣⱼ, 𝒻ᵢᵣ
end 
 
function sᵢⱼ_is(ϑ, sᵢᵣⱼ, 𝒻ᵢᵣ, oᵢᵣ)
    ℊᵢᵣ = pdf(MvNormal([ϑ[1], ϑ[2]], T(ϑ)), oᵢᵣ)    
    𝓈ᵢᵣⱼ = sᵢᵣⱼ .* ℊᵢᵣ ./ 𝒻ᵢᵣ
    return dropdims(mean(reshape(𝓈ᵢᵣⱼ, R, N, J + 1), dims = 1), dims = 1)
end 

lnℒ_is(ϑ, yᵢⱼ, sᵢᵣⱼ, 𝒻ᵢᵣ, oᵢᵣ) = - sum( log.( sum( yᵢⱼ .* sᵢⱼ_is(ϑ, sᵢᵣⱼ, 𝒻ᵢᵣ, oᵢᵣ), dims = 2 ) ) )

function importancesampling(proposal_theta)
  
    rθ = round.([proposal_theta[1:2]... diag(T(proposal_theta))...], digits = 3)
    print("Proposal distribution parameters: $rθ.\n")

    oᵢᵣ = rand(MvNormal([proposal_theta[1], proposal_theta[2]], T(proposal_theta)), N * R) # draws from proposal density
  
    sᵢᵣⱼ, 𝒻ᵢᵣ = sᵢᵣⱼ_is(proposal_theta, oᵢᵣ)
      
    IS = optimize(ϑ -> lnℒ_is(ϑ, yᵢⱼ, sᵢᵣⱼ, 𝒻ᵢᵣ, oᵢᵣ), [ -0.25 0.25 -.7 -.7 ], NelderMead(), 
                  Optim.Options(g_tol = 1e-12, time_limit = 100, show_trace = true, show_every = 10))
    ISt = Optim.minimizer(IS)

    return [ISt[1:2]... diag(T(ISt))...]
end 

IS1 = importancesampling([ 0.0 0.0 0.0 0.0 ])
IS2 = importancesampling([ -0.25 0.25 log(0.5) log(0.5) ])
IS3 = importancesampling([μ... log.(diag(cholesky(Σ).factors))...])

[μ... diag(Σ)...]
IS1 
IS2   # a good proposal distribution is extremely important for the performance of the estimator.
IS3 

######################################################
### iteratively updating the proposal distribution ###
######################################################

function importancesampling_iterative(proposal_theta; time_limit = 100000.0, nupdates = 1)
  
    time_run = 0.0 
    g_converged = false
    θ = proposal_theta
    iter = 0

    uᵢᵣ = rand(MvNormal([0.0, 0.0], [1.0 0.0; 0.0 1.0]), N * R) # draws from N(0, I) to be held fixed  
    
    while (time_run < time_limit) & (g_converged == false) & (iter < nupdates)

        rθ = round.([θ[1:2]... diag(T(θ))...], digits = 3)
        print("Proposal distribution updated to $rθ. Runtime: $time_run. \n")

        oᵢᵣ = θ[1:2] .+ T(θ) * uᵢᵣ # updating draws from proposal density using affine transformation
    
        sᵢᵣⱼ, 𝒻ᵢᵣ = sᵢᵣⱼ_is(θ, oᵢᵣ)
            
        IS = optimize(ϑ -> lnℒ_is(ϑ, yᵢⱼ, sᵢᵣⱼ, 𝒻ᵢᵣ, oᵢᵣ), [ -0.25 0.25 -.7 -.7 ], NelderMead(), 
                      Optim.Options(g_tol = 1e-12, time_limit = time_limit, show_trace = true, show_every = 5))
        θ = Optim.minimizer( IS )
        
        g_converged = IS.g_converged
        time_run += Optim.time_run(IS)
        iter += 1
    end 

    # after nupdates resets, let the last optimization run to conclusion 
    rθ = round.([θ[1:2]... diag(T(θ))...], digits = 3)
    print("Proposal distribution finally updated to $rθ. Runtime: $time_run. \n")

    oᵢᵣ = θ[1:2] .+ T(θ) * uᵢᵣ

    sᵢᵣⱼ, 𝒻ᵢᵣ = sᵢᵣⱼ_is(θ, oᵢᵣ)
    
    IS = optimize(ϑ -> lnℒ_is(ϑ, yᵢⱼ, sᵢᵣⱼ, 𝒻ᵢᵣ, oᵢᵣ), θ, NelderMead(), 
                  Optim.Options(g_tol = 1e-12, time_limit = 2 * time_limit, show_trace = true, show_every = 5))
    θ = Optim.minimizer( IS )

    return θ, [θ[1:2]... diag(T(θ))...]
end 

ISI1, ISI1t = importancesampling_iterative([ 0.0 0.0 0.0 0.0 ])
ISI3, ISI3t = importancesampling_iterative([μ... log.(diag(cholesky(Σ).factors))...])

ISI1t
ISI3t