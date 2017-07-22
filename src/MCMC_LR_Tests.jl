module MCMC_LR_Tests

using ArgCheck
using Distributions
using MCMCDiagnostics

export mean_LR_statistic, mean_LR_distribution, cov_LR_statistic, cov_LR_distribution

"""
    column_ρs(X)

For observations in the rows of `X`, return the estimated `ρ=ESS/n` for each column (variable).

Provides the default value to [`mean_LR_statistic`](@ref), otherwise only useful for debugging.
"""
column_ρs(X::AbstractMatrix) = vec(mapslices(first ∘ ess_factor_estimate, X, 1))

"""
    mean_LR_statistic(μ, X, [ρ])

Likelihood ratio statistic for testing that the columnwise mean of `X` are equal to the vector `μ`, if the rows of `X` are normally distributed (with unknown variance, estimated from the sample).

Assume that some form of a CLT holds, so

``n(x̄(n)-x̄) → N(0, Σ)`` in distribution, where `x̄(n)` is the running mean.

The test statistic is

``(x̄(n)-μ)' Σ⁻¹ (x̄(n)-μ) ∼ T²(k, n-1)``

where `T²` is [Hotelling's T²-distribution](https://en.wikipedia.org/wiki/Hotelling%27s_T-squared_distribution), if `x̄ = μ` (Dootika et al 2017).

For IID samples, `Σ` can be estimated from the sample. However, MCMC draws are autocorrelated, so this is more difficult. Let ρᵢ = ESSᵢ/n denote the scaling factor to obtain the effective sample size. For a crude approximation, we rescale the variance by ρᵢ, and covariances by √(ρᵢρⱼ).

A rescaled version of the test statistic is returned. Compare it to the result of [`mean_LR_distribution`](@ref).
"""
function mean_LR_statistic(μ::AbstractVector, X::AbstractMatrix, ρ = column_ρs(X))
    N, K = size(X)
    @argcheck K == length(μ) == length(ρ)
    R = Diagonal(.√ρ)           # correction by ρ = ESS/n
    S = cov(X, 1, false)        # sample covariance
    δ = vec(mean(X, 1)) - μ
    h = (δ' * (R * inv(S) * R) * δ) * N # has Hotelling T²(K, N-1) distribution
    h*((N-K)/((N-1)*K))                 # convert to F(K, N-K) distribution
end

"""
    mean_LR_distribution(N, K)

Return a distribution for [`mean_LR_statistic`](@ref) under the null hypothesis `x̄=μ`, where `N` is the sample size and `K` is the dimension.

Use of a Kolmogorov-Smirnov test is recommended, note however that the distribution is slightly different because of the ESS rescaling, so be conservative with the p-values in unit tests.
"""
mean_LR_distribution(N::Int, K::Int) = FDist(K, N-K)

"""
    cov_LR_statistic(Σ, X)

Likelihood ratio statistic for comparing the covariance of `X` (with observations in rows) to `Σ`. See [`cov_LR_distribution`](@ref).
"""
function cov_LR_statistic(Σ::AbstractMatrix, X::AbstractMatrix)
    N, K = size(X)
    @argcheck size(Σ) == (K, K)
    U = chol(Symmetric(Σ))
    S = cov(X / U, 1)
    L = trace(S) - logdet(S) - K
    N*L
end

"""
   cov_LR_distribution(N, K)

Distribution for the test statistic of [`cov_LR_statistic`](@ref), where the data has dimension `K` and `N` observations.

Note that this is an asymptotic statistic, so it has better coverage when `N ≫ K`.
"""
cov_LR_distribution(N::Int, K::Int) = Chisq(K*(K+1)/2)

end # module
