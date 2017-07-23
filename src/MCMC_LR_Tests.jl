module MCMC_LR_Tests

using ArgCheck
using Distributions
using MCMCDiagnostics

export mean_LR_pvalue, cov_LR_pvalue

"""
    column_ρs(X)

For observations in the rows of `X`, return the estimated `ρ=ESS/n` for each column (variable).

Provides the default value to [`mean_LR_statistic`](@ref), otherwise only useful for debugging.
"""
column_ρs(X::AbstractMatrix) = vec(mapslices(first ∘ ess_factor_estimate, X, 1))

"""
    mean_LR_pvalue(μ, X, [ρ])

P-value for a likelihood ratio statistic for testing that the columnwise mean of `X` are equal to the vector `μ`, if the rows of `X` are normally distributed (with unknown variance, estimated from the sample).

Assume that some form of a CLT holds, so

``n(x̄(n)-x̄) → N(0, Σ)`` in distribution, where `x̄(n)` is the running mean.

The test statistic is

``(x̄(n)-μ)' Σ⁻¹ (x̄(n)-μ) ∼ T²(k, n-1)``

where `T²` is [Hotelling's T²-distribution](https://en.wikipedia.org/wiki/Hotelling%27s_T-squared_distribution), if `x̄ = μ` (Dootika et al 2017).

For IID samples, `Σ` can be estimated from the sample. However, MCMC draws are autocorrelated, so this is more difficult. Let ρᵢ = ESSᵢ/n denote the scaling factor to obtain the effective sample size. For a crude approximation, we rescale the variance by ρᵢ, and covariances by √(ρᵢρⱼ).
"""
function mean_LR_pvalue(μ::AbstractVector, X::AbstractMatrix, ρ = column_ρs(X))
    N, K = size(X)
    @argcheck K == length(μ) == length(ρ)
    R = Diagonal(.√ρ)           # correction by ρ = ESS/n
    S = cov(X, 1, false)        # sample covariance
    δ = vec(mean(X, 1)) - μ     # deviations from μ
    h = (δ' * (R*inv(S)*R) * δ) * N    # has Hotelling T²(K, N-1) distribution
    ccdf(FDist(K, N-K), h*((N-K)/((N-1)*K))) # convert to F(K, N-K) distribution
end

"""
    cov_LR_pvalue(Σ, X)

P-value for likelihood ratio statistic for comparing the covariance of `X` (with observations in rows) to `Σ`.

Note that this is an *asymptotic* statistic, so it has better coverage when `N ≫ K`. Coverage may be very bad when `N` is small, see Bai et al (2009).
"""
function cov_LR_pvalue(Σ::AbstractMatrix, X::AbstractMatrix)
    N, K = size(X)
    @argcheck size(Σ) == (K, K)
    U = chol(Symmetric(Σ))
    S = cov(X / U, 1)
    L = trace(S) - logdet(S) - K
    ccdf(Chisq(K*(K+1)/2), N*L)
end

end # module
