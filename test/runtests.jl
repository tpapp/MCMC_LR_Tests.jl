using MCMC_LR_Tests
using Base.Test
using Distributions
using HypothesisTests

# consistent testing
srand(UInt32[0x7f712df6, 0xc9c492e6, 0x4935b672, 0x946d94f9])

@testset "IID normal mean" begin
    N = 10000
    dist = MvNormal(Float64[1, -1, 2], [1,2,3])
    LR_dist = mean_LR_distribution(N, 3)
    ss = [mean_LR_statistic(mean(dist), rand(dist, N)') for _ in 1:100]
    t = ExactOneSampleKSTest(ss, mean_LR_distribution(N, 3))
    @test pvalue(t) > 0.05
end

@testset "IID normal cov" begin
    N = 10000
    A = [0.1 0.1 0.0;
         0.0 0.5 0.1;
         0.0 0.0 0.7]
    Σ = A'*A
    dist = MvNormal(Float64[1, -1, 2], Σ)
    LR_dist = cov_LR_distribution(N, 3)
    ss_max = quantile([cov_LR_statistic(Σ, rand(dist, N)') for _ in 1:100], 0.9)
    @test ss_max ≤ quantile(LR_dist, 0.90)
end

"Simulate `N` draws of a vector AR(1) process `x = A x + ϵ`, where `ϵ ∼ MvNormal(0, Σ)`."
function simulate_VAR1(A, Σ, N)
    K = Base.LinAlg.checksquare(Σ)
    @assert K == Base.LinAlg.checksquare(A)
    @assert maximum(abs.(eigvals(A))) < 1
    ϵ = MvNormal(zeros(K), Σ)
    x = (I - A) \ rand(ϵ)
    X = similar(x, (N, K))
    for i in 1:N
        x .= A*x .+ rand(ϵ)
        X[i, :] .= x
    end
    X
end
        
@testset "VAR(1)" begin
    A = [0.1 0.1 0.0;
         0.0 0.5 0.1;
         0.0 0.0 0.7]
    Σ = A'*A
    μ = zeros(3)
    N = 10000
    ss_max = quantile([mean_LR_statistic(μ, simulate_VAR1(A, Σ, N)) for _ in 1:100], 0.9)
    @test ss_max ≤ quantile(mean_LR_distribution(N, 3), 0.90)
end
