using MCMC_LR_Tests
using Base.Test
using Distributions
using HypothesisTests

# consistent testing
srand(UInt32[0x7f712df6, 0xc9c492e6, 0x4935b672, 0x946d94f9])

@testset "IID normal mean" begin
    N = 10000
    dist = MvNormal(Float64[1, -1, 2], [1,2,3])
    ps = [mean_LR_pvalue(mean(dist), rand(dist, N)') for _ in 1:100]
    t = ExactOneSampleKSTest(ps, Uniform(0,1))
    @test pvalue(t) > 0.05
end

@testset "IID normal cov" begin
    N = 10000
    A = [0.1 0.1 0.0;
         0.0 0.5 0.1;
         0.0 0.0 0.7]
    Σ = A'*A
    dist = MvNormal(Float64[1, -1, 2], Σ)
    @test quantile([cov_LR_pvalue(Σ, rand(dist, N)') for _ in 1:100], 0.1) ≥ 0.1
end

@testset "IID normal cov 2" begin
    Σ = [6.366588528827408 3.25842202827805 0.9729856844746686 3.858746597584623;
         3.25842202827805 3.938997699450515 -2.809133292920283 1.5566006096285196;
         0.9729856844746686 -2.809133292920283 7.344095307429978 -1.2388040483400584;
         3.858746597584623 1.5566006096285196 -1.2388040483400584 9.227552651298026]
    dist = MvNormal(zeros(4), Σ)
    @test quantile([cov_LR_pvalue(Σ, rand(dist, N)') for _ in 1:100], 0.1) ≥ 0.1
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
    p10 = quantile([mean_LR_pvalue(μ, simulate_VAR1(A, Σ, N)) for _ in 1:100], 0.1)
    @test p10 ≥ 0.10
end
