# MCMC\_LR\_Test

[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)
[![Build Status](https://travis-ci.org/tpapp/MCMC_LR_Tests.jl.svg?branch=master)](https://travis-ci.org/tpapp/MCMC_LR_Tests.jl)
[![Coverage Status](https://coveralls.io/repos/tpapp/MCMC_LR_Tests.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/tpapp/MCMC_LR_Tests.jl?branch=master)
[![codecov.io](http://codecov.io/github/tpapp/MCMC_LR_Tests.jl/coverage.svg?branch=master)](http://codecov.io/github/tpapp/MCMC_LR_Tests.jl?branch=master)

Compare the mean and covariance of a MCMC sample to a known value using a likelihood ratio test with a p-value.

For unit testing *MCMC software*. The returned p-values should be treated as a distance measure that has a distribution that is not too far away from the uniform, but the latter is not guaranteed. Simple testing can assert that p-values are bigger than a certain threshold, more complex testing can compare quantiles. See the docstrings of the exported functions `mean_LR_pvalue` and `cov_LR_pvalue`, and the unit tests.

## References

Bai, Z., Jiang, D., Yao, J. F., & Zheng, S. (2009). Corrections to LRT on large-dimensional covariance matrix by RMT. The Annals of Statistics, 3822-3840.

Vats, D., Flegal, J. M., & Jones, G. L. (2015). Multivariate output analysis for Markov chain Monte Carlo. arXiv preprint [arXiv:1512.07713](https://arxiv.org/abs/1512.07713).
