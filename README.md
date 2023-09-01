# RCPsurv

RCPsurv (which stands for <ins>**Random**</ins> <ins>**Change-Point**</ins> for <ins>**surv**</ins>ival outcome) is a package that performs semiparametric estimation and inference for right-censored data with a random change-point using the method proposed by Lee and Wong (2023+).

# How to import the Functions #
> install.packages("devtools")<br />
> library(devtools) <br /> 
> source_url("https://github.com/lcyjames/RCPsurv/blob/73b2b72df2abc89c62ed1618bfe11186524c8749/RCPsurv.R?raw=TRUE")



# Usage #
The package contains 2 functions:
|Functions  | Description|
|------------- | -------------|
RCPsurvSIM  | Generate a data set according to the simulation study in Lee and Wong. (2023+)
RCPsurvEST  | Perform the semiparametric estimation methods of Lee and Wong. (2023+)

<ins>**RCPsurvSIM**</ins>

```
RCPsurvSIM(seed=NA, n, gamma, beta, alpha1, alpha2, mu, sigma)
```
This function generates a data set according to the model of the simulation study in Lee and Wong (2023+) that takes the arguments:
>- `n` is the sample size
>- `gamma` is the coefficient of X
>- `beta` is the 'baseline' coefficient of Z
>- `alpha1` is an intercept added to the regression when Z is greater than the random change-point
>- `alpha2` is the coefficient of Z added to the regression when Z is greater than the random change-point
>- `mu` is the mean of the change-point distribution
>- `sigma` is the standard deviation of the change-point distribution

Example:
```
#This is the setting in Scenario I
Data <- RCPsurvSIM(seed = 1234, n = 500, gamma = 0.5, beta = -1, alpha1 = 2, alpha2 = 1.5, mu = 1.5, sigma = 0.5)
head(Data)

#   id          Yi cen            X        Z
# 1  1 0.004496932   0  1.334176034 1.435799
# 2  2 0.015786129   1 -1.121502497 3.819158
# 3  3 0.026820595   1  1.674326510 3.091647
# 4  4 0.026929615   0 -0.005234058 0.540405
# 5  5 0.030748594   1  0.011925135 3.015005
# 6  6 0.039647808   1  0.129085401 3.971592
```

This data structure is as follows:
>- `id` is the sample identifier
>- `Yi` is the exact failure time or censoring time
>- `cen` is the right-censoring indicator
>- `X` is the non-change-point covariate, which can have multiple columns
>- `Z` is the univariate change-point variable

<ins>**RCPsurvEST**</ins>

```
RCPsurvEST(data, P, m=10, tolerance=10^{-3}, gamma0=NA, beta0=NA, alpha10=NA, alpha20=NA, mu0=NA, sigma0=NA, TRACE=FALSE)
```
This function performs the semiparametric estimation methods of Lee and Wong (2023+). The details of the arguments are as follows:
>- `data` is a data.frame object shown in the above, with columns `id`, `Yi`, `cen`, `X[1]`,...,`X[P]`, `Z`
>- `P` is the dimension of covariate X, which is also equal to the dimension of gamma0
>- `m` is the number of nodes used in the Gaussian quadrature rule for truncated normal distributions
>- `tolerance` is the stopping criterion for the EM algorithm, set to 10^{-3} by default
>- `gamma0` is a vector of constants of size `P` for the initial values of parameter gamma, set to be rep(0,P) by default (gamma0=NA)
>- `beta0` is a constant for the initial value of parameter beta, set to be 0 by default (beta0=NA)
>- `alpha10` is a constant for the initial value of parameter alpha1, set to be 0 by default (alpha10=NA)
>- `alpha20` is a constant for the initial value of parameter alpha2, set to be 0 by default (alpha20=NA)
>- `mu0` is a constant for the initial value of parameter mu, set to be the median of `Z` in `data` by default (mu0=NA)
>- `sigma0` is a constant for the initial value of parameter sigma, set to be 2 by default (sigma0=NA)
>- `TRACE` is an option for tracking the converging path of the parameter estimation, set to be FALSE by default

Example:
```
Data<-RCPsurvSIM(seed = 1234, n = 500, gamma = 0.5, beta = -1, alpha1 = 2, alpha2 = 1.5, mu = 1.5, sigma = 0.5)
Result <-RCPsurvEST(data = Data, P = 1, gamma0 = 0.5, beta0 = -1, alpha10 = 2, alpha20 = 1.5, mu0 = 1.5, sigma0 = 0.5,TRACE = F)
Result

# $loglik
# [1] -2456.817
# 
# $gamma.hat
# [1] 0.4504628
# 
# $beta.hat
# [1] -0.788148
# 
# $alpha1.hat
# [1] 1.751743
# 
# $alpha2.hat
# [1] 1.418709
# 
# $mu.hat
# [1] 1.603879
# 
# $sigma.hat
# [1] 0.3963745
# 
# $gamma.hat.se
# [1] 0.06527598
# 
# $beta.hat.se
# [1] 0.205753
# 
# $alpha1.hat.se
# [1] 0.2864374
# 
# $alpha2.hat.se
# [1] 0.2920645
# 
# $mu.hat.se
# [1] 0.1028265
# 
# $sigma.hat.se
# [1] 0.08603625
# 
# $conv
# [1] 0
```

# Contact #
Lee Chun Yin, James <<james-chun-yin.lee@polyu.edu.hk>>

# Reference #
Lee, C. Y.,& Wong, K. Y. (2023+). Survival analysis with a random change-point. Statistical Methods in Medical Research [online], DOI: [10.1177/09622802231192946](https://doi.org/10.1177/09622802231192946).
