# RCPsurv

RCPsurv (which stands for <ins>**Random**</ins> <ins>**Change-Point**</ins> for <ins>**surv**</ins>ival outcome) is a package that performs semiparametric estimation and inference for right-censored data with a random change-point using the method proposed by Lee and Wong (202X) <DOI: [xx.xxxx/xxxx](https://doi.org/xxxx/xxxx)>.

# How to import the Functions #
> install.packages("devtools")<br />
> library(devtools) <br /> 
> source_url("https://github.com/lcyjames/JointCSsurv/blob/50057f050334e59f038596c898d492be1dbb708a/JointCSsurv.R?raw=TRUE")



# Usage #
The package contains 2 functions:
|Functions  | Description|
|------------- | -------------|
RCPsurvSIM  | Generate a data set according to the simulation study in Lee and Wong. (202X)
RCPsurvEST  | Perform the semiparametric estimation methods of Lee and Wong. (202X)

<ins>**RCPsurvSIM**</ins>

```
RCPsurvSIM(seed=NA, n, gamma, beta, alpha1, alpha2, mu, sigma)
```
This function generates a data set according to the model under scenario I of the simulation study in Lee and Wong (202X) that takes the arguments:
>- `n` is the sample size
>- `gamma` is the coefficient of X
>- `beta` is the 'baseline' coefficient of Z
>- `alpha1` is an intercept added to the regression when Z is greater than the random change-point
>- `alpha2` is the coefficient of Z added to the regression when Z is greater than the random change-point
>- `mu` is the mean of the change-point distribution
>- `sigma` is the standard deviation of the change-point distribution

Example:
```
data <- RCPsurvSIM(seed = 1234, n = 500, gamma = 0.5, beta = -1, alpha1 = 2, alpha2 = 1.5, mu = 1.5, sigma = 0.5)
head(data)

#   id Ti           cen X           Z
# 1  1 0.004496932   0  1.334176034 1.435799
# 2  2 0.015786129   1 -1.121502497 3.819158
# 3  3 0.026820595   1  1.674326510 3.091647
# 4  4 0.026929615   0 -0.005234058 0.540405
# 5  5 0.030748594   1  0.011925135 3.015005
# 6  6 0.039647808   1  0.129085401 3.971592
```

This data structure is as follows:
>- `id` is the sample identifier
>- `cs` is the size within a specific cluster
>- `Lij` is the left endpoint of an observed interval, which takes the value NA for right-censored observations
>- `Rij` is the right endpoint of an observed interval, which takes the value NA for left-censored observations
>- `DL` is the left censoring indicator
>- `DI` is the interval censoring indicator
>- `X` is a covariate in the proportional hazards model, which can have multiple columns
>- `Z` is a covariate in the binomial model without an intercept, which can have multiple columns


<ins>**JointCSsurvEST**</ins>

```
JointCSsurvEST(data, K = 7, P, Q, deg = 3, max.m, M = 20, tolerance = 10^{-3}, 
               gam_0 = NA, beta_0 = NA, alpha_0 = NA, kappa_0 = NA, sigma_0 = NA, TRACE = FALSE)
```
This function performs the semiparametric estimation methods of Lee et al. (2022). The details of the arguments are as follows:
>- `data` is a data.frame object shown in the above, with columns `id`, `cs`, `Lij`, `Rij`, `DL`,`DI`,`X[1]`,...,`X[P]`,`Z[1]`,...,`Z[Q-1]`
>- `K` is the dimension of parameter gamma; the order of the I-splines equal to (`K`-`deg`+1); `K` must be greater than or equal to `deg`
>- `P` is the dimension of covariate X in the proportional hazards model
>- `Q` is the dimension of covariate Z (without an intercept) plus 1 in the binomial model 
>- `deg` is the degree of polynomial used in the I-splines, set to 3 by default
>- `max.m` is the maximum cluster size in the binomial distribution
>- `M` is the number of nodes used in adaptive Gauss-Hermite quadrature, set to 20 by default
>- `tolerance` is the stopping criterion for the EM algorithm, set to 10^{-3} by default
>- `gam_0` is a vector of positive constants of size `K` for the initial values of gamma, set to be rep(2,K) by default (gam_0=NA)
>- `beta_0` is a vector of constants of size `P` for the initial values of parameter beta, set to be rep(0,P) by default (beta_0=NA)
>- `alpha_0` is a vector of constants of size `Q` for the initial values of parameter alpha, set to be rep(0,Q) by default (alpha_0=NA)
>- `kappa_0` is a constant for the initial values of parameter kappa, set to be 0 by default (kappa_0=NA)
>- `sigma_0` is a constant for the initial values of parameter sigma, set to be 2 by default (sigma_0=NA)
>- `TRACE` is an option for tracking the converging path of the parameter estimation, set to be FALSE by default

Example:
```
Dataset<-JointCSsurvSIM(seed = 1234, n = 50, m = 10, beta = 1, alpha = c(1,log(2)), kappa = -0.5, sigma= 1)
Result <-JointCSsurvEST(data = Dataset, K = 7, P = 1, Q = 2, deg = 3, max.m = 10, tolerance = 10^{-3}, M = 20, TRACE = FALSE)
Result

# $loglik
# [1] -943.2061
# 
# $gam.hat
# [1] 0.1006968 0.0913659 0.4115768 0.5069621 0.8006600 0.3812705 0.5546317
# 
# $alpha.hat
# [1] 1.09652 0.53590
# 
# $kappa.hat
# [1] -0.5807313
# 
# $beta.hat
# [1] 0.901442
# 
# $sigma.hat
# [1] 0.8442931
# 
# $alpha.hat.se
# [1] 0.1303666 0.1338147
# 
# $kappa.hat.se
# [1] 0.173564
# 
# $beta.hat.se
# [1] 0.07883315
# 
# $sigma.hat.se
# [1] 0.1251456
```

# Contact #
Lee Chun Yin, James <<james-chun-yin.lee@polyu.edu.hk>>

# Reference #
Lee, C. Y., Wong, K. Y., Lam, K. F.,& Bandyopadhyay, D. (2022). A semiparametric joint model for cluster size and subunit-specific interval-censored
outcomes. Biometrics [online], DOI: [10.1111/biom.13795](https://doi.org/10.1111/biom.13795).
