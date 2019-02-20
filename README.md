
<!-- README.md is generated from README.Rmd. Please edit that file -->

# IntegrativeCI

Implements integrative analyses for the average treatment effect
combining main data with unmeasured confounders and validation data with
supplementary information on these
confounders.

### Visit the [package website](https://github.com/shuyang1987/IntegrativeCI)

## Installation with `devtools`:

``` r
devtools::install_github("shuyang1987/IntegrativeCI")
```

## Description

Under the unconfoundedness assumption with completely observed
confounders, the smaller validation data allow for constructing
consistent estimators for causal effects, but the big main data can only
give error-prone estimators in general.

The integrative estimator leverages the information in the big main data
to improve the estimation efficiencies yet preserve the consistencies of
the initial estimators based solely on the validation data.
Specifically, it exploits the difference of the error-prone estimators
applied to the main data and the validation data, which is consistent
for zero assuming that the main data and the validation data are
representative of the same target population.

The framework applies to asymptotically normal estimators, including the
commonly-used regression imputation, (augmented) weighting, and matching
estimators.

### Main Paper: Yang and Ding (2018)

Yang, S. and Ding, P. (2018). Combining multiple observational data
sources to estimate causal effects. (<https://arxiv.org/abs/1801.0080>)

## Use

The function implementing the integrative average treatment effect is
`IntegrativeCI::IntegrativeATE()`.

## Example

``` r

set.seed(1)

n<-1000  # the (combined) main sample size (Samples 1 and 2)
n2<-500  # the validation sample size (Sample 2)

## generate covariates (x,u)

x<-runif(n,0,2)
u<-2*sin(x)+runif(n,0,1)

## generate treatment A

lps<-cbind(1,x,u)%*%c(0,-1,1)
ps<-exp(lps)/(1+exp(lps))
A<-rbinom(n,1,ps)

## generate potential and observed outcomes

loc.a1<-which(A==1)
loc.a0<-which(A==0)
y1<- -x+2*u+rnorm(n,0,1)
y0<- -x+1*u+rnorm(n,0,1)
y<-y1*A+y0*(1-A)
true<-mean(y1-y0)

## generate indicator of membership of the validation sample (Sample 2)

I<-rep(1,n)
I[((1+(n-n2)):n)]<-0

## u is not observed for Sample 1
loc1<-which(I==1)
u[loc1]<-NA

method_val <-c("reg")
method_ep<-c("reg","ipw")

true
#> [1] 1.839584
out<-IntegrativeCI::IntegrativeATE(I,x,u,y,A,method_val,method_ep,nboot=50)
out$est
#>          [,1]
#> [1,] 1.860912
out$ve
#>            [,1]
#> [1,] 0.00753725
out$ve_boot
#>             [,1]
#> [1,] 0.007406224
```
