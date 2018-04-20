
Sparse Subspace Constrained Partial Least Squares
=================================================

`sscpls` is an R package that provides an implementation of an ADMM algorithm for sparse subspace constrained PLS.

### Install instructions

Use devtools to install:

``` r
#library(devtools)
#install_github("sscpls", "matt-sutton")
```

#### Example 1 (Orthogonal components)

The PLS direction vectors are sparse using ℓ<sub>1</sub> penalty term and must lie in a subspace which ensurse orthogonality of components. See below for example and comparison to mixOmics (uses naive SPLS-NIPALS deflation).

``` r
set.seed(1)
library(sscpls)
n<-50 # number of observations
p<-100 # number of variables
X<-matrix(rnorm(n*p),ncol=p)
X <- scale(X)
beta <- c(1:5, rep(0, p-5))
y<-X%*%beta + matrix(rnorm(n))

fitcv <- cv_sscpls(X, y, K = 1:4, lambda = 1:9/10, fold = 10)
fit <- sscpls(X, y, ncomp = 4, lambda = 0.9, scale = F)

Beta_comps <- matrix(unlist(fit$Beta), ncol = 4) # Get Matrix of beta estimates at each component
show_nonzero(cbind(Beta_comps, beta))
```

    ##                                          beta
    ## [1,] 0.000000 0.000000 0.000000 1.229254    1
    ## [2,] 0.000000 0.000000 1.835830 1.887444    2
    ## [3,] 1.747846 1.872779 1.872779 3.055013    3
    ## [4,] 0.000000 4.074708 4.112529 4.112529    4
    ## [5,] 5.480903 5.480903 5.409163 4.873663    5

``` r
# sum(cov(X%*%fit$Uorg, y)) # objective maximised
# round(crossprod(fit$x.scores),digits = 4)  # constraint satisfied
```

#### Example 2 (Compositional data)

PLS direction vectors for compositional data. Following the process for non-sparse PLS on compositional data we use the centered log transform of the data.

``` r
##-- centered log transform --#
clrt <- function(x){
  X <- t(log(x))
  X <- t(scale(X, center = T, scale = F))
  return(X)
}
```

Simulate sparse data as discribed in Li 2014.

``` r
##-- Short simulation as in Li --#
set.seed(1)
p <- 60; n<- 50

Sigma <- matrix(0.5, nrow = p, ncol = p)
for(i in 1:p) { for(j in 1:p){Sigma[i,j] <- 0.5^abs(i-j)} }

theta <- matrix(0, nrow = 1, ncol = p)
theta[1:5] <- log(0.5*p)

W <- MASS::mvrnorm(n = n, mu = theta, Sigma = Sigma)
X <- exp(W)
X <- t(apply(X, 1, function(x) x/sum(x)))
# rowSums(X) # All lie in the simplex >0 and sum to 1

Z <- t(apply(X, 1, function(z) log(z/z[p])))

Betastar <- matrix(c(1:5, rep(0, p - 5)))
Betastar[p] <- -sum(Betastar)
error <- matrix(rnorm(n, sd = 0.5))

y <- Z%*%Betastar + error
```

Standard compositional PLS (see Hinkle 1995). Implemented without sparsity.

``` r
Xtilde <- clrt(X)
library(mixOmics)
fit_mix <- pls(Xtilde, y, ncomp = 5, scale = F, mode = "regression")
fit_B <- predict(fit_mix, newdata = fit_mix$X)$B.hat

cbind(fit_B[,,4], Betastar)[1:15,] # only show some all are nonzero
```

    ##            [,1] [,2]
    ## X1   1.10265919    1
    ## X2   1.00540521    2
    ## X3   2.95145633    3
    ## X4   2.90012744    4
    ## X5   3.80178650    5
    ## X6   1.79847952    0
    ## X7   0.24279385    0
    ## X8  -0.62521063    0
    ## X9  -1.45060643    0
    ## X10 -0.25444222    0
    ## X11  0.68327214    0
    ## X12  0.06776761    0
    ## X13 -0.23007247    0
    ## X14 -0.59405261    0
    ## X15  0.22940184    0

``` r
sum(fit_B[,,4]) # Compositional PLS satisfies sum(B) = 0 
```

    ## [1] 3.849698e-14

Now using sparsity

``` r
fit <- sscpls(Xtilde, y, ncomp = 2, lambda = 0.45, scale = F, compositional = T, max_itter = 10^4, abstol = 1e-7, reltol = 1e-7)

Beta_comps <- matrix(unlist(fit$Beta), ncol = 2) # Get Matrix of beta estimates at each component
show_nonzero(cbind(Beta_comps,Betastar))
```

    ##             [,1]        [,2] [,3]
    ## [1,]   0.0000000   0.0000000    1
    ## [2,]   0.0000000   3.0390226    2
    ## [3,]   0.3300458   2.2246210    3
    ## [4,]   3.1765315   3.1765315    4
    ## [5,]   6.8961250   6.8961250    5
    ## [6,]   2.9082911  -0.5731075    0
    ## [7,] -13.3111062 -14.7632668  -15

``` r
sum(fit$Beta[[2]]) # Compositional PLS satisfies sum(B) = 0  
```

    ## [1] -7.429346e-05

Alternative Sparse PLS methods are not able to retain the sum of beta is zero constraint.

Bibliography
------------

Lin, Wei, Pixu Shi, Rui Feng, and Hongzhe Li. 2014. “Variable Selection in Regression with Compositional Covariates.” Biometrika 101 (4). Oxford University Press: 785–97.

Hinkle, John, and William Rayens. 1995. “Partial Least Squares and Compositional Data: Problems and Alternatives.” Chemometrics and Intelligent Laboratory Systems 30 (1): 159–72.
