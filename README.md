
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
p<-50 # number of variables
X<-matrix(rnorm(n*p),ncol=p)
X <- scale(X)
beta <- c(1:5, rep(0, p-5))
y<-X%*%beta + matrix(rnorm(n))

fit <- sscpls(X, y, ncomp = 4, lambda = 0.7, scale = F)

Beta_comps <- matrix(unlist(fit$Beta), ncol = 4) # Get Matrix of beta estimates at each component
show_nonzero(cbind(Beta_comps, beta))
```

    ##                                          beta
    ## [1,] 0.000000 0.000000 0.957490 0.957490    1
    ## [2,] 0.000000 1.387000 1.996100 1.996100    2
    ## [3,] 3.381773 3.381773 3.381773 3.381773    3
    ## [4,] 1.076479 3.865946 3.865946 3.865946    4
    ## [5,] 5.459896 4.758949 4.758949 4.758949    5

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

Betastar <- matrix(c(1, -0.8, 0.6, 0, 0,-1.5,0.5, 1.2, rep(0, p - 8)))
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
    ## X1   0.54536013  1.0
    ## X2  -0.28020479 -0.8
    ## X3   0.42767293  0.6
    ## X4   0.12172604  0.0
    ## X5  -0.33034597  0.0
    ## X6  -0.96705289 -1.5
    ## X7   0.53214948  0.5
    ## X8   0.87387024  1.2
    ## X9   0.02515786  0.0
    ## X10 -0.07873966  0.0
    ## X11 -0.03711091  0.0
    ## X12  0.21280074  0.0
    ## X13 -0.28873248  0.0
    ## X14 -0.26195036  0.0
    ## X15 -0.20454224  0.0

``` r
sum(fit_B[,,4]) # Compositional PLS satisfies sum(B) = 0 
```

    ## [1] 2.709638e-15

Now using sparsity

``` r
fit <- sscpls(Xtilde, y, ncomp = 3, lambda = 0.85, scale = F, compositional = T)

Beta_comps <- matrix(unlist(fit$Beta), ncol = 3) # Get Matrix of beta estimates at each component
show_nonzero(cbind(Beta_comps,Betastar))
```

    ##            [,1]       [,2]        [,3] [,4]
    ##  [1,]  0.000000  0.5653049  0.56530494  1.0
    ##  [2,]  0.000000  0.0000000 -0.37275514 -0.8
    ##  [3,]  0.000000  0.0000000  0.00000000  0.6
    ##  [4,]  0.000000 -1.1358419 -1.14603966 -1.5
    ##  [5,]  0.000000  0.0000000  0.04524048  0.5
    ##  [6,]  1.034775  1.2002163  1.20021626  1.2
    ##  [7,]  0.000000  0.0000000  0.33773541  0.0
    ##  [8,]  0.000000  0.4051827  0.40518270  0.0
    ##  [9,] -1.034879 -1.0348790 -1.03487895 -1.0

``` r
sum(fit$Beta[[3]]) # Compositional PLS satisfies sum(B) = 0  
```

    ## [1] 6.034395e-06

Alternative Sparse PLS methods are not able to retain the sum of beta is zero constraint.

Bibliography
------------

Lin, Wei, Pixu Shi, Rui Feng, and Hongzhe Li. 2014. “Variable Selection in Regression with Compositional Covariates.” Biometrika 101 (4). Oxford University Press: 785–97.

Hinkle, John, and William Rayens. 1995. “Partial Least Squares and Compositional Data: Problems and Alternatives.” Chemometrics and Intelligent Laboratory Systems 30 (1): 159–72.
