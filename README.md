
Sparse Subspace Constrained Partial Least Squares
=================================================

`sscpls` is an R package that provides an implementation of an ADMM algorithm for sparse subspace constrained PLS.

### Install instructions

Use devtools to install:

``` r
#library(devtools)
#install_github("matt-sutton/sgspls")
```

#### Example 1 (Orthogonal components)

The PLS direction vectors are sparse using ℓ<sub>1</sub> penalty term and must lie in a subspace which ensurse orthogonality of components. See below for an example.

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

The objective and orthogonality constraints:

``` r
M <- crossprod(X, y)
# objective maximised
diag(crossprod(t(M)%*%fit$x.weights)) 
```

    ## [1] 1741.6414  825.1240  164.8108  151.3136

``` r
# constraint satisfied
round(crossprod(fit$x.scores),digits = 4)  
```

    ##      [,1] [,2] [,3] [,4]
    ## [1,]    1    0    0    0
    ## [2,]    0    1    0    0
    ## [3,]    0    0    1    0
    ## [4,]    0    0    0    1

### Replication of simulation study

Replication of the simulation study of the subspace constrained PLS paper. The simulation uses the package simrel for generating data from a sparse subspace.

See the subfolder replication.

<img src="replication/results/figure.png" width="3529" />
