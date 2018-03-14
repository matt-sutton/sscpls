#' Sparse Subspace Constrained Partial Least Squares
#'
#' Provides an implementation of an ADMM algorithm for computing the sparse
#' subspace constrained PLS estimator. Computes PLS regression where the
#' direction vectors are constrained to be both sparse and lie within a given subspace.
#'
#' @param X A matrix of regressors (n x p). By default the matrix will be
#'   centered to have mean zero.
#' @param Y A matrix of continuous responses (n x q). By default the matrix will
#'   be centered to have mean zero.
#' @param ncomp The number of components to include in the model.
#' @param lambda Thresholding parameter for sparsity
#' @param scale Scale predictors and response by dividing each variable by its sample standard deviation
#' @param center Center predictors and response by subtracting the column means from the variables
#' @param compositional Use either compositional projection or SIMPLS orthogonal projection
#' @param abstol Absolute tolerance for the ADMM algoithm (default is 1e-06)
#' @param reltol Relative tolerance for the ADMM algoithm (default is 1e-06)
#' @param max_itter Maximum number of itterations for the ADMM method
#'
#' @export
#'
#' @return \code{simpls} returns a fit PLS with:
#' \item{weights}{a list containing the X and Y pls weights.}
#' \item{scores}{a list containing the X and Y pls scores.}
#' \item{beta}{a regression coefficient}
#'
#' @examples
#'n<-50 # number of observations
#'p<-50 # number of variables
#'X<-matrix(rnorm(n*p),ncol=p)
#'X <- scale(X)
#'y<-X[,1:5]%*%1:5 + matrix(rnorm(n))
#'
#'res <- sscpls(X, y, ncomp = 5, lambda = 0.6, scale = F, abstol = 1e-07, reltol = 1e-07)
#'
#'# abs and rel tol give a cutoff for the orthogonality constriant, i.e. here precision to 1e-06
#'round(crossprod(res$x.scores), 5)
#'

sscpls <- function(X, Y, ncomp=2, lambda, scale = F, center=T, compositional=F, abstol = 1e-06, reltol= 1e-06, max_itter = 10^5) {

  #-- Internal funciton for vector normalisation --#
  normalise <- function(x, xnorm){
    xnorm <- drop(xnorm + 1*(xnorm == 0))
    return(x/(xnorm))
  }

  #-- Internal funciton to get projection given basis X --#
  get_proj <- function(X){
    X <- as.matrix(X)
    if(all(X == 0)) return(tcrossprod(X))

    L <- t(suppressWarnings(chol(crossprod(X), pivot = TRUE)))
    r <- attr(L, "rank")
    piv <- attr(L, "pivot")
    Qt <- forwardsolve(L, t(X[, piv]), r)

    ## P = QQ'
    H <- crossprod(Qt)
    return(H)
  }

  #-- Data parameters --#
  n <- nrow(X);  p <- ncol(X);   q <- ncol(Y)

  #-- Scaling --#
  X <- scale(X,center = center, scale = scale)
  Y <- scale(Y,center = center, scale = scale)

  #-- Initialisation --#
  Vorg <- V <- matrix(0, p, ncomp);  U <- matrix(0, q, ncomp)
  Z <- matrix(0, n, ncomp);  Z_center <- rep(0, ncomp); Z_scale <- rep(1, ncomp)

  lambdas <- rnorm <- snorm <- NULL

  #-- Compute a projection matrix onto the column space of S --#
  proj_matrices <- array(0, dim = c(p,p,ncomp+1))
  proj_matrices[,,1] <- proj <- diag(1, p,p)

  #-- Compute the svd of M --#
  M <- crossprod(X, Y)
  Mt <- t(M)
  svdm <- svd(Mt, nu = 1, nv = 1)

  #-- initalise c = Mu and lambda --#
  uold <- svdm$u;  d <- svdm$d[1]

  lambda_max <- max(abs(svdm$v))*d
  lambda_h <- lambda_max*lambda
  lambda <- lambda_h/d

  for(h in 1:ncomp){

    #-- find lambda (Witten style Tuning) --#
    d <- svd(proj%*%M)$d[1]
    lambda_h <- d * lambda
    lambdas <- c(lambdas, lambda_h)

    for (i in 1:500){

      Mu <- M%*%uold
      #-- get v direction vector --#
      admm_u <- get_u(Mu, proj, lambda_h, rho=1, abstol = abstol, reltol= reltol, max_itter= max_itter, compositional=compositional)
      rnorm <- c(rnorm, admm_u$rnorm)
      snorm <- c(snorm, admm_u$snorm)
      v <- admm_u$u

      #-- get u direction vector --#
      u <- normalise(Mt%*%v, norm(Mt%*%v,"e"))

      if(norm(u - uold, type = "e") < 1e-10){
        break
      }
      uold <- u
    }
    U[,h] <- u

    #-- compute score --#
    z <- X%*%v

    #-- save away scaled loadings --#
    z_c <- mean(z)
    z <- z - z_c;
    z_s <- drop(sqrt(crossprod(z)))
    Z[,h] <- z <- normalise(z, z_s)
    V[,h] <- v <- normalise(v, z_s)

    #-- construct the unwanted space --#
    cx <- t(X)%*%Z[,1:h, drop = F]

    #-- construct wanted space --#
    proj <- diag(1,p,p) - get_proj(cx)
    proj_matrices[,,h+1] <- proj

  }

  Beta <- sapply(1:ncomp, function(h) tcrossprod(V[,1:h, drop = F])%*%M, simplify = F)

  res = list(x.scores = Z,
             x.weights = round(V,9),
             y.weights = round(U,9),
             Beta = Beta,
             Px = proj_matrices,
             lambdas = lambdas,
             rnorm = rnorm,
             snorm = snorm,
             Vorg= Vorg)
  return(res)
}


