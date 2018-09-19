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

sscpls <- function(X, Y, ncomp=2, lambda, scale = F, center=T, compositional=F, abstol = 1e-06, reltol= 1e-06, max_itter = 500) {

  #-- Data parameters --#
  n <- nrow(X);  p <- ncol(X);   q <- ncol(Y)

  #-- Scaling --#
  X <- scale(X,center = center, scale = scale)
  Y <- scale(Y,center = center, scale = scale)

  #-- Initialisation --#
  V <- Vscaled <- matrix(0, p, ncomp);  U <- matrix(0, q, ncomp)
  Z <- matrix(0, n, ncomp);  Z_center <- rep(0, ncomp); Z_scale <- rep(1, ncomp)

  #-- Compute a projection matrix onto the column space of S --#
  proj_matrices <- array(0, dim = c(p,p,ncomp+1))
  proj_matrices[,,1] <- proj <- diag(1, p,p)

  #-- Compute the svd of M --#
  M <- crossprod(X, Y)/(n)
  Mt <- t(M)
  svdm <- svd(Mt, nu = 1, nv = 1)

  #-- initalise c = Mu and lambda --#
  uold <- svdm$u;
  vold <- svdm$v;
  convg <- list()

  for(h in 1:ncomp){

    for(iter in 1:100){

      lambda_h <- lambda*max(abs(proj%*%M%*%uold))

      uv <- get_uv(M, uold, vold, lambda_h, proj, abstol = abstol, reltol= reltol,
                   max_itter= max_itter, compositional=compositional)

      u <- uv$u; v <- uv$v

      if(norm(uold - u,"e") || norm(u) < 1e-6 || q == 1) {
        break
      }
    }
    convg$rnorm[h] <- uv$rnorm
    convg$snorm[h] <- uv$snorm
    convg$niter[h] <- uv$niter
    convg$lambda[h] <- lambda_h

    U[,h] <- u
    V[,h] <- v

    #-- compute score --#
    z <- X%*%v

    #-- save away scaled loadings --#
    z_c <- mean(z)
    z <- z - z_c;
    z_s <- drop(sqrt(crossprod(z)))
    Z[,h] <- z <- normalise(z, z_s)
    Vscaled[,h] <- v <- normalise(v, z_s)

    #-- construct the unwanted space --#
    cx <- t(X)%*%Z[,1:h, drop = F]

    #-- construct wanted space --#
    proj <- diag(1,p,p) - get_proj(cx)
    proj_matrices[,,h+1] <- proj

  }

  Beta <- sapply(1:ncomp, function(h) tcrossprod(Vscaled[,1:h, drop = F])%*%M*n, simplify = F)

  res = list(x.scores = Z,
             x.weights = Vscaled,
             y.weights = U,
             Beta = Beta,
             Px = proj_matrices,
             convg = convg,
             V = V)
  return(res)
}


