show_nonzero <- function(X){
  X <- as.matrix(X)
  rownames(X) <- 1:nrow(X)
  return(X[which(rowSums(abs(X)) > 0),] )
}


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


get_uv <- function(M, uold, lambda, proj, abstol = abstol, reltol= reltol,
             max_itter= max_itter, compositional=compositional){
  q <- ncol(M)

  if( q == 1 ) {
    admm_v <- get_v(M, proj, lambda, rho=max(abs(proj%*%M)),
                    abstol = abstol, reltol= reltol,max_itter= max_itter, compositional=compositional)
    res <- admm_v;    res$u <- 1
    return(res)

  } else{

    for (i in 1:500){

      Mu <- M%*%uold
      admm_v <- get_v(Mu, proj, lambda_h, rho=max(abs(proj%*%Mu)),
                      abstol = abstol, reltol= reltol,max_itter= max_itter, compositional=compositional)

      #-- get v direction vector --#
      v <- admm_v$v

      #-- get u direction vector --#
      u <- normalise(Mt%*%v, norm(Mt%*%v,"e"))

      if(norm(u - uold, type = "e") < 1e-10 || norm(u) < 1e-6 ){
        break
      }
      uold <- u
    }
    res <- admm_v;    res$u <- u
    return(res)
  }
}


