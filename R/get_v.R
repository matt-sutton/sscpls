get_uv <- function(XZ, M, uold, vold, lambda, proj, abstol = abstol, reltol= reltol,
                   max_itter= max_itter, compositional=compositional, keepX){

  Mu <- M%*%uold
  if(!keepX) lambda <- lambda*max(abs(proj%*%Mu))

  #-- get v direction vector --#
  res <- get_v(Mu, proj, lambda, rho=max(abs(proj%*%Mu)),
                  abstol = abstol, reltol= reltol,max_itter= max_itter, compositional=compositional, keepX = keepX)

  ## Continuity correction of Moghaddam
  inds <- which(abs(res$v)>0)
  if(length(inds) > 0){
    proj <- diag(1,length(inds),length(inds)) - get_proj(XZ[inds,, drop = F])$H
    pM <- proj%*%M[inds,]
    #cbind(pM%*%uold/drop(sqrt(crossprod(pM%*%uold))),show_nonzero(res$v))
    res$v[inds] <- normalise(pM%*%uold, norm(pM%*%uold,"e"))
  }

  #-- get u direction vector --#
  res$u <- normalise(t(M)%*%res$v, norm(t(M)%*%res$v,"e"))

  return(res)
}


#-- ADMM implementation of get_v --#
get_v <- function(ell, proj, lambda, rho=1, abstol = 1e-04, reltol= 1e-02, max_itter = 10^3, xi_init=NULL, compositional, keepX) {

  #-- Projection operator --#
  if(compositional){
    proj_S <- proj_compositional
  } else {
    proj_S <- proj_orthogonal
  }

  p <- length(ell)

  z <- xi <- rep(0,p)

  xi <- proj%*%ell/rho
  l <- lambda

  for( iter in 1:max_itter ){

    zold <- z

    x <- proj_S(z - xi + ell/rho, proj = proj)

    if(keepX) l <- sort(abs(x+xi)*rho, decreasing = T)[lambda+1]
    z <- prox_l2(prox_soft(x + xi, l/rho))

    xi <- xi + x - z

    #-- terminal checks (taken from Boyd matlab example)--#
    history_r_norm <- norm(x - z,"f");
    history_eps_pri <- sqrt(p)*abstol + reltol*max(norm(x), norm(-z))
    history_s_norm <- norm(rho*(z - zold));
    history_eps_dual <- sqrt(p)*abstol + reltol*norm(xi)
    zold <- z

    if(history_r_norm < history_eps_pri &&
       history_s_norm < history_eps_dual){
      break
    }

    # Update penalty (Boyd 2010)
    if(history_r_norm > 10*history_s_norm){
      rho <- rho * 2
      xi <- xi/2
    }
    if(history_s_norm > 10*history_r_norm){
      rho <- rho / 2
      xi <- xi*2
    }
  }

  res <- list(v = z, rnorm = history_r_norm, snorm = history_s_norm, niter = iter)
  return( res )
}

