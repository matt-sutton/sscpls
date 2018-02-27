
#-- ADMM implementation of get_v --#
get_u <- function(ell, proj, lambda, rho=1, abstol = 1e-04, reltol= 1e-02, max_itter = 10^3, xi_init=NULL, compositional) {

  #-- Some simple internal funcitons --#
  prox_soft <- function(x, lambda){
    return(sign(x)*pmax(0, abs(x) - lambda))
  }

  prox_l2 <- function(u){
    nrm <- drop(sqrt(crossprod(u)))
    if(nrm  < 1 ){
      return(u)
    }  else {
      return(u/nrm)
    }
  }

  #-- Projection operator --#
  if(compositional){
    proj_S <- function(u, proj){
      if(all(u == 0)) {return(u)}
      simplex_root<- function(eta, v){
        res <- proj%*%(v - eta) +10e-10
        return(sum(res))
      }
      eta <- uniroot(simplex_root, interval = c(-max(abs(u)), max(abs(u))), v=u)$root
      return(proj%*%(u - eta))
    }
  } else {
    proj_S <- function(u, proj){
      return(proj%*%u)
    }
  }

  p <- length(ell)

  xi <- if(is.null(xi_init)) rep(0,p) else xi_init
  zold <- z <- rep(0, p)

  ell <- ell/rho; lambda <- lambda/rho

  for( iter in 1:max_itter ){

    u <- prox_l2(prox_soft(z - xi + ell, lambda))

    z <- proj_S(xi + u, proj = proj)

    xi <- xi + u - z

    #-- terminal checks (taken from Boyd matlab example)--#
    history_r_norm <- crossprod(u - z)^0.5;
    history_eps_pri <- sqrt(p)*abstol + reltol*max(norm(u), norm(-z))
    history_s_norm <- crossprod(rho*(round(z, 8) - round(zold, 8)))^0.5;
    history_eps_dual <- sqrt(p)*abstol + reltol*norm(xi)
    zold <- z

    if(history_r_norm < history_eps_pri &&
       history_s_norm < history_eps_dual){
      break
    }
  }
  #cat(iter, history_r_norm < history_eps_pri, history_s_norm < history_eps_dual)
  #cat("\n\n",u,"\n",z,"\n")
  res <- list(u = u, rnorm = history_r_norm, snorm = history_s_norm)
  return( res )
}

