
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

  z <- xi <- rep(0,p)

  for( iter in 1:max_itter ){

    zold <- z

    x <- proj_S(z - xi + ell/rho, proj = proj)

    z <- prox_l2(prox_soft(x + xi, lambda/rho))

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

  res <- list(u = z, rnorm = history_r_norm, snorm = history_s_norm, niter = iter)
  return( res )
}

