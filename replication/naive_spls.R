
naive_spls <- function(X, Y, ncomp=2, lambda, scale = F, center=T, deflate="simpls") {
  
  #-- Internal funciton for vector normalisation --#
  normalise <- function(x, xnorm){
    xnorm <- drop(xnorm + 1*(xnorm == 0))
    return(x/(xnorm))
  }
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
  
  #-- Data parameters --#
  n <- nrow(X);  p <- ncol(X);   q <- ncol(Y)
  
  #-- Scaling --#
  X <- scale(X,center = center, scale = scale)
  Y <- scale(Y,center = center, scale = scale)
  
  #-- Initialisation --#
  V <- matrix(0, p, ncomp);  U <- matrix(0, q, ncomp)
  Z_s <- Z <- matrix(0, n, ncomp);  Z_center <- rep(0, ncomp); Z_scale <- rep(1, ncomp)
  
  lambdas <- rnorm <- snorm <- NULL
  
  #-- Compute a projection matrix onto the column space of S --#
  proj_matrices <- array(0, dim = c(p,p,ncomp+1))
  proj_matrices[,,1] <- proj <- diag(1, p,p)
  
  #-- Compute the svd of M --#
  M_k <- M <- crossprod(X, Y)
  X_k <- X; Y_k <- Y
  
  for(h in 1:ncomp){
    
    #-- find lambda (Witten style Tuning) --#
    
    #-- initalise c = Mu and lambda --#
    svdm <- svd(t(M_k), nu = 1, nv = 1)
    uold <- svdm$u;
    
    lambda_max <- max(abs(M_k%*%uold))
    lambda_h <- lambda_max*lambda
    lambdas <- c(lambdas, lambda_h)
    
    for (i in 1:500){
      
      #-- get v direction vector --#
      v <- prox_soft(M_k%*%uold, lambda_h)
      v <- if(!all(v == 0)) v/norm(v,"e")
      
      #-- get u direction vector --#
      u <- prox_l2(t(M_k)%*%v)
      u <- if(!all(u == 0)) u/norm(u,"e")
      
      if(norm(u - uold, type = "e") < 1e-10){
        break
      }
      uold <- u
    }
    U[,h] <- u
    
    #-- compute score --#
    Z[,h] <- z <- X_k%*%v
    
    #-- save away scaled loadings --#
    z_c <- mean(z)
    z <- z - z_c;
    z_s <- drop(sqrt(crossprod(z)))
    
    #-- Deflate methods NIPALS or SIMPLS --#
    if(deflate == "simpls"){
      V[,h] <- normalise(v, z_s)
      Z_s[,h] <- normalise(z, z_s)
      R <- t(X)%*%Z_s[,1:h, drop = F]
      M_k <- (diag(1,p,p) - get_proj(R)) %*% M
    }
    
    if(deflate == "nipals"){
      Px <- diag(1,n,n) - get_proj(Z)
      X_k <- Px%*%X; Y_k <- Px%*%Y
      M_k <- crossprod(X_k, Y_k)
      V[,h] <- v
      Z_s[,h] <- normalise(z, z_s^2)
    }
    
  }
  if(deflate == "simpls"){
    Beta <- sapply(1:ncomp, function(h) V[,1:h, drop = F]%*%MASS::ginv(crossprod(Z_s[,1:h, drop = F]))%*%t(V[,1:h, drop = F])%*%M, simplify = F)
  }
  if(deflate == "nipals"){
    R <- crossprod(X, Z_s[,1:h, drop = F])
    V_star <- V%*%MASS::ginv(t(R)%*%V)
    Beta <- sapply(1:ncomp, function(h) V_star[,1:h, drop = F]%*%MASS::ginv(crossprod(Z[,1:h, drop = F]))%*%t(V_star[,1:h, drop = F])%*%M, simplify = F)
  }
  
  res = list(x.scores = Z,
             x.weights = V,
             x.loads = t(X)%*%Z_s[,1:h, drop = F],
             y.weights = U,
             Beta = Beta,
             Px = proj_matrices,
             lambdas = lambdas)
  return(res)
}