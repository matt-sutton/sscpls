##-- Simrel Package for data simulation --#
library(simrel)
library(ggplot2)
library(ggpubr)
library(sscpls)
library(dplyr)

## --Input: A list of factors for sim study --#
simlist <- list(n =   c(40),
                p =   c(1000, 5000),
                gamma = c(0.4),
                q = c(20),
                m=c(6))


param <- expand.grid(simlist)
nsim <- nrow(param)
nrepeat <- 50

#- Function to get parameter values --#

get_param <- function(paramters){
  return(list(
    n = as.numeric(paramters[1]), p = as.numeric(paramters[2]),
    gamma = as.numeric(paramters[3]), q = as.numeric(paramters[4]),
    m = as.numeric(paramters[5]))
  )
}


#-- Load important packages and R code --#

source('performance_metrics.R')
source('naive_spls.R')
source('cv_naive_spls.R')

##-- Constant factors across sims --#

R2 <- 0.85
ntest <- 10^3
lambda_cv <- seq(0.1, 0.9, by = 0.1)

res <- vector("list", nsim)

library(doSNOW)
cl <- makeCluster(4)
registerDoSNOW(cl)

for(i in 1:nsim){

  cat("\n--------------------------------\n")
  print(param[i,])
  cat("\n--------------------------------\n\n")

  jstart <- 1

  sensitivity <- specificity <- MC <- ME <- MSE <- MSEP <- data.frame( spls_nipals = rep(0,nrepeat),
                                                                       spls_simpls = rep(0,nrepeat),
                                                                       ssc_simpls = rep(0,nrepeat))

  for(j in jstart:nrepeat){
    cat(j,"... ")
    set.seed(j)
    parameters <- get_param(param[i,])
    n <- parameters$n; p <- parameters$p; m <- parameters$m; gamma <- parameters$gamma; q <- parameters$q

    mydata <- simrel(n = n, p = p, m = m, q = q, relpos = 1:m, gamma = gamma, R2 = R2, ntest = ntest)

    #--- Tuning parameters for different methods using 10-fold cross validation -----#
    # naive SPLS-NIPALS
    set.seed(j)
    cvfit1 <- cv_naive_spls(mydata$X, mydata$Y, fold = 10, K = 1:m, lambda = lambda_cv, deflate = "nipals")
    fit1 <- naive_spls(mydata$X, mydata$Y, cvfit1$K, lambda = cvfit1$lambda.1se, deflate = "nipals")

    # naive SPLS-SIMPLS
    set.seed(j)
    cvfit2 <- cv_naive_spls(mydata$X, mydata$Y, fold = 10, K = 1:m, lambda = lambda_cv, deflate = "simpls")
    fit2 <- naive_spls(mydata$X, mydata$Y, cvfit2$K, lambda = cvfit2$lambda.1se, deflate = "simpls")

    # SSC-SIMPLS
    set.seed(j)
    cvfit3 <- cv_sscpls(mydata$X, mydata$Y, fold = 10, K = 1:m, lambda = lambda_cv)
    fit3 <- sscpls(mydata$X, mydata$Y, ncomp = cvfit3$K, lambda = cvfit3$lambda.1se, scale = F)

    #----------------------------#
    #-- Perf Metrics --#
    #
    Bhat <- cbind(fit1$Beta[[cvfit1$K]], fit2$Beta[[cvfit2$K]], fit3$Beta[[cvfit3$K]])

    # Selected variables
    Sel_features <- apply(Bhat, MARGIN = 2, FUN = feature_selection, B_True = mydata$beta)
    sensitivity[j, ] <- Sel_features[1,]
    specificity[j, ] <- Sel_features[2,]
    MC[j, ]  <- Sel_features[3,]

    # MSE estimate performance
    V <- t(mydata$Rotation)%*%mydata$Sigma[-1,-1]%*%mydata$Rotation
    ME[j,] <- apply(Bhat, MARGIN = 2, FUN = Cal_ME, B_True = mydata$beta, V = V)

    # MSEP Prediction performance
    MSEP[j,] <- apply(Bhat, MARGIN = 2, FUN = Cal_MSPE, x=mydata$TESTX, y=mydata$TESTY, oldX=mydata$X)
  }

  df <- data.frame(sensitivity = unlist(sensitivity), MC = unlist(MC), specificity = unlist(specificity), ME = unlist(ME), MSEP = unlist(MSEP), method =
                     factor(rep(c("spls_nipals","spls_simpls", "ssc_simpls"), each = nrepeat)))


  p1 <- ggplot(df, aes(x = method, y = sensitivity)) +
    geom_boxplot()
  p2 <- ggplot(df, aes(x = method, y = specificity)) +
    geom_boxplot()
  p3 <- ggplot(df, aes(x = method, y = ME)) +
    geom_boxplot()
  p4 <- ggplot(df, aes(x = method, y = MSEP)) +
    geom_boxplot()
  fig <- ggarrange(p1,p2,p3,p4)
  fig

  fig_name <- paste0("n_",parameters$n,"p_",parameters$p,"gam_",parameters$gamma,
                     "q_",parameters$q, "m_",m,"R2_",R2,".rda")

  save(fig, file = paste0("results/fig_",fig_name))
  save(df, file = paste0("results/data_",fig_name))

  get_boot_median <- function(x){
    boots <- boot::boot(x, function(x, i){median(x[i])}, R = 500)
    return(data.frame(median=boots$t0, se = sd(boots$t)))
  }

  res[[i]]$sensitivity <-
    df %>%
    group_by(method) %>%
    do(get_boot_median(.$sensitivity))

  res[[i]]$specificity <-
    df %>%
    group_by(method) %>%
    do(get_boot_median(.$specificity))

  res[[i]]$ME <-
    df %>%
    group_by(method) %>%
    do(get_boot_median(.$ME))

  res[[i]]$MSEP <-
    df %>%
    group_by(method) %>%
    do(get_boot_median(.$MSEP))

  res[[i]]$MC <-
    df %>%
    group_by(method) %>%
    do(get_boot_median(.$MC))

  save(res, file = "results/data_summary.rda")
}
stopCluster(cl)
