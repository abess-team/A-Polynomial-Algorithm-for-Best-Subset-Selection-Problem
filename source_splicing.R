# source("utilities.R")
library(glmnet)
library(ncvreg)
library(L0Learn)
library(BeSS)
library(abess)
library(mvtnorm)
library(EvaluationMeasures)

repmat = function(X, m, n) {
  ## R equivalent of repmat (matlab)
  X = as.matrix(X)
  mx = dim(X)[1]
  nx = dim(X)[2]
  matrix(t(matrix(X, mx, nx * n)), mx * m, nx * n, byrow = T)
}

getdf = function(coef.beta) {
  apply(abs(coef.beta) > 1e-10, 2, sum)
}

loglik <- function(X, y, beta, family) {
  K = dim(beta)[2]
  link = cbind(1, X) %*% beta
  yrep = repmat(y, 1, K)
  if (family == "gaussian")
    return(apply((yrep - link)^2, 2, sum))
  if (family == "poisson")
    return(apply(exp(link) - yrep * link, 2, sum))
  if (family == "binomial")
    return(apply(log(1 + exp(link)) - yrep * link, 2, sum))
}

classification_measure <- function(actual, predict) {
  actual <- as.vector(actual)
  predict <- as.vector(predict)
  tpr <- EvaluationMeasures.TPR(Real = actual, Predicted = predict)
  fpr <- EvaluationMeasures.FPR(Real = actual, Predicted = predict)
  mcc <- EvaluationMeasures.MCC(Real = actual, Predicted = predict)
  c(tpr, fpr, mcc)
}

beta_error_measure <- function(actual, predict) {
  actual <- as.vector(actual)
  predict <- as.vector(predict)
  absolute_error <- sum((actual - predict)^2)
  relative_error <- absolute_error / sum(actual^2)
  relative_error
}

gendata = function(n = 100, p = 100, Nonzero, SNR = 100, neg = TRUE){
  x = matrix(rnorm(n*p), n, p)
  Tbeta = rep(0,p)
  if(neg){
    Tbeta[Nonzero] = c(rep(1, round(length(Nonzero)/2)),
                       rep(-1, length(Nonzero) - round(length(Nonzero)/2)))
  } else {
    Tbeta[Nonzero] = 1
  }
  epsilon = rnorm(n)
  Nnorm = drop(t(Tbeta)%*%Tbeta)/SNR
  epsilon = epsilon*sqrt(Nnorm)
  y = drop(x%*%Tbeta + epsilon)
  return(list(x=x, y=y, Tbeta=Tbeta))
}

# type 2 and 3 correlated
gendata2 = function(n = 100, p = 100, Nonzero, SNR = 100, sigma = 1, neg = TRUE, rho = 0.5){
  cor_struct = min(p, 1000)
  if(sigma == 1){
    Sigma = matrix(rep(rho, cor_struct * cor_struct), cor_struct, cor_struct)
    diag(Sigma) = 1
  }else{
    Sigma = rho^abs(matrix(1:cor_struct, cor_struct, cor_struct)
                    - matrix(1:cor_struct, cor_struct, cor_struct, byrow = TRUE))
  }
  X_1 = rmvnorm(n, rep(0, cor_struct),sigma = Sigma)
  if(p<=1000){
    x = X_1
  } else {
    X_2 = matrix(rnorm(n * (p-cor_struct)), n, p-cor_struct)
    x = cbind(X_1[, 1:(floor(cor_struct/2))],
              X_2, X_1[, (floor(cor_struct/2)+1):cor_struct])
  }
  Tbeta = rep(0,p)
  if(neg){
    Tbeta[Nonzero] = c(rep(1, round(length(Nonzero)/2)),
                       rep(-1, length(Nonzero) - round(length(Nonzero)/2)))
  } else {
    Tbeta[Nonzero] = 1
  }
  epsilon = rnorm(n)
  x = scale(x)*sqrt(n)/sqrt(n-1)
  if(p<=1000){
    tmp = length(Nonzero)
  } else {
    tmp = drop(t(Tbeta[1:500])%*%Sigma[1:500, 1:500]%*%Tbeta[1:500]+
                 sum(Tbeta[501:(p-500)]^2)
               + t(Tbeta[(p-499):p]%*%Sigma[501:1000, 501:1000]%*%Tbeta[(p-499):p]))
  }
  Nnorm = tmp/SNR
  epsilon = epsilon*sqrt(Nnorm)
  y = drop(x%*%Tbeta + epsilon)
  return(list(x=x, y=y, Tbeta=Tbeta))
}

# decay
gendata2_1 <- function(n = 100, p = 100, Nonzero, SNR = 100, sigma = 1, neg = TRUE, rho = 0.5){
  Sigma = rho^abs(matrix(1:p, p, p) - matrix(1:p, p, p, byrow = TRUE))
  diag(Sigma) = 1
  x = rmvnorm(n, sigma = Sigma)
  Sigma[Sigma < 1e-10] <- 0
  Tbeta = rep(0,p)
  if(neg){
    Tbeta[Nonzero] = c(rep(1, round(length(Nonzero)/2)),
                       rep(-1, length(Nonzero) - round(length(Nonzero)/2)))
  } else {
    Tbeta[Nonzero] = 1
  }
  epsilon = rnorm(n)
  x = scale(x)*sqrt(n)/sqrt(n-1)
  snr_numerator <- t(Tbeta) %*% Sigma %*% as.matrix(Tbeta)
  snr_numerator <- snr_numerator[1, 1]
  Nnorm = snr_numerator / SNR
  epsilon = epsilon*sqrt(Nnorm)
  y = drop(x%*%Tbeta + epsilon)
  return(list(x=x, y=y, Tbeta=Tbeta))
}

# constant
gendata2_2 <- function(n = 100, p = 100, Nonzero, SNR = 100, sigma = 1, neg = TRUE, rho = 0.5){
  Sigma = matrix(rep(rho, p * p), p, p)
  diag(Sigma) = 1
  x = rmvnorm(n, sigma = Sigma)

  Tbeta = rep(0,p)
  if(neg){
    Tbeta[Nonzero] = c(rep(1, round(length(Nonzero)/2)),
                       rep(-1, length(Nonzero) - round(length(Nonzero)/2)))
  } else {
    Tbeta[Nonzero] = 1
  }
  epsilon = rnorm(n)
  x = scale(x)*sqrt(n)/sqrt(n-1)
  snr_numerator <- t(Tbeta) %*% Sigma %*% as.matrix(Tbeta)
  snr_numerator <- snr_numerator[1, 1]
  Nnorm = snr_numerator / SNR
  epsilon = epsilon*sqrt(Nnorm)
  y = drop(x%*%Tbeta + epsilon)
  return(list(x=x, y=y, Tbeta=Tbeta))
}

gendata_low_1 <- function(n = 40, p = 8, sigma = 3, rho = 0.5) {
  Sigma = rho^abs(matrix(1:p, p, p) - matrix(1:p, p, p, byrow = TRUE))
  diag(Sigma) = 1
  x = rmvnorm(n, sigma = Sigma)
  Tbeta = c(3, 1.5, 0, 0, 2, 0, 0, 0)
  if (p <= 8) {
    Tbeta <- Tbeta[1:p]
  } else {
    Tbeta <- c(Tbeta, rep(0, p - 8))
  }
  epsilon = rnorm(n) * sigma
  y = drop(x %*% as.matrix(Tbeta) + epsilon)
  return(list(x=x, y=y, Tbeta=Tbeta))
}

beta_generator <- function() {
  strong_num <- 3
  moderate_num <- 4
  weak_num <- 3
  signal_num <- strong_num + moderate_num + weak_num

  strong_signal <- rnorm(n = strong_num, sd = 10)
  moderate_signal <- rnorm(n = moderate_num, sd = 5)
  weak_signal <- rnorm(n = weak_num, sd = 2)
  beta_value <- c(strong_signal, moderate_signal, weak_signal)

  beta_value <- beta_value[sample(1:signal_num, size = signal_num, replace = FALSE)]
  beta_value
}

sparse_beta_generator <- function(p, Nonzero) {
  Tbeta = rep(0, p)
  beta_value <- beta_generator()
  Tbeta[Nonzero] <- beta_value
  Tbeta
}

# independence
gendata_high_1 <- function(n = 100, p = 100, Nonzero, SNR = 1){
  x = matrix(rnorm(n*p), n, p)
  Tbeta <- sparse_beta_generator(p, Nonzero)

  epsilon = rnorm(n)
  Nnorm = drop(t(Tbeta) %*% Tbeta)/SNR
  y = drop(x%*%Tbeta + epsilon)
  return(list(x=x, y=y, Tbeta=Tbeta))
}

# constant
gendata_high_2 <- function(n = 100, p = 100, Nonzero, SNR = 1, rho){
  Sigma = matrix(rep(rho, p * p), p, p)
  diag(Sigma) = 1

  x = rmvnorm(n, sigma = Sigma)
  Tbeta <- sparse_beta_generator(p, Nonzero)

  Sigma[Sigma < 1e-10] <- 0
  x = scale(x)*sqrt(n)/sqrt(n-1)
  snr_numerator <- t(Tbeta) %*% Sigma %*% as.matrix(Tbeta)
  snr_numerator <- snr_numerator[1, 1]
  Nnorm = snr_numerator / SNR

  epsilon = rnorm(n)
  y = drop(x%*%Tbeta + epsilon)
  return(list(x=x, y=y, Tbeta=Tbeta))
}

# L0Learn
sim.l0learn = function(x, y, Tbeta, npara = 20, penalty = "L0", algorithm = "CD")
{
  p = ncol(x)
  Nonzero = which(Tbeta!=0)

  time.L0CDSIC = system.time(fit.L0CDSIC <- L0Learn.fit(x, y, penalty = penalty,
                                                          algorithm = algorithm, nLambda = npara,
                                                          screenSize = p))
  fit.L0CDSIC$loss = "SquareError"
  df = unlist(fit.L0CDSIC$suppSize)
  time_2 <- system.time(pred <- predict(fit.L0CDSIC, newx = x, lambda = unlist(fit.L0CDSIC$lambda)))
  time_3 <- system.time(MSE <- apply((pred - y)^2, 2, mean))
  time_4 <- system.time(SIC <- nrow(x)*log(MSE) + log(p)*log(log(nrow(x)))*df)
  time.L0CDSIC <- time.L0CDSIC + time_2 + time_3 + time_4
  time.L0CDSIC <- time.L0CDSIC[3]

  beta.L0CDSIC_temp = coef(fit.L0CDSIC, lambda = unlist(fit.L0CDSIC$lambda)[which.min(SIC)])
  beta.L0CDSIC = rep(0, p)
  beta.L0CDSIC[beta.L0CDSIC_temp@i[-1]] = beta.L0CDSIC_temp@x[-1]

  PE.L0CDSIC = sum((beta.L0CDSIC - Tbeta)^2)/sum(Tbeta^2)
  TP.L0CDSIC = length(intersect(Nonzero, which(beta.L0CDSIC!=0)))
  ms.L0CDSIC = length(which(beta.L0CDSIC!=0))
  FP.L0CDSIC = ms.L0CDSIC - TP.L0CDSIC

  Nonzero <- rep(0, ncol(x))
  Nonzero[Tbeta != 0] <- 1
  beta_index <- ifelse(abs(beta.L0CDSIC) > 1e-6, 1, 0)
  class_res <- classification_measure(Nonzero, beta_index)
  k <- sum(beta_index)

  return(c(time.L0CDSIC = time.L0CDSIC,
           TP.L0CDSIC = TP.L0CDSIC, FP.L0CDSIC = FP.L0CDSIC,
           PE.L0CDSIC = PE.L0CDSIC,
           class_res, k))
}

sim.l0learn.cv = function(x, y, Tbeta, npara = 20, penalty = "L0", algorithm = "CD")
{
  p = ncol(x)
  Nonzero = which(Tbeta!=0)
  if(algorithm == "CD"){
    time.L0CDcv = system.time(fit.L0CDSIC <- L0Learn.cvfit(x, y, penalty = penalty, algorithm = algorithm, nLambda = npara, screenSize = p))[3]
    cv_error <- unlist(fit.L0CDSIC$cvMeans)
    cv_lambda <- fit.L0CDSIC$fit$lambda[[1]]
    beta.L0CDcv <- coef(fit.L0CDSIC, lambda = cv_lambda[which.min(cv_error)])
    beta.L0CDcv <- as.vector(beta.L0CDcv)[-1]

    PE.L0CDcv = sum((beta.L0CDcv - Tbeta)^2)/sum(Tbeta^2)
    TP.L0CDcv = length(intersect(Nonzero, which(beta.L0CDcv!=0)))
    ms.L0CDcv = length(which(beta.L0CDcv!=0))
    FP.L0CDcv = ms.L0CDcv - TP.L0CDcv

    Nonzero <- rep(0, ncol(x))
    Nonzero[Tbeta != 0] <- 1
    beta_index <- ifelse(abs(beta.L0CDcv) > 1e-6, 1, 0)
    class_res <- classification_measure(Nonzero, beta_index)
    k <- sum(beta_index)

    return(c(time.L0CDcv = time.L0CDcv, TP.L0CDcv = TP.L0CDcv, FP.L0CDcv = FP.L0CDcv, PE.L0CDcv = PE.L0CDcv, class_res, k))
  }
  else{
    time.L0CDPSIcv = system.time(fit.L0CDPSISIC <- L0Learn.cvfit(x, y, penalty = penalty, algorithm = algorithm, nLambda = npara, screenSize = p))[3]
    cv_error <- unlist(fit.L0CDPSISIC$cvMeans)
    cv_lambda <- fit.L0CDPSISIC$fit$lambda[[1]]
    beta.L0CDPSIcv <- coef(fit.L0CDPSISIC, lambda = cv_lambda[which.min(cv_error)])
    beta.L0CDPSIcv <- as.vector(beta.L0CDPSIcv)[-1]

    PE.L0CDPSIcv = sum((beta.L0CDPSIcv - Tbeta)^2)/sum(Tbeta^2)
    TP.L0CDPSIcv = length(intersect(Nonzero, which(beta.L0CDPSIcv!=0)))
    ms.L0CDPSIcv = length(which(beta.L0CDPSIcv!=0))
    FP.L0CDPSIcv = ms.L0CDPSIcv - TP.L0CDPSIcv

    Nonzero <- rep(0, ncol(x))
    Nonzero[Tbeta != 0] <- 1
    beta_index <- ifelse(abs(beta.L0CDPSIcv) > 1e-6, 1, 0)
    class_res <- classification_measure(Nonzero, beta_index)
    k <- sum(beta_index)

    return(c(time.L0CDPSIcv = time.L0CDPSIcv, TP.L0CDPSIcv = TP.L0CDPSIcv, FP.L0CDPSIcv = FP.L0CDPSIcv, PE.L0CDPSIcv = PE.L0CDPSIcv, class_res, k))
  }
}

# BeSS
# if flag = TRUE, use abess; otherwise, sdar
sim.bess = function(x, y, Tbeta, s.max = 20, warm.start = FALSE, exchange.num = 5, flag = TRUE, criterion = "cv")
{
  Nonzero = which(Tbeta!=0)
  if(flag){
    if (low_dimension_run_time) {
      time.abess = mean(microbenchmark(fit.abess <- abess(x, y, tune.path = "sequence", c.max = exchange.num,
                                                          support.size = 1:s.max, warm.start = warm.start),
                                       times = 20)[["time"]], trim = 0.05) / 10^9
    } else {
      time.abess = microbenchmark(fit.abess <- abess(x, y, tune.path = "sequence", c.max = exchange.num,
                                                     support.size = 1:s.max, warm.start = warm.start),
                                  times = 1)[["time"]] / 10^9
    }
    beta.abess = coef(fit.abess, support.size = fit.abess$best.size)[-1]
    PE.abess = sum((beta.abess - Tbeta)^2)/sum(Tbeta^2)
    TP.abess = length(intersect(Nonzero, which(beta.abess!=0)))
    ms.abess = length(which(beta.abess!=0))
    FP.abess = ms.abess - TP.abess

    Nonzero <- rep(0, ncol(x))
    Nonzero[Tbeta != 0] <- 1
    beta_index <- ifelse(abs(beta.abess) > 1e-6, 1, 0)
    class_res <- classification_measure(Nonzero, beta_index)
    k <- sum(beta_index)

    return(c(time.abess = time.abess, TP.abess = TP.abess, FP.abess = FP.abess, PE.abess = PE.abess, class_res, k))
  } else {
    time.sdar = system.time(fit.sdar <- BeSS::bess(x, y, method = "sequential", s.list = 1:s.max))[3]
    num <- nrow(x)
    p <- ncol(x)
    if (criterion == "sic") {
      sic <- log(fit.sdar[["mse"]]) * nrow(x) + log(log(nrow(x))) * fit.sdar[["s.list"]] * log(dim(x)[2])
      min_criterion <- which.min(sic)
    } else {
      if (num > p) {
        criterion <- "AIC"
      } else {
        criterion <- "EBIC"
      }
      min_criterion = which.min(fit.sdar[[criterion]])
    }
    beta.sdar = fit.sdar$beta[, min_criterion]

    PE.sdar = sum((beta.sdar - Tbeta)^2)/sum(Tbeta^2)
    TP.sdar = length(intersect(Nonzero, which(beta.sdar!=0)))
    ms.sdar = length(which(beta.sdar!=0))
    FP.sdar = ms.sdar - TP.sdar

    Nonzero <- rep(0, ncol(x))
    Nonzero[Tbeta != 0] <- 1
    beta_index <- ifelse(abs(beta.sdar) > 1e-6, 1, 0)
    class_res <- classification_measure(Nonzero, beta_index)
    k <- sum(beta_index)

    return(c(time.sdar = time.sdar, TP.sdar = TP.sdar, FP.sdar = FP.sdar, PE.sdar = PE.sdar, class_res, k))
  }
}

sim.glmnet = function(x, y, Tbeta, npara = 20, one.se = TRUE, criterion = "cv")
{
  Nonzero = which(Tbeta!=0)

  if (criterion == "cv") {
    time.L1cv = system.time(fit.L1cv <- cv.glmnet(x, y, nlambda = npara))[3]
    if (one.se) {
      markcv = which(fit.L1cv$lambda.1se == fit.L1cv$lambda) + 1
    } else {
      markcv = which(fit.L1cv$lambda.min == fit.L1cv$lambda)
    }
    beta.L1cv = fit.L1cv$glmnet.fit$beta[,markcv]
  } else {
    time_1 <- system.time(fit.L1cv <- glmnet(x, y, nlambda = npara))
    time_2 <- system.time(dev <- deviance(fit.L1cv)) ## dev is MSE
    time_3 <- system.time(reg.df <- fit.L1cv$df)
    if (criterion == "ebic") {
     time_4 <- system.time(ic <- nrow(x) * log(dev / nrow(x)) + log(nrow(x)) * reg.df + log(choose(dim(x)[2], reg.df)))
    } else if (criterion == "sic") {
      time_4 <- system.time(ic <- log(dev / nrow(x)) * nrow(x) + log(log(nrow(x))) * reg.df * log(dim(x)[2]))
    }  else if (criterion == "bic") {
      time_4 <- system.time(ic <- nrow(x) * log(dev / nrow(x)) + log(nrow(x)) * reg.df)
    }
    time.L1cv <- time_1 + time_2 + time_3 + time_4
    time.L1cv <- time.L1cv[3]
    markcv = which.min(ic)
    beta.L1cv = fit.L1cv$beta[, markcv]
  }

  PE.L1cv = sum((beta.L1cv - Tbeta)^2)/sum(Tbeta^2)
  TP.L1cv = length(intersect(Nonzero, which(beta.L1cv!=0)))
  ms.L1cv = length(which(beta.L1cv!=0))
  FP.L1cv = ms.L1cv - TP.L1cv

  Nonzero <- rep(0, ncol(x))
  Nonzero[Tbeta != 0] <- 1
  beta_index <- ifelse(abs(beta.L1cv) > 1e-6, 1, 0)
  class_res <- classification_measure(Nonzero, beta_index)
  k <- sum(beta_index)

  return(c(time.L1cv = time.L1cv, TP.L1cv = TP.L1cv, FP.L1cv = FP.L1cv, PE.L1cv = PE.L1cv, class_res, k))
}

sim.scad <- function(x, y, Tbeta, npara = 20, criterion = "cv") {
  Nonzero = which(Tbeta!=0)

  if (criterion == "cv") {
    time.L1cv = system.time(fit.L1cv <- cv.ncvreg(x, y, penalty = "SCAD",
                                                  nlambda = npara, warn = FALSE))[3]
    beta <- coef(fit.L1cv)
    beta.L1cv <- beta[-1]
  } else {
    time_1 <- system.time(fit.L1cv <- ncvreg(x, y, nlambda = npara, penalty = "SCAD", warn = FALSE))
    time_2 <- system.time(dev <- loglik(x, y, fit.L1cv[["beta"]], family = "gaussian"))
    time_3 <- system.time(reg.df <- getdf(fit.L1cv[["beta"]][-1, , drop = FALSE]))
    if (criterion == "ebic") {
      time_4 <- system.time(ic <- nrow(x) * log(dev / nrow(x)) + log(nrow(x)) * reg.df + log(choose(dim(x)[2], reg.df)))
    } else if (criterion == "sic") {
      time_4 <- system.time(ic <- log(dev / nrow(x)) * nrow(x) + log(log(nrow(x))) * reg.df * log(dim(x)[2]))
    } else if (criterion == "bic") {
      time_4 <- system.time(ic <- nrow(x) * log(dev / nrow(x)) + log(nrow(x)) * reg.df)
    }
    time.L1cv <- time_1 + time_2 + time_3 + time_4
    time.L1cv <- time.L1cv[3]
    markcv = which.min(ic)
    beta.L1cv = fit.L1cv$beta[, markcv][-1]
  }

  PE.L1cv = sum((beta.L1cv - Tbeta)^2)/sum(Tbeta^2)
  index <- ifelse(abs(beta.L1cv) > 1e-10, 1, 0)
  TP.L1cv = length(intersect(Nonzero, which(index!=0)))
  ms.L1cv = length(which(index!=0))
  FP.L1cv = ms.L1cv - TP.L1cv

  Nonzero <- rep(0, ncol(x))
  Nonzero[Tbeta != 0] <- 1
  beta_index <- ifelse(abs(beta.L1cv) > 1e-6, 1, 0)
  class_res <- classification_measure(Nonzero, beta_index)
  k <- sum(beta_index)

  return(c(time.L1cv = time.L1cv, TP.L1cv = TP.L1cv, FP.L1cv = FP.L1cv, PE.L1cv = PE.L1cv, class_res, k))
}

sim.mcp <- function(x, y, Tbeta, npara = 20, criterion = "cv") {
  Nonzero = which(Tbeta!=0)

  if (criterion == "cv") {
    time.L1cv = system.time(fit.L1cv <- cv.ncvreg(x, y, penalty = "MCP",
                                                  nlambda = npara, warn = FALSE))[3]
    beta <- coef(fit.L1cv)
    beta.L1cv <- beta[-1]
  } else {
    time_1 <- system.time(fit.L1cv <- ncvreg(x, y, nlambda = npara, penalty = "MCP", warn = FALSE))
    time_2 <- system.time(dev <- loglik(x, y, fit.L1cv[["beta"]], family = "gaussian"))
    time_3 <- system.time(reg.df <- getdf(fit.L1cv[["beta"]][-1, , drop = FALSE]))
    if (criterion == "ebic") {
      time_4 <- system.time(ic <- nrow(x) * log(dev / nrow(x)) + log(nrow(x)) * reg.df + log(choose(dim(x)[2], reg.df)))
    } else if (criterion == "sic") {
      time_4 <- system.time(ic <- log(dev / nrow(x)) * nrow(x) + log(log(nrow(x))) * reg.df * log(dim(x)[2]))
    } else if (criterion == "bic") {
      time_4 <- system.time(ic <- nrow(x) * log(dev / nrow(x)) + log(nrow(x)) * reg.df)
    }
    time.L1cv <- time_1 + time_2 + time_3 + time_4
    time.L1cv <- time.L1cv[3]
    markcv = which.min(ic)
    beta.L1cv = fit.L1cv$beta[, markcv][-1]
  }

  PE.L1cv = sum((beta.L1cv - Tbeta)^2)/sum(Tbeta^2)
  index <- ifelse(abs(beta.L1cv) > 1e-10, 1, 0)
  TP.L1cv = length(intersect(Nonzero, which(index!=0)))
  ms.L1cv = length(which(index!=0))
  FP.L1cv = ms.L1cv - TP.L1cv

  Nonzero <- rep(0, ncol(x))
  Nonzero[Tbeta != 0] <- 1
  beta_index <- ifelse(abs(beta.L1cv) > 1e-6, 1, 0)
  class_res <- classification_measure(Nonzero, beta_index)
  k <- sum(beta_index)

  return(c(time.L1cv = time.L1cv, TP.L1cv = TP.L1cv, FP.L1cv = FP.L1cv, PE.L1cv = PE.L1cv, class_res, k))
}

sim.exhaust <- function(x, y, Tbeta, criterion = "Cp") {
  Nonzero <- which(Tbeta!=0)

  time.L1cv <- microbenchmark(cv_fit <- summary(leaps::regsubsets(x = x, y = y, intercept = FALSE,
                                                                  method = "exhaustive",
                                                                  nvmax = length(Tbeta))), times = 1)[["time"]] / 10^9
  cv_fit[["Cp"]] <- cv_fit[["cp"]]
  if (criterion != "Cp") {
    time_2 <- system.time(mse <- apply(cv_fit[["which"]], 1, function(index) {
      x_fit <- x[, index]
      dat <- cbind.data.frame("x" = x_fit, "y" = y)
      mse_value <- mean((lm(y ~ ., data = dat)[["residuals"]])^2)
      mse_value
    }))[3]
    if (criterion == "sic") {
      time_3 <- system.time(cv_fit[["Cp"]] <- log(mse) * nrow(x) +
                              log(log(nrow(x))) * rowSums(cv_fit[["which"]]) * log(dim(x)[2]))[3]
    } else {
      time_3 <- system.time(cv_fit[["Cp"]] <- log(mse) * nrow(x) +
                              log(nrow(x)) * rowSums(cv_fit[["which"]]))[3]
    }
    time.L1cv <- time.L1cv + time_2 + time_3
  }
  time_4 <- microbenchmark(min_gic <- which.min(cv_fit$Cp), times = 1)[["time"]] / 10^9
  time.L1cv <- time.L1cv + time_4

  index <- ifelse(cv_fit$which[min_gic, , drop = TRUE], 1, 0)

  TP.L1cv = length(intersect(Nonzero, which(index != 0)))
  ms.L1cv = length(which(index!=0))
  FP.L1cv = ms.L1cv - TP.L1cv

  x_fit <- x[, which(index != 0)]
  dat <- cbind.data.frame("x" = x_fit, "y" = y)
  beta.L1cv <- lm(y ~ ., data = dat)[["coefficients"]][-1]
  PE.L1cv = sum((beta.L1cv - Tbeta[which(index != 0)])^2)/sum(Tbeta^2)

  Nonzero <- as.numeric(Tbeta != 0)
  beta_index <- index
  class_res <- classification_measure(Nonzero, beta_index)
  k <- sum(beta_index)

  return(c(time.L1cv = time.L1cv, TP.L1cv = TP.L1cv, FP.L1cv = FP.L1cv, PE.L1cv = PE.L1cv, class_res, k))
}

sim_once = function(seed, SNR = 1, type = 1, n, p, k = 5, sigma = 1, npara = 20, neg = TRUE,
                    method = c("sdar", "abess", "L0cv", "L0PSIcv", "L1", "scad", "mcp"),
                    warm.start = TRUE, rho = 0.5, exchange.num = 5){
  set.seed(seed)
  if(type == 1){
    Nonzero = 1:k
    data = gendata(n, p, Nonzero, SNR, neg)
    validata = gendata(n/2, p, Nonzero, SNR, neg)
  }  else {
    if(p<=1000){
      Nonzero = sample(1:p, k)
    } else {
      Nonzero = sample(c(1:500, (p-499):p), k)
    }
    data = gendata2(n, p, Nonzero, SNR, sigma, neg, rho)
  }
  x = data$x
  y = data$y
  Tbeta = data$Tbeta

  default_res <- rep(0, 4)

  #bess
  if("sdar" %in% method){
    res.sdar = sim.bess(x, y, Tbeta, npara, warm.start = warm.start, flag = FALSE)
  } else {
    res.sdar = default_res
  }

  #abess
  if("abess" %in% method){
    res.abess = sim.bess(x, y, Tbeta, npara, warm.start = warm.start, exchange.num = exchange.num, flag= TRUE)
  } else {
    res.abess = default_res
  }

  #L0Learn L0
  if("L0" %in% method){
    res.L0.CD = sim.l0learn(x, y, Tbeta, npara, algorithm = "CD")
  } else {
    res.L0.CD = default_res
  }

  if("L0cv" %in% method){
    res.L0.CD.cv = sim.l0learn.cv(x, y, Tbeta, npara, algorithm = "CD")
  } else {
    res.L0.CD.cv = default_res
  }
  #L0Learn L0PSI
  if("L0PSI" %in% method){
    res.L0.PSI = sim.l0learn(x, y, Tbeta, npara, algorithm = "CDPSI")
  } else {
    res.L0.PSI = default_res
  }

  if("L0PSIcv" %in% method){
    res.L0.PSI.cv = sim.l0learn.cv(x, y, Tbeta, npara, algorithm = "CDPSI")
  } else {
    res.L0.PSI.cv = default_res
  }
  #glmnet L1
  if("L1" %in% method){
    res.L1= sim.glmnet(x, y, Tbeta, npara)
  } else {
    res.L1 = default_res
  }

  # SCAD
  if ("scad" %in% method) {
    res.scad= sim.scad(x, y, Tbeta, npara)
  } else {
    res.scad = default_res
  }

  # SCAD
  if ("mcp" %in% method) {
    res.mcp <- sim.mcp(x, y, Tbeta, npara)
  } else {
    res.mcp <- default_res
  }

  res <- rbind(res.sdar, res.abess, res.L0.CD, res.L0.PSI,
               res.L0.CD.cv, res.L0.PSI.cv, res.L1, res.scad, res.mcp)
  res
}

## low-dimensional comparison
sim_once_low <- function(seed, SNR = 1, type = 1, n, p, k = 5, sigma = 1, npara = 20, neg = TRUE,
                         method = c("exhaust", "exhaust-sic",
                                    "sdar", "sdar-sic", "abess",
                                    "L0cv", "L0-sic", "L0PSIcv", "L0PSI-sic",
                                    "L1", "L1-sic", "scad", "scad-sic", "mcp", "mcp-sic"),
                         warm.start = TRUE, rho = 0.5, exchange.num = 5) {
  set.seed(seed)
  if (type == 1) {
    data <- gendata_low_1(n = n, sigma = sigma, p = p)
  } else if (type == 2) {
    data <- gendata_low_1(n = n, sigma = sigma)
  } else {
    data <- gendata_low_1(n = n, sigma = sigma)
  }
  x = data$x
  y = data$y
  Tbeta = data$Tbeta
  npara <- pmin(npara, ncol(x))

  default_res <- rep(0, 4)

  # oracle
  if ("exhaust" %in% method) {
    res.exhaust = sim.exhaust(x, y, Tbeta)
  } else {
    res.exhaust <- default_res
  }

  if ("exhaust-sic" %in% method) {
    res.exhaust.sic = sim.exhaust(x, y, Tbeta, criterion = "sic")
  } else {
    res.exhaust.sic <- default_res
  }

  if ("exhaust-bic" %in% method) {
    res.exhaust.bic = sim.exhaust(x, y, Tbeta, criterion = "bic")
  } else {
    res.exhaust.bic <- default_res
  }

  #bess
  if("sdar" %in% method){
    res.sdar = sim.bess(x, y, Tbeta, npara, warm.start = warm.start, flag = FALSE)
  } else {
    res.sdar = default_res
  }

  if("sdar-sic" %in% method){
    res.sdar.sic = sim.bess(x, y, Tbeta, npara, warm.start = warm.start, flag = FALSE, criterion = "sic")
  } else {
    res.sdar.sic = default_res
  }

  #abess
  if("abess" %in% method){
    res.abess = sim.bess(x, y, Tbeta, npara, warm.start = warm.start, exchange.num = exchange.num, flag= TRUE)
  } else {
    res.abess = default_res
  }

  #L0Learn L0
  if("L0cv" %in% method){
    res.L0.CD.cv = sim.l0learn.cv(x, y, Tbeta, npara, algorithm = "CD")
  } else {
    res.L0.CD.cv = default_res
  }

  if("L0-sic" %in% method){
    res.L0.CD.cv.sic = sim.l0learn(x, y, Tbeta, npara, algorithm = "CD")
  } else {
    res.L0.CD.cv.sic = default_res
  }

  if("L0PSIcv" %in% method){
    res.L0.PSI.cv = sim.l0learn.cv(x, y, Tbeta, npara, algorithm = "CDPSI")
  } else {
    res.L0.PSI.cv = default_res
  }

  if("L0PSI-sic" %in% method){
    res.L0.PSI.cv.sic = sim.l0learn(x, y, Tbeta, npara, algorithm = "CDPSI")
  } else {
    res.L0.PSI.cv.sic = default_res
  }

  #glmnet L1
  if("L1" %in% method){
    res.L1= sim.glmnet(x, y, Tbeta, npara, FALSE)
  } else {
    res.L1 = default_res
  }

  if("L1-sic" %in% method){
    res.L1.sic= sim.glmnet(x, y, Tbeta, npara, FALSE, criterion = "sic")
  } else {
    res.L1.sic = default_res
  }

  # SCAD
  if ("scad" %in% method) {
    res.scad= sim.scad(x, y, Tbeta, npara)
  } else {
    res.scad = default_res
  }

  if ("scad-sic" %in% method) {
    res.scad.sic = sim.scad(x, y, Tbeta, npara, criterion = "sic")
  } else {
    res.scad.sic = default_res
  }

  # MCP
  if ("mcp" %in% method) {
    res.mcp <- sim.mcp(x, y, Tbeta, npara)
  } else {
    res.mcp <- default_res
  }

  if ("mcp" %in% method) {
    res.mcp.sic <- sim.mcp(x, y, Tbeta, npara, criterion = "sic")
  } else {
    res.mcp.sic <- default_res
  }

  res <- rbind(res.abess,
               res.exhaust,
               res.exhaust.sic,
               res.exhaust.bic,
               res.sdar,
               res.sdar.sic,
               res.L0.CD.cv,
               res.L0.CD.cv.sic,
               res.L0.PSI.cv,
               res.L0.PSI.cv.sic,
               res.mcp,
               res.mcp.sic,
               res.scad,
               res.scad.sic,
               res.L1,
               res.L1.sic
  )
  res
}

## High-dimensional comparison
sim_once_high <- function(seed, SNR = 1, type = 1, n, p, k = 15, sigma = 1, npara = 20, neg = TRUE,
                          method = c("sdar", "sdar-sic", "abess",
                                     "L0cv", "L0-sic", "L0PSIcv", "L0PSI-sic",
                                     "L1", "L1-sic", "L1-ebic", "L1-bic",
                                     "scad", "scad-sic", "scad-ebic", "scad-bic",
                                     "mcp", "mcp-sic", "mcp-ebic", "mcp-bic"),
                          warm.start = TRUE, rho = 0.5, exchange.num = 5) {
  set.seed(seed)
  Nonzero = sample(1:p, k)
  if (type == 1) {
    data <- gendata_high_1(n = n, p = p, Nonzero = Nonzero, SNR = SNR)
  } else if (type == 2) {
    data <- gendata_high_2(n = n, p = p, Nonzero = Nonzero, rho = rho, SNR = SNR)
  } else {
  }
  x = data$x
  y = data$y
  Tbeta = data$Tbeta
  npara <- pmin(npara, ncol(x))

  default_res <- rep(0, 4)

  #bess
  if("sdar" %in% method){
    res.sdar = sim.bess(x, y, Tbeta, npara, warm.start = warm.start, flag = FALSE)
  } else {
    res.sdar = default_res
  }

  if("sdar-sic" %in% method){
    res.sdar.sic = sim.bess(x, y, Tbeta, npara, warm.start = warm.start, flag = FALSE, criterion = "sic")
  } else {
    res.sdar.sic = default_res
  }

  #abess
  if("abess" %in% method){
    res.abess = sim.bess(x, y, Tbeta, npara, warm.start = warm.start, exchange.num = exchange.num, flag= TRUE)
  } else {
    res.abess = default_res
  }

  #L0Learn L0
  if("L0cv" %in% method){
    res.L0.CD.cv = sim.l0learn.cv(x, y, Tbeta, npara, algorithm = "CD")
  } else {
    res.L0.CD.cv = default_res
  }

  if("L0-sic" %in% method){
    res.L0.CD.cv.sic = sim.l0learn(x, y, Tbeta, npara, algorithm = "CD")
  } else {
    res.L0.CD.cv.sic = default_res
  }

  if("L0PSIcv" %in% method){
    res.L0.PSI.cv = sim.l0learn.cv(x, y, Tbeta, npara, algorithm = "CDPSI")
  } else {
    res.L0.PSI.cv = default_res
  }

  if("L0PSI-sic" %in% method){
    res.L0.PSI.cv.sic = sim.l0learn(x, y, Tbeta, npara, algorithm = "CDPSI")
  } else {
    res.L0.PSI.cv.sic = default_res
  }

  #glmnet L1
  if("L1" %in% method){
    res.L1= sim.glmnet(x, y, Tbeta, npara, FALSE)
  } else {
    res.L1 = default_res
  }

  if("L1-sic" %in% method){
    res.L1.sic= sim.glmnet(x, y, Tbeta, npara, FALSE, criterion = "sic")
  } else {
    res.L1.sic = default_res
  }

  if("L1-ebic" %in% method){
    res.L1.ebic= sim.glmnet(x, y, Tbeta, npara, FALSE, criterion = "ebic")
  } else {
    res.L1.ebic = default_res
  }

  if("L1-bic" %in% method){
    res.L1.bic= sim.glmnet(x, y, Tbeta, npara, FALSE, criterion = "bic")
  } else {
    res.L1.bic = default_res
  }

  # SCAD
  if ("scad" %in% method) {
    res.scad= sim.scad(x, y, Tbeta, npara)
  } else {
    res.scad = default_res
  }

  if ("scad-sic" %in% method) {
    res.scad.sic = sim.scad(x, y, Tbeta, npara, criterion = "sic")
  } else {
    res.scad.sic = default_res
  }

  if("scad-ebic" %in% method){
    res.scad.ebic= sim.scad(x, y, Tbeta, npara, criterion = "ebic")
  } else {
    res.scad.ebic = default_res
  }

  if("scad-bic" %in% method){
    res.scad.bic= sim.scad(x, y, Tbeta, npara, criterion = "bic")
  } else {
    res.scad.bic = default_res
  }

  # MCP
  if ("mcp" %in% method) {
    res.mcp <- sim.mcp(x, y, Tbeta, npara)
  } else {
    res.mcp <- default_res
  }

  if ("mcp-sic" %in% method) {
    res.mcp.sic <- sim.mcp(x, y, Tbeta, npara, criterion = "sic")
  } else {
    res.mcp.sic <- default_res
  }

  if ("mcp-ebic" %in% method) {
    res.mcp.ebic <- sim.mcp(x, y, Tbeta, npara, criterion = "ebic")
  } else {
    res.mcp.ebic <- default_res
  }

  if ("mcp-bic" %in% method) {
    res.mcp.bic <- sim.mcp(x, y, Tbeta, npara, criterion = "bic")
  } else {
    res.mcp.bic <- default_res
  }

  res <- rbind(res.abess,
               res.sdar,
               res.sdar.sic,
               res.L0.CD.cv,
               res.L0.CD.cv.sic,
               res.L0.PSI.cv,
               res.L0.PSI.cv.sic,
               res.mcp,
               res.mcp.sic,
               res.mcp.ebic,
               res.mcp.bic,
               res.scad,
               res.scad.sic,
               res.scad.ebic,
               res.scad.bic,
               res.L1,
               res.L1.sic,
               res.L1.ebic,
               res.L1.bic
  )
  res
}
