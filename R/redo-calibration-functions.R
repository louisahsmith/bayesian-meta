myLogLikelihoodErrorModel <- function(theta, logRr, seLogRr, trueLogRr) {
  
  estimateLl <- function(i) {
    mean <- theta[1] + (1 + theta[2]) * trueLogRr[i]
    sd <- theta[3] + theta[4] * abs(trueLogRr[i])
    if (sd < 0) {
      return(Inf)
    }
    else if (sd < 1e-06) {
      return(-dnorm(logRr[i], mean, seLogRr[i], log = TRUE))
    }
    else {
      return(-log(EmpiricalCalibration:::gaussianProduct(logRr[i], mean, seLogRr[i], 
                                                         sd)))
    }
  }
  result <- sum(sapply(1:length(logRr), estimateLl))
  if (is.infinite(result) || is.na(result))
    result <- 9999999999999999999
  result
}

# if trueLogRr > 1
gr_ind_pos <- deriv(expression(-log( (2 * pi)^(-1/2) * 
                                       (seLogRr^2 + (theta3 + theta4 * trueLogRr)^2)^(-1/2) * 
                                       exp(-(logRr - (theta1 + (1 + theta2) * trueLogRr))^2 /
                                             (2 * (seLogRr^2 + (theta3 + theta4 * trueLogRr)^2))))), 
                    c("theta1", "theta2", "theta3", "theta4"),
                    function(theta1, theta2, theta3, theta4, logRr, seLogRr, trueLogRr){})

gr_ind_neg <- deriv(expression(-log((2 * pi)^(-1/2) * 
                                      (seLogRr^2 + (theta3 + theta4 * (-trueLogRr))^2)^(-1/2) * 
                                      exp(-(logRr - (theta1 + (1 + theta2) * trueLogRr))^2 /
                                            (2 * (seLogRr^2 + (theta3 + theta4 * (-trueLogRr))^2))))), 
                    c("theta1", "theta2", "theta3", "theta4"),
                    function(theta1, theta2, theta3, theta4, logRr, seLogRr, trueLogRr){})

gr_lk <- function(theta, logRr, seLogRr, trueLogRr) {
  gr_i <- function(i) {
    if (trueLogRr[i] >= 0) {
      gr_ind_pos(theta[1], theta[2], theta[3], theta[4], logRr[i], seLogRr[i], trueLogRr[i])
    } else {
      gr_ind_neg(theta[1], theta[2], theta[3], theta[4], logRr[i], seLogRr[i], trueLogRr[i])
    }
  }
  result <- rowSums(sapply(1:length(logRr), function(i) attr(gr_i(i),"gradient")))
  result
}

mySystematicErrorModel <- function(logRr, seLogRr, trueLogRr, estimateCovarianceMatrix = FALSE, 
                                   legacy = FALSE, gradient = TRUE, trace = 0, maxit = 1000, method = "BFGS") {
  if (any(is.infinite(seLogRr))) {
    warning("Estimate(s) with infinite standard error detected. Removing before fitting error model")
    trueLogRr <- trueLogRr[!is.infinite(seLogRr)]
    logRr <- logRr[!is.infinite(seLogRr)]
    seLogRr <- seLogRr[!is.infinite(seLogRr)]
  }
  if (any(is.infinite(logRr))) {
    warning("Estimate(s) with infinite logRr detected. Removing before fitting error model")
    trueLogRr <- trueLogRr[!is.infinite(logRr)]
    seLogRr <- seLogRr[!is.infinite(logRr)]
    logRr <- logRr[!is.infinite(logRr)]
  }
  if (any(is.na(seLogRr))) {
    warning("Estimate(s) with NA standard error detected. Removing before fitting error model")
    trueLogRr <- trueLogRr[!is.na(seLogRr)]
    logRr <- logRr[!is.na(seLogRr)]
    seLogRr <- seLogRr[!is.na(seLogRr)]
  }
  if (any(is.na(logRr))) {
    warning("Estimate(s) with NA logRr detected. Removing before fitting error model")
    trueLogRr <- trueLogRr[!is.na(logRr)]
    seLogRr <- seLogRr[!is.na(logRr)]
    logRr <- logRr[!is.na(logRr)]
  }
  theta <- c(0, 0, 0.1, 0)
  parscale <- c(1, 1, 1, 1)
  if (gradient) {
    fit <- optim(theta, myLogLikelihoodErrorModel, gr = gr_lk, logRr = logRr, seLogRr = seLogRr, 
                 trueLogRr = trueLogRr, method = method, hessian = TRUE, 
                 control = list(parscale = parscale, trace = trace, maxit = maxit))
    
  } else {
    fit <- optim(theta, myLogLikelihoodErrorModel, logRr = logRr, seLogRr = seLogRr, 
                 trueLogRr = trueLogRr, method = method, hessian = TRUE, 
                 control = list(parscale = parscale, trace = trace, maxit = maxit))
  }
  if (fit$convergence != 0) warning(paste0("Model did not converge with error code ", fit$convergence))
  model <- fit$par
  names(model) <- c("meanIntercept", "meanSlope", "sdIntercept", "sdSlope")
  class(model) <- "systematicErrorModel"
  model
}

mycalibrateConfidenceInterval <- function(logRr, seLogRr, model, ciWidth = 0.95) {
  result <- data.frame(logRr = rep(0, length(logRr)), logLb95Rr = 0, 
                       logUb95Rr = 0, seLogRr = 0)
  for (i in 1:nrow(result)) {
    if (is.infinite(logRr[i]) || is.na(logRr[i]) || is.infinite(seLogRr[i]) || 
        is.na(seLogRr[i])) {
      result$logRr[i] <- NA
      result$logLb95Rr[i] <- NA
      result$logUb95Rr[i] <- NA
      result$seLogRr[i] <- NA
    }
    else {
      # truth should be obs - beta = obs - (a + b*true) = true
      # true + a + b*true = obs
      # a + (1 + b)*true = obs
      # true = (obs - a) / (1 + b)
      # 0 = (obs - a) / (1 + b) - true
      result$logRr[i] <- (logRr[i] - model[1]) / (1 + model[2])
      # bounds should be true +/ 1.96 * sqrt(seLogRr[i]^2 + (c + d*abs(true))^2)
      result$seLogRr[i] <- sqrt(seLogRr[i]^2 + (model[3] + model[4]*abs(result$logRr[i]))^2)
      result$logLb95Rr[i] <- result$logRr[i] + qnorm((1 - ciWidth)/2)*result$seLogRr[i]
      result$logUb95Rr[i] <- result$logRr[i] + qnorm(1 - (1 - ciWidth)/2)*result$seLogRr[i]
    }
  }
  return(result)
}

mycalibrateCiLoo <- function(subset, leaveOutId) {
  subsetMinusOne <- subset[subset$oldOutcomeId != leaveOutId, ]
  one <- subset[subset$oldOutcomeId == leaveOutId, ]
  model <- mySystematicErrorModel(logRr = subsetMinusOne$logRr,
                                  seLogRr = subsetMinusOne$seLogRr,
                                  trueLogRr = log(subsetMinusOne$trueEffectSize),
                                  estimateCovarianceMatrix = FALSE)
  calibratedCi <- mycalibrateConfidenceInterval(logRr = one$logRr,
                                                seLogRr = one$seLogRr,
                                                model = model)
  one$calibratedRr <- exp(calibratedCi$logRr)
  one$calibratedCi95Lb <- exp(calibratedCi$logLb95Rr)
  one$calibratedCi95Ub <- exp(calibratedCi$logUb95Rr)
  one$calibratedLogRr <- calibratedCi$logRr
  one$calibratedSeLogRr <- calibratedCi$seLogRr
  return(one)
}

mycalibrate <- function(subset) {
  
  pcs <- subset %>%
    filter(.data$targetEffectSize > 1 & !is.na(.data$seLogRr))
  
  if (nrow(pcs) > 5) {
    subset <- purrr::map_dfr(unique(subset$oldOutcomeId), mycalibrateCiLoo, subset = subset)
  } else {
    subset$calibratedRr <- rep(NA, nrow(subset))
    subset$calibratedCi95Lb <- rep(NA, nrow(subset))
    subset$calibratedCi95Ub <- rep(NA, nrow(subset))
    subset$calibratedLogRr <- rep(NA, nrow(subset))
    subset$calibratedSeLogRr <- rep(NA, nrow(subset))
  }
  return(subset)
}
