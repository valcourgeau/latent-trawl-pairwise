
PropagationProbability <- function(x_now, x_future, h, xi, sigma, kappa, rho, type="exp", eps=1e-7){
  B_funcs <- GetTrawlFunctions(type)
  B1_func <- B_funcs[[1]]
  B2_func <- B_funcs[[2]]
  beta <- sigma / abs(xi)
  alpha <- 1/xi
  x_now <- (kappa+x_now)/beta
  x_future <- (kappa+x_future)/beta
  return(
      (1+x_now)^{alpha-alpha*B1_func(rho, h)/(B1_func(rho, h)+B2_func(rho, h))} * 
      (1+x_now+x_future)^{-alpha*B2_func(rho, h)/(B1_func(rho, h)+B2_func(rho, h))} * 
      (1+x_future)^{-alpha*B1_func(rho, h)/(B1_func(rho, h)+B2_func(rho, h))}
  )
}

PropagationProbability(x_now = 0,
                       x_future = 0,
                       h = c(1:20),
                       xi = 0.083,
                       sigma = 0.589,
                       kappa = 1.56,
                       rho = 0.17)


EmpiricalPropagation <- function(data, threshold, h=50){
  decider <- data > threshold
  n <- length(decider)
  counters <- vapply(0:h, function(i){
    tmp <- rbind(decider[1:(n-i)], decider[(i+1):n])
    return(
      sum(
        as.numeric(
          apply(tmp, 2, all)
        )
      )
    )
  }, 2.0)
 
  return(counters/counters[1])
}

EmpiricalPropagation(pollution_data[,1], 0, 10)
