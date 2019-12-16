
PropagationProbability <- function(x_now, x_future, h, xi, sigma, kappa, rho, type="exp", eps=1e-7){
  B_funcs <- GetTrawlFunctions(type)
  B1_func <- B_funcs[[1]]
  B2_func <- B_funcs[[2]]
  beta <- sigma * abs(xi)
  x_now <- (kappa+x_now)/beta
  x_future <- (kappa+x_future)/beta
  print(B1_func)
  print(B2_func)
  return(
    (1+x_future/(x_now+eps))^{-B2_func(param=rho, h)} * (1.0+x_future)^{-B1_func(param=rho, h)}
  )
}

PropagationProbability(x_now = 0.1,
                       x_future = 0,
                       h = 1:10,
                       xi = 0.3,
                       sigma = 3,
                       kappa = 3,
                       rho = 0.3)
