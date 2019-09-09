library('zeallot')

CheckAllNonpositive <- function(elems){
  all(elems <= 0.0)
}

CheckAllPositive <- function(elems){
  all(elems > 0.0)
}

StandTrawlTerms <- function(alpha, elems){
  A <- sum(elems[1:2])
  return(-alpha*elems/A)
}

CaseZeroZero <- function(alpha, beta, kappa, B1, B2, B3){
  stopifnot(CheckAllPositive(c(beta, kappa)))
  # c(b1, b2, b3) %<-%(StandTrawlTerms(alpha, c(B1, B2, B3)))
  b_values <- StandTrawlTerms(alpha, c(B1, B2, B3))
  b1 <- b_values[1]
  b2 <- b_values[2]
  b3 <- b_values[3]
  
  stopifnot(CheckAllNonpositive(c(b1, b2, b3)))
  
  tmp <- 1 - 2* (1+kappa/beta)^(-alpha) + (1+kappa/beta)^(b1+b3)*(1+2*kappa/beta)^(b2)
  return(tmp)
}


CaseOneZero <- function(xs, alpha, beta, kappa, B1, B2, B3){
  stopifnot(CheckAllPositive(c(beta, kappa)))
  
  xs <- (xs[1] == 0.0) * xs[2:1] + (xs[1] != 0.0) * xs
  
  x <- xs[1]
  # c(b1, b2, b3) %<-%(StandTrawlTerms(alpha, c(B1, B2, B3)))
  b_values <- StandTrawlTerms(alpha, c(B1, B2, B3))
  b1 <- b_values[1]
  b2 <- b_values[2]
  b3 <- b_values[3]
  
  CheckAllNonpositive(c(b1, b2, b3))
  tmp <- alpha/beta*(1+kappa/beta)^(-alpha-1)
  tmp <- tmp + 1/beta*(1+(kappa+x)/beta)^(b1-1)*(1+(2*kappa+x)/beta)^(b2-1)*(1+kappa/beta)^(b1)*(-alpha*(1+(kappa+x)/beta)+b1*kappa/beta)
  return(tmp)
}

CaseOneOne <- function(xs, alpha, beta, kappa, B1, B2, B3){
  stopifnot(CheckAllPositive(c(beta, kappa)))
  # c(b1, b2, b3) %<-%(StandTrawlTerms(alpha, c(B1, B2, B3)))
  b_values <- StandTrawlTerms(alpha, c(B1, B2, B3))
  b1 <- b_values[1]
  b2 <- b_values[2]
  b3 <- b_values[3]
  
  stopifnot(CheckAllNonpositive(c(b1, b2, b3)))
  
  # c(x1, x2) %<-%(xs)
  x1 <- xs[1]
  x2 <- xs[2]
  
  tmp_1 <- 1/beta^2*(1+(2*kappa+x1+x2)^(b2-2))*(1+(kappa+x1))^(b1-1)*(1+(kappa+x2))^(b3-1)
  tmp_2 <- (b2*(b2-1.0)*(1+(kappa+x1)/beta)*(1+(kappa+x2)/beta)+b1^2*(1+(2*kappa+x1+x2)/beta)^2)
  tmp_3 <- b2*b1*(1+(2*kappa+x1+x2))*((1+(kappa+x1)/beta)+(1+(kappa+x2)/beta))
  
  return(tmp_1 * (tmp_2 + tmp_3))
}

CaseSeparator <-  function(xs, alpha, beta, kappa, B1, B2, B3){
  stopifnot(length(xs) == 2)
  stopifnot(CheckAllNonpositive(-xs))
  # 
  # tmp <- (xs[1] == 0.0)*(xs[2] == 0.0)*CaseZeroZero(alpha, beta, kappa, B1, B2, B3)
  # tmp <- tmp + ((xs[1] != 0.0)*(xs[2] == 0.0) | (xs[1] == 0.0)*(xs[2] != 0.0))*CaseOneZero(xs, alpha, beta, kappa, B1, B2, B3)
  # tmp <- tmp + CheckAllPositive(xs)*CaseOneOne(xs, alpha, beta, kappa, B1, B2, B3)
  #   
  # return(tmp)
  if(CheckAllNonpositive(xs)){
    return(CaseZeroZero(alpha, beta, kappa, B1, B2, B3))
  }else{
    if(prod(xs) == 0.0){
      return(CaseOneZero(xs, alpha, beta, kappa, B1, B2, B3))
    }else{
      return(CaseOneOne(xs, alpha, beta, kappa, B1, B2, B3))
    }
  }
}

TrawlExpB1 <- function(param, h){
  stopifnot(length(param) == 1)
  stopifnot(CheckAllPositive(c(param, h)))
  
  return(exp(-param*h)/param)
}

TrawlExpB2 <- function(param, h){
  stopifnot(length(param) == 1)
  stopifnot(CheckAllPositive(c(param, h)))
  
  return((1.0- exp(-param*h))/param)
}

TrawlExpB3 <- function(param, h){
  return(TrawlExpB1(param, h))
}

GetTrawlFunctions <- function(type){
  return(switch (type,
                 'exp' = c(TrawlExpB1, TrawlExpB2, TrawlExpB3)
  ))
}

PairPDFConstructor <- function(params_noven, type='exp'){
  # params is (alpha, beta, kappa, trawl_params)
  c(B1_func, B2_func, B3_func) %<-% GetTrawlFunctions(type)
  alpha <- params_noven[1]
  beta <- params_noven[2]
  kappa <- params_noven[3]
  trawl_params <- params_noven[4:length(params_noven)]
  
  return(function(xs, h){
    CaseSeparator(xs,
                  alpha = alpha,
                  beta = beta,
                  kappa = kappa,
                  B1 = B1_func(trawl_params, h),
                  B2 = B2_func(trawl_params, h),
                  B3 = B3_func(trawl_params, h)
                  )
  })
}

library(compiler)
library(parallel)

PLConstructor <- function(depth, pair_likehood){
  # returns function implementing Consecutive PL with depth depth
  stopifnot(depth > 1)
  pl_f <- function(data){
    n_sample <- length(data)
    this_pl <- cmpfun(pair_likehood)
    cores <- detectCores(logical = TRUE)
    
    cl <- makeCluster(cores)
    clusterExport(cl, c('CaseSeparator',
                        'data',
                        'CheckAllNonpositive',
                        'CaseOneOne',
                        'CaseOneZero',
                        'CaseZeroZero',
                        'CheckAllPositive',
                        'StandTrawlTerms',
                        'TrawlExpB1',
                        'TrawlExpB2',
                        'TrawlExpB3'))
    
    log_pl_per_depth <- vapply(1:depth, # loop through depths
               FUN = function(k){
                  xs_stack <- cbind(data[1:(n_sample-k)], data[(k+1):(n_sample)])
                  return(
                      sum(
                        unlist(
                          parApply(
                            cl,
                            X = xs_stack, 
                            MARGIN = 1, 
                            FUN = function(xs){
                                  tmp_pl_val <-log(this_pl(xs, h=k))
                                   return(log(this_pl(xs, h = k)))
                                   })
                        )
                    ))
                },
               FUN.VALUE = 1.0)
    return(sum(log_pl_per_depth))
  }
  
  return(pl_f)
}

TrawlPLStandard <- function(params, depth, type='exp'){
  # param with (xi, sigma, kappa, trawl_params)
  c(B1_func, B2_func, B3_func) %<-% GetTrawlFunctions(type)
  params_tmp <- params
  params_tmp[1] <- 1.0/params[1]
  params_tmp[2] <- params[2]/abs(params[1]) - params[3]
  
  pair_likehood_f <- PairPDFConstructor(params_noven = params_tmp, type = type) # yields a function of (xs, h)
  wrapper_with_jacobian <- function(data, h){
    return(pair_likehood_f(data, h)*abs(params[1])^(-3))
  }
  return(PLConstructor(depth = depth, pair_likehood = wrapper_with_jacobian))
}

TrawlPLNoven <- function(params, depth, type='exp'){
  # param with (xi, sigma, kappa, trawl_params)
  c(B1_func, B2_func, B3_func) %<-% GetTrawlFunctions(type)
  params_tmp <- params
  
  pair_likehood_f <- PairPDFConstructor(params_noven = params_tmp, type = type)
  wrapper_with_jacobian <- function(data, h){
    return(pair_likehood_f(data, h)*abs(params[1])^(-3))
  }
  
  return(PLConstructor(depth = depth, pair_likehood = wrapper_with_jacobian))
}

  
TrawlPLFunctional <- function(params, depth, type='exp', parametrisation='standard'){
  PLOperator <- switch(parametrisation,
    'standard' = TrawlPLStandard(params = params, depth = depth, type = type),
    'noven' =  TrawlPLNoven(params = params, depth = depth, type = type)
  )
  
  # multiply by - 1
  return(function(data){
    return(-1*PLOperator(data))})
}

TrawlPL <- function(data, depth, type='exp', parametrisation='standard'){
  return(function(params){
    pl_functional <- TrawlPLFunctional(params = params,
                        depth = depth,
                        type=type,
                        parametrisation=parametrisation)
    return(pl_functional(data))
  })
}

ok <- TrawlPLFunctional(params = c(1,10,1,1), depth = 4, parametrisation='noven')

#set up as function of params
set.seed(42)
rdm_data <- pmax(rnorm(1000), 0.0)
ok(rdm_data)/length(rdm_data)

library(profvis)
profvis({
  ok <- TrawlPLFunctional(params = c(1,10,1,1), depth = 4, parametrisation='noven')

  #set up as function of params
  set.seed(42)
  rdm_data <- pmax(rnorm(1000), 0.0)
  ok(rdm_data)/length(rdm_data)
})
profvis({
  test_params <- c(1,10,1,1)
  #set up as function of params
  set.seed(42)
  rdm_data <- pmax(rnorm(100000), 0.0)
  dac <- TrawlPL(data = rdm_data, depth = 4, parametrisation='noven')

  print(dac(test_params)/length(rdm_data))
})




