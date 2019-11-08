library('zeallot')
library(compiler)
library(parallel)
library(assertthat)
library(magrittr)
# Rcpp::sourceCpp('cpp_core_utils.cpp')
library('ev.trawl.cpp')
source('custom_mle.R')
source('trawl_exponential.R')
source('trawl_gamma.R')
source('trawl_sup_ig.R')
source('trawl_sum_exponential.R')
source('trawl_utils.R')
source('simulation.R')
source('sim_study_utils.R')
source('utils.R')
source('dupuis_regression.R')
source('pairwise_transformation.R')


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


PairPDFConstructor <- function(params_noven, type='exp'){
  # params is (alpha, beta, kappa, trawl_params)
  B_funcs <- GetTrawlFunctions(type)
  B1_func <- B_funcs[[1]]
  B2_func <- B_funcs[[2]]
  B3_func <- B_funcs[[3]]
  alpha <- params_noven[1]
  beta <- params_noven[2]
  kappa <- params_noven[3]
  trawl_params <- params_noven[4:length(params_noven)]
  assertthat::assert_that(CheckAllPositive(trawl_params))
  return(function(xs, h){
    return(ev.trawl.cpp::CppCaseSeparator(xs,
                            alpha = alpha,
                            beta = beta,
                            kappa = kappa,
                            B1 = B1_func(trawl_params, h),
                            B2 = B2_func(trawl_params, h),
                            B3 = B3_func(trawl_params, h)
    ))
  })
}


PLConstructor <- function(depth, pair_likehood, parallel=TRUE, jacob_transform=NULL){
  # returns function implementing Consecutive PL with depth depth
  stopifnot(depth >= 1)
  pl_f <- function(data){
    n_sample <- length(data)
    this_pl <- cmpfun(pair_likehood)
    # this_pl <- pair_likehood
    
    if(parallel){
      cores <- detectCores(logical = TRUE)
      cl <- makeCluster(cores)
      clusterExport(cl, c('CaseSeparator',
                          'CppCaseSeparator',
                          'data',
                          'CheckAllNonpositive',
                          'CppCaseOneOne',
                          'CppCaseOneZero',
                          'CppCaseZeroZero',
                          'CheckAllPositive',
                          'StandTrawlTerms',
                          'TransformationMapInverse',
                          'TransformationMap',
                          'TransformationJacobian',
                          'ParametrisationTranslator',
                          GetTrawlEnvsList()))
      # clusterEvalQ(cl, c(ExponentialTrawl, SumExponential))
  
      log_pl_per_depth <- vapply(1:depth, # loop through depths
                 FUN = function(k){
                    xs_stack <- cbind(data[1:(n_sample-k)], data[(k+1):(n_sample)])
                    return(
                        sum(
                          unlist(
                            parallel::parApply(
                              cl,
                              X = xs_stack, 
                              MARGIN = 1, 
                              FUN = function(xs){
                                pl_val <- this_pl(xs, h=k)
                                if(is.nan(pl_val)){cat('NA', xs, '\n'); return(-10)}
                                if(pl_val < 0.0){
                                  # warning(paste('negative PL', this_pl(xs, h=k), '\n'))
                                  return(pl_val)
                                }else{
                                  return(log(max(pl_val, 1e-9)))
                                }})
                          )
                      ))
                  },
                 FUN.VALUE = 1.0)
    
      parallel::stopCluster(cl)
    }else{
      log_pl_per_depth <- vapply(1:depth, # loop through depths
                                 FUN = function(k){
                                   xs_stack <- cbind(data[1:(n_sample-k)], data[(k+1):(n_sample)])
                                   return(
                                     sum(
                                         apply(
                                           X = xs_stack, 
                                           MARGIN = 1, 
                                           FUN = function(xs){
                                             pl_val <- this_pl(xs, h=k)
                                             if(is.nan(pl_val)){print(xs);return(-0.5)}
                                             # print(pl_val)
                                             if(pl_val < 0.0){
                                               # warning(paste('negative PL', this_pl(xs, h=k), '\n'))
                                               return(pl_val)
                                             }else{
                                               return(log(max(pl_val, 1e-9)))
                                             }})
                                     ))
                                 },
                                 FUN.VALUE = 1.0)
    }
    
    if(!is.null(jacob_transform)){
      log_pl_per_depth <- log_pl_per_depth + vapply(1:depth, # loop through depths
                                                    FUN = function(k){
                                                      xs_stack <- cbind(data[1:(n_sample-k)], data[(k+1):(n_sample)])
                                                      return(
                                                        sum(
                                                          log(jacob_transform(xs_stack))
                                                        ))
                                                    },
                                                    FUN.VALUE = 1.0)
    }
    return(sum(log_pl_per_depth))
  }
  
  return(pl_f)
}

TrawlPLStandard <- function(params, depth, type='exp', parallel=TRUE){
  # param with (xi, sigma, kappa, trawl_params)
  B_funcs <- GetTrawlFunctions(type)
  B1_func <- B_funcs[[1]]
  B2_func <- B_funcs[[2]]
  B3_func <- B_funcs[[3]]
  
  # TODO add translator here
  params_tmp <- params
  params_tmp[1] <- 1.0/params[1]
  params_tmp[2] <- params[2]/abs(params[1]) - params[3]
  cat('Standard params:', params, '\n')
  cat('Noven params:', params_tmp, '\n')
  
  pair_likehood_f <- PairPDFConstructor(params_noven = params_tmp, type = type) # yields a function of (xs, h)
  wrapper_with_jacobian <- function(data, h){
    if(params_tmp[2] <= 0.0){
      return(-1000)
    }else{
      return(pair_likehood_f(data, h)*abs(params[1])^(3))
    }
  }
  
  cat('TrawlPLStandard', parallel, '\n')
  
  return(PLConstructor(depth = depth, pair_likehood = wrapper_with_jacobian, parallel=parallel))
}



TrawlPLStandardTrf <- function(params, depth, type='exp', parallel=TRUE, target_alpha=3.0){
  # param with (xi, sigma, kappa, trawl_params)
  B_funcs <- GetTrawlFunctions(type)
  B1_func <- B_funcs[[1]]
  B2_func <- B_funcs[[2]]
  B3_func <- B_funcs[[3]]
  
  # TODO add translator here
  params_tmp <- params
  params_tmp[1] <- 1.0/params[1]
  params_tmp[2] <- params[2]/abs(params[1]) - params[3]
  
  params_tmp <- ParametrisationTranslator(params = params, parametrisation = 'standard', target = 'noven')
  params_trf <- ParametrisationTranslator(params = params, parametrisation = 'standard', target = 'transform', target_alpha = target_alpha)
  
  cat('Standard params:', params, '\n')
  cat('Noven params:', params_tmp, '\n')
  cat('Trf params:', params_trf, '\n')
  
  
  pair_likehood_f <- PairPDFConstructor(params_noven = params_trf, type = type) # yields a function of (xs, h)
  jacob <- TransformationJacobian(params_std = params, params_trf = params_trf, target_alpha = target_alpha) # yields a function of x
  
  wrapper_with_jacobian_trf <- function(data, h){
    if(params_tmp[2] <= 0.0){
      return(-1000)
    }else{
      return(pair_likehood_f(TransformationMap(data, params_std = params, params_trf = params_trf), h)*abs(params[1])^(3))
    }
    
  }
  return(PLConstructor(depth = depth, pair_likehood = wrapper_with_jacobian_trf, parallel=parallel, jacob_transform = jacob))
}

TrawlPLNoven <- function(params, depth, type='exp', parallel=TRUE){
  # param with (xi, sigma, kappa, trawl_params)
  B_funcs <- GetTrawlFunctions(type)
  B1_func <- B_funcs[[1]]
  B2_func <- B_funcs[[2]]
  B3_func <- B_funcs[[3]]
  params_tmp <- params
  
  pair_likehood_f <- PairPDFConstructor(params_noven = params_tmp, type = type)
  # wrapper_with_jacobian <- function(data, h){
  #   return(pair_likehood_f(data, h)*abs(params[1])^(-3))
  # }
  # 
  return(PLConstructor(depth = depth, pair_likehood = pair_likehood_f, parallel=parallel))
}

  
TrawlPLFunctional <- function(params, depth, type='exp', parametrisation='standard',
                              parallel=TRUE, target_alpha=3.0){
  PLOperator <- switch(parametrisation,
    'standard' = TrawlPLStandard(params = params, depth = depth, type = type, parallel = parallel),
    'std_trf' = TrawlPLStandardTrf(params = params, depth = depth, type = type, parallel = parallel, target_alpha = target_alpha),
    'noven' =  TrawlPLNoven(params = params, depth = depth, type = type, parallel = parallel)
  )
  
  # multiply by - 1
  return(function(data){
    return((-1)*PLOperator(data))})
}

TrawlPL <- function(data, depth, type='exp', parametrisation='standard', parallel=TRUE){
  return(function(params){
    pl_functional <- TrawlPLFunctional(params = params,
                        depth = depth,
                        type=type,
                        parametrisation=parametrisation,
                        parallel=parallel) # returns a function of data
    return(pl_functional(data))
  })
}



