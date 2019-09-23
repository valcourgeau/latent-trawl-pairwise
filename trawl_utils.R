GetTrawlFunctions <- function(type){
  select_env <- switch (type,
                 'exp' = ExponentialTrawl,
                 'sum_exp' = SumExponentialTrawl
  )
  
  return(c(select_env$TrawlB1, select_env$TrawlB2, select_env$TrawlB3))
}

GetTrawlFunctions('exp')[1][[1]]
GetTrawlFunctions('sum_exp')[1][[1]]



GetTrawlParamsConfig <- function(type){
  # returns triplet list (n_params, lower, upper) for each trawl type.
  return(switch (type,
                 'exp' = list(n_params=1, lower=0.05, upper=1.0),
                 'sum_exp' = list(n_params=5, lower=1e-5, upper=1.0)
  ))
}

GetTrawlEnvs <- function(){
  return(c(
    'ExponentialTrawl',
    'SumExponentialTrawl'
  ))
}
