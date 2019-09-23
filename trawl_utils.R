GetTrawlEnv <- function(type){
  return(switch(type,
                'exp' = ExponentialTrawl,
                'sum_exp' = SumExponentialTrawl,
                'sup_ig' = SupIGTrawl,
                'gamma' = GammaTrawl
  ))
}

GetTrawlFunctions <- function(type){
  select_env <- GetTrawlEnv(type)
  
  return(c(select_env$TrawlB1, select_env$TrawlB2, select_env$TrawlB3))
}

GetTrawlFunctions('exp')[1][[1]]
GetTrawlFunctions('sum_exp')[1][[1]]


GetTrawlParamsConfig <- function(type){
  # returns triplet list (n_params, lower, upper) for each trawl type.
  return(GetTrawlEnv(type)$Config())
}


GetTrawlEnvsList <- function(){
  return(c(
    'ExponentialTrawl',
    'SumExponentialTrawl',
    'SupIGTrawl',
    'GammaTrawl'
  ))
}
