GetTrawlFunctions <- function(type){
  return(switch (type,
                 'exp' = c(ExponentialTrawl$TrawlB1, ExponentialTrawl$TrawlB2, ExponentialTrawl$TrawlB3)
  ))
}

GetTrawlParamsConfig <- function(type){
  # returns triplet list (n_params, lower, upper) for each trawl type.
  return(switch (type,
                 'exp' = list(n_params=1, lower=1e-5, upper=1.0)
  ))
}