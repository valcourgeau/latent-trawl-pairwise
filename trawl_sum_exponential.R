SumExponentialTrawl <- new.env()

SumExponentialTrawl$SumExp <- function(param, h, trawl_f){
  n_sup <- length(param)
  exp_trawl_params <- seq(from=1/n_sup, to=1, length.out = n_sup)
  trawl_val <- vapply(exp_trawl_params, function(par){
    trawl_f(param = par, h)
  }, rep(1.0, length(h)))    
  trawl_val <- trawl_val * matrix(rep(param/sum(param), nrow(trawl_val)), ncol=n_sup, byrow = T)
  return(rowSums(trawl_val))
}

SumExponentialTrawl$TrawlB1 <- function(param, h){
  return(SumExponentialTrawl$SumExp(param, h, ExponentialTrawl$TrawlB1))
}

# must have 
# ExponentialTrawl$TrawlB1(1.0, 0:10)
# equal to 
# SumExponentialTrawl$TrawlB1(c(1.0), h=0:10)


SumExponentialTrawl$TrawlB2 <- function(param, h){
  return(SumExponentialTrawl$SumExp(param, h, ExponentialTrawl$TrawlB2))
}

# ExponentialTrawl$TrawlB2(1.0, 0:10)
# SumExponentialTrawl$TrawlB2(c(1.0), h=0:10)

SumExponentialTrawl$TrawlB3 <- function(param, h){
  return(SumExponentialTrawl$SumExp(param, h, ExponentialTrawl$TrawlB3))
}

# ExponentialTrawl$TrawlB3(1.0, 0:10)
# SumExponentialTrawl$TrawlB3(c(1.0), h=0:10)
