SupIGTrawl <- new.env()

# param is (gamma, delta)

SupIGTrawl$LebA <- function(param){
  return(param[1]/param[2])
}

SupIGTrawl$TrawlB1 <- function(param, h){
  return(SupIGTrawl$LebA(param)*(1-exp(param[1]*param[2]*(1-sqrt(1+2*h/param[1]^2)))))
}


SupIGTrawl$TrawlB2 <- function(param, h){
  return(exp(param[1]*param[2]*(1-sqrt(1+2*h/param[1]^2)))*SupIGTrawl$LebA(param))
}


SupIGTrawl$TrawlB3 <- function(param, h){
  return(SupIGTrawl$TrawlB1(param, h))
}

SupIGTrawl$Config <- function(){
  return(list(n_params=2, lower=rep(1e-5, 2), upper=rep(10.0,2)))
}

