GammaTrawl <- new.env()

# param is (alpha, H)

GammaTrawl$LebA <- function(param){
  return(param[1]/(param[2]-1))
}

GammaTrawl$TrawlB1 <- function(param, h){
  return(SupIGTrawl$LebA(param)-GammaTrawl$TrawlB2(param, h))
}


GammaTrawl$TrawlB2 <- function(param, h){
  return(exp((1-param[2])*log(1+h/param[1]))*GammaTrawl$LebA(param))
}


GammaTrawl$TrawlB3 <- function(param, h){
  return(SupIGTrawl$TrawlB1(param, h))
}

GammaTrawl$Config <- function(){
  return(list(n_params=2, lower=c(0.5, 1.7), upper=c(3,4)))
}

