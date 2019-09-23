ExponentialTrawl <- new.env()

ExponentialTrawl$TrawlB1 <- function(param, h){
  stopifnot(length(param) == 1)
  # cat('param B1', param, 'h', h, '\n')
  # assertthat::assert_that(CheckAllPositive(c(param, h)), msg = cat('params', param, '/ h', h,'\n'))
  return((1.0- exp(-param*h))/param)
}

ExponentialTrawl$TrawlB2 <- function(param, h){
  stopifnot(length(param) == 1)
  # cat('param B2', param, 'h', h, '\n')
  # assertthat::assert_that(CheckAllPositive(c(param, h)), msg = cat('params', param, '/ h', h,'\n'))
  return(exp(-param*h)/param)
}

ExponentialTrawl$TrawlB3 <- function(param, h){
  return(ExponentialTrawl$TrawlB1(param, h))
}