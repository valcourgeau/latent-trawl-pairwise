example_params_noven <- c(6.33, 20.12, 12.18, 0.1) #
o <- ExceedancesSimulation(params = example_params_noven, parametrisation='noven', n = 5000, vanishing_depth = 100, type = 'exp')
# TODO Check trawl last one generated
(o>0) %>% mean
acf(o)
lines(0:50, acf_trawl_num_approx(0:50, 6.33, 20.12, 12.18, 0.1)*sqrt(0.05))

acf(as.numeric(o>0))
lines(0:50, acf_trawl_num_approx(0:50, 6.33, 20.12, 12.18, 0.1))
plot(o)
is_pos <- as.numeric(o>0)
flags_chunk_2 <- zoo::rollapply(is_pos, width=2, FUN=function(chunk){sum(chunk)*(chunk[1])*prod(chunk)}) > 0
abline(v=which(flags_chunk_2>0))
flags_chunk_3 <- zoo::rollapply(is_pos, width=3, FUN=function(chunk){sum(chunk)*(chunk[1])*prod(chunk)}) > 0
abline(v=which(flags_chunk_3>0), col='red')
acf(o)
acf(is_pos)

ev_fit <- EVTrawlFit(o, depth = 5, parametrisation = 'standard', type = 'exp')
ev_fit
