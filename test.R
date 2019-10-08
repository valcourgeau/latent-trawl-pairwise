example_params_noven <- c(6.33, 20.12, 12.18, 0.27) 
set.seed(42)
o <- ExceedancesSimulation(params = example_params_noven, parametrisation='noven', n = 5000,
                           vanishing_depth = 200, type = 'exp')
# TODO Check trawl last one generated
(o[,1]>0) %>% mean
acf(o[,1])
acf(as.numeric(o[,1]>0))
lines(0:10, ((o[,1]>0) %>% mean)^{1-exp(-0.27*0:10)}, type='b')
lines(0:10, (0.05)^{1-exp(-0.27*0:10)})
lines(0:50, acf_trawl_num_approx(0:50, 6.33, 20.12, 12.18, 0.27))

acf(as.numeric(o>0))
lines(0:50, acf_trawl_num_approx(0:50, 6.33, 20.12, 12.18, 0.1), col='red')
plot(o)

is_pos <- as.numeric(o>0)
flags_chunk_2 <- zoo::rollapply(is_pos, width=2, FUN=function(chunk){(chunk[1])*prod(chunk)}) > 0
# flags_chunk_2 <- as.numeric(flags_chunk_2)
abline(v=which(flags_chunk_2>0))
flags_chunk_3 <- zoo::rollapply(is_pos, width=3, FUN=function(chunk){(chunk[1])*prod(chunk)}) > 0
# flags_chunk_3 <- as.numeric(flags_chunk_3)
abline(v=which(flags_chunk_3>0), col='red')
acf(o)
acf(is_pos)

ev_fit <- EVTrawlFit(o, depth = 5, parametrisation = 'std_transform', type = 'exp')
ev_fit


# Testing the correlation between latent and observed
set.seed(42)
sim_vals_and_latent <- ExceedancesSimulation(params = example_params_noven,
                                             parametrisation='noven', n = 20000,
                                             vanishing_depth = 100, type = 'exp')
latent <- sim_vals_and_latent[,2]
obs <- sim_vals_and_latent[,1]
mean(obs>0)
which_sim_pos <- which(obs>0)

# is_pos <- as.numeric(which_sim_pos)
is_pos <- as.numeric(obs>0)
flags_chunk_2 <- zoo::rollapply(is_pos, width=2, FUN=function(chunk){sum(chunk)*(chunk[1])*prod(chunk)}) > 0
flags_chunk_3 <- zoo::rollapply(is_pos, width=3, FUN=function(chunk){sum(chunk)*(chunk[1])*prod(chunk)}) > 0


par(mfrow=c(3,4))
for(shift in c(0:11)){
  plot(latent[which_sim_pos[1:(length(which_sim_pos)-shift)]],
       obs[which_sim_pos[(1+shift):length(which_sim_pos)]] %>% log,
       xlab='Latent', ylab='log-Observed', main=paste('Shift of', shift, 'corr', round(cor(
         latent[which_sim_pos[1:(length(which_sim_pos)-shift)]],
         obs[which_sim_pos[(1+shift):length(which_sim_pos)]] %>% log
       ), 2)))
}

par(mfrow=c(1,3))
# with zeros
plot(diff(latent)[which_sim_pos], diff(obs)[which_sim_pos], main = paste('First elem of increment is positive: corr', round(cor(
  diff(latent)[which_sim_pos], diff(obs)[which_sim_pos], use = 'pairwise'
), digits = 2)), lwd=1.4, cex.axis=1.3, cex.lab=1.3,cex.main=1.5, cex=1.5,
xlab='Diff in latent', ylab='Diff in obs')
cor(diff(latent)[which_sim_pos], diff(obs)[which_sim_pos], use = 'pairwise')
lm(diff(obs)[which_sim_pos] ~ diff(latent)[which_sim_pos], na.action = na.omit) %>% summary


plot(diff(latent)[flags_chunk_2], diff(obs)[flags_chunk_2], main = paste('2 Consecutive extremes: corr', round(cor(
  diff(latent)[flags_chunk_2], diff(obs)[flags_chunk_2]
), digits = 2)), lwd=1.4, cex.axis=1.3, cex.lab=1.3, cex.main=1.5, cex=1.5,
xlab='Diff in latent', ylab='Diff in obs')
cor(diff(latent)[flags_chunk_2], diff(obs)[flags_chunk_2])
lm(diff(obs)[flags_chunk_2] ~ diff(latent)[flags_chunk_2]) %>% summary

plot(diff(latent)[flags_chunk_3], diff(obs)[flags_chunk_3], main = paste('3 Consecutive extremes: corr', round(cor(
  diff(latent)[flags_chunk_3], diff(obs)[flags_chunk_3]
), digits = 2)), lwd=1.4, cex.axis=1.3, cex.lab=1.3,cex.main=1.5, cex=1.5,
xlab='Diff in latent', ylab='Diff in obs')
cor(diff(latent)[flags_chunk_3], diff(obs)[flags_chunk_3])
lm(diff(obs)[flags_chunk_3] ~ diff(latent)[flags_chunk_3]) %>% summary


# test with complement probability to simulate
example_params_noven <- c(6.33, 20.12, 0.15, 0.1) 

set.seed(42)
o <- ExceedancesSimulation(params = example_params_noven, parametrisation='noven', n = 5000, vanishing_depth = 40, type = 'exp')
mean(as.numeric(o[,1]>0))
plot(o[,1])
acf(o[,1])
lines(0:20, acf_trawl_num_approx_inv(0:20, alpha = 6.33, beta = 20.12, kappa = 0.15, 0.1, cov = F))

which_pos <- which(o[,1]>0)
plot(o[which_pos,1], o[which_pos,2])
plot(density(o[which_pos,1]))
lines(0:500/10, eva::dgpd(0:500/10, shape = 1/6.33, scale = (20.12+12.18)/6.33))
