data_path <- "C:/Users/Valentin/Documents/GitHub/multi-trawl-extremes/data/hourly_bloomsbury_2000_2017.csv"
pdbl <- read.csv(data_path)

# Reducing sample size
pdbl <- pdbl

library(evir)
library(forecast)
library(lubridate)
PLOT.IT <- FALSE
set.seed(42)

raw_data <- pdbl[,-c(1:3)]

if(PLOT.IT){
  stl_tmp <- stl(ts(pdbl[,4], frequency = 24), s.window = 24*7, t.window = 24*30*3) # monthly seasonality, quarterly trend
  plot(stl_tmp)
}

# monthly seasonality, quarterly trend
stl_clean_fct <- function(uni_data){stl(ts(uni_data, frequency = 24), s.window = 24*30, t.window = 24*30*3)$time.series[,3]}
clean_data <- apply(raw_data, MARGIN = 2, FUN = stl_clean_fct)
jittered_clean_data <-  apply(clean_data, MARGIN = 2, FUN =  function(x){
  return(x + rnorm(n = length(x), mean = 0, sd = 0.1 * sd(x)))
})
std_clean_data <- apply(jittered_clean_data, MARGIN = 2, FUN = function(x){x/sd(x)})


# taking exceedances with 95% quantile
final_clean_data <- apply(std_clean_data, MARGIN = 2, FUN = function(x){return(pmax(x-quantile(x, 0.95), 0.0))})
gpd_fit_p_values <- apply(final_clean_data[1:10000,], MARGIN = 2, FUN = function(x){gPdtest::gpd.test(x[x>0], J = 2000)$p.values[1:2]})
# add exponential case
gpd_fit_p_values_95 <- rbind(gpd_fit_p_values, apply(final_clean_data[1:10000,], MARGIN = 2, FUN = function(x){Renext::gofExp.test(x[x > 0])$p.value}))
row.names(gpd_fit_p_values_95) <- c('xi < 0', 'xi > 0', 'xi = 0')
print(gpd_fit_p_values_95)
# CCL: all good except N0

# taking exceedances with 98% quantile
final_clean_data <- apply(std_clean_data, MARGIN = 2, FUN = function(x){return(pmax(x-quantile(x, 0.98), 0.0))})
gpd_fit_p_values <- apply(final_clean_data[1:10000,], MARGIN = 2, FUN = function(x){gPdtest::gpd.test(x[x>0], J = 2000)$p.values[1:2]})
# add exponential case
gpd_fit_p_values_98 <- rbind(gpd_fit_p_values, apply(final_clean_data[1:10000,], MARGIN = 2, FUN = function(x){Renext::gofExp.test(x[x > 0])$p.value}))
row.names(gpd_fit_p_values_98) <- c('xi < 0', 'xi > 0', 'xi = 0')
print(gpd_fit_p_values_98)
# CCL: NO is xi > 0

# FINAL CCL:    
# O3   (exp at 95)   
# CO   (xi > 0 at 95)    
# NO   (exp at 98)     
# NO2  (xi > 0 at 95)    
# PM10 (xi > 0 at 95)    
# SO2  (x > 0  at 95)

# Mean-Excess plot
if(PLOT.IT){
  par(mfrow=c(2,3))
  apply(std_clean_data, MARGIN = 2, FUN = function(x){
    evir::meplot(x)
    abline(v=quantile(x, c(0.95, 0.98)), col=c('lightblue', 'blue'), lwd = 2)
  })
}

clean_data_to_save <- apply(std_clean_data, MARGIN = 2, FUN = function(x){return(pmax(x-quantile(x, 0.95), 0.0))})
clean_data_to_save[,3] <- pmax(std_clean_data[,3]-quantile(std_clean_data[,3], 0.98), 0.0)
write.csv(clean_data_to_save, file = 'data/clean_pollution_data.csv')
write.csv(gpd_fit_p_values_95, file = 'data/gpd_fit_p_values_95.csv')
write.csv(gpd_fit_p_values_98, file = 'data/gpd_fit_p_values_98.csv')
