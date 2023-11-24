
library(ssdtools)
library(ssddata)
# library(ggplot2)
# library(MATA)
# library(tidyverse)
use_dists <-  c("gamma", "lgumbel", "llogis", "lnorm", "lnorm_lnorm", "weibull")
set.seed(10)
n_vals <- c(8, 16, 24, 32, 40)
n_sims <- 500
n_boot <- 5000
percent_vals <- c(1, 5, 10, 20)

data("ccme_data")
ccme_datasets <- unique(ccme_data$Chemical)

burr3_sim_datasets <- ssd_fits |> dplyr::filter(Distribution == "BurrIII") |> 
  dplyr::filter(is.na(Filter)) |> 
  dplyr::select(Dataset) |> 
  unique() |> 
  unlist()
# Original model fits ----------------------------------------------------------
burr3_orig_fits <- lapply(burr3_sim_datasets, FUN = function(x){
  data <- get_ssddata(x)
  dist <- ssd_fit_dists(data, 
                        left = 'Conc', dists = 'burrIII3', 
                        silent = TRUE, reweight = FALSE, min_pmix = 0, nrow = 6L, 
                        computable = FALSE, at_boundary_ok = TRUE, rescale = FALSE)
  return(dist)})
names(burr3_orig_fits) <- burr3_sim_datasets  

# Generate random data --------------------------------------------------------- 

# all_burr_sim_dat <- lapply(1:n_sims, FUN = function(s){
#   burr_sim_dat <- lapply(n_vals, FUN = function(n) {
#     burr_sim_dat_n <- lapply(burr3_orig_fits, FUN = function(x){
#       dat_out <- data.frame(Conc = ssdtools::ssd_rburrIII3(n, 
#                                     scale = exp(x$burrIII3$pars$log_scale), 
#                                     shape1 = exp(x$burrIII3$pars$log_shape1),
#                                     shape2 = exp(x$burrIII3$pars$log_shape2)))      
#       return(dat_out)  
#     })
#     names(burr_sim_dat_n) <- burr3_sim_datasets  
#     
#     return(burr_sim_dat_n)
#   })
#   names(burr_sim_dat) <- paste("N", n_vals, sep = "_") 
#   return(burr_sim_dat)
# })
# save(all_burr_sim_dat, file = "sim_data/all_burr_sim_dat.Rdata")
load("sim_data/all_burr_sim_dat.Rdata")

# Fit all simdat using the standard distribution set ---------------------------
# all_sim_fits <- lapply(all_burr_sim_dat, FUN = function(x){
#   fits.x <- lapply(x, FUN = function(y){
#     fits.y <- lapply(y, ssd_fit_dists, 
#                      left = 'Conc', dists = c('gamma', 'lgumbel', 'llogis', 'lnorm',  'weibull'), 
#                      silent = TRUE, reweight = FALSE, min_pmix = 0, nrow = 6L, 
#                      computable = TRUE, at_boundary_ok = FALSE, rescale = FALSE)  
#   })  
# })
# save(all_sim_fits, file = "sim_data/all_burr_sim_fits.Rdata")
load("sim_data/all_burr_sim_fits.Rdata")



# arithmetic weighted mean
ar_mean <- ssd_hc(dist, ci = TRUE, nboot = n_boot, average = TRUE, 
       averaging_method = "arithmetic", percent = percent_vals)
# geometric weighted mean
gm_mean <- ssd_hc(dist, ci = TRUE, nboot = n_boot, average = TRUE, 
       averaging_method = "geometric", percent = percent_vals)

# uniroot bootstrap
uniroot <- ssd_hc(dist, ci = TRUE, nboot = n_boot, average = TRUE, 
       averaging_method = "uniroot", percent = percent_vals)

# weighed_sample
w_sample <- ssd_hc(dist, ci = TRUE, nboot = n_boot, average = TRUE, 
        averaging_method = "weighted_sample", percent = percent_vals)

# multiroot bootstrap
multiroot <- ssd_hc(dist, ci = TRUE, nboot = n_boot, average = TRUE, 
                   averaging_method = "multiroot", percent = percent_vals)


