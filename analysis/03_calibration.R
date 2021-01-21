################################################################################ 
# This script calibrates the Sick-Sicker state-transition model (STM) to       #
# epidemiological targets using a Bayesian approach with the Incremental       #
# Mixture Importance Samping (IMIS) algorithm                                 #
#                                                                              # 
# Authors:                                                                     #
#     - Fernando Alarid-Escudero, PhD, <fernando.alarid@cide.edu>              # 
#     - Eline Krijkamp, MS                                                     #
#     - Petros Pechlivanoglou, PhD                                             #
#     - Hawre Jalal, MD, PhD                                                   #
#     - Eva A. Enns, PhD                                                       # 
################################################################################ 
# The structure of this code is according to the DARTH framework               #
# https://github.com/DARTH-git/darthpack                                       #
################################################################################ 

rm(list = ls()) # to clean the workspace

#### 03.1 Load packages, data and functions ####
#### 03.1.1 Load packages and functions ####
# Dependencies have been loaded with 'darthpack'
library(logitnorm)

#### 03.1.2 Load inputs ####
l_params_inmit_calib <- load_params_init(n_age_init = 75, 
                                         n_age_max = 80)
l_params_all_calib <- load_all_params(l_params_init = l_params_inmit_calib)

#### 03.1.3 Load functions ####
# no required functions

#### 03.1.4 Load calibration targets ####
data("03_calibration_targets")

#### 03.2 Visualize targets ####
### TARGET 1: Survival ("Surv")
plotrix::plotCI(x    = SickSicker_targets$Surv$Time, 
                y    = SickSicker_targets$Surv$value, 
                ui   = SickSicker_targets$Surv$ub,
                li   = SickSicker_targets$Surv$lb,
                ylim = c(0, 1), 
                xlab = "Time", ylab = "Pr(Alive)")

### TARGET 2: Prevalence ("Prev")
plotrix::plotCI(x    = SickSicker_targets$Prev$Time, 
                y    = SickSicker_targets$Prev$value, 
                ui   = SickSicker_targets$Prev$ub,
                li   = SickSicker_targets$Prev$lb,
                ylim = c(0, 1), 
                xlab = "Time", ylab = "Pr(Sick+Sicker)")

### TARGET 3: Proportion who are Sicker ("PropSicker"), among all those 
###           afflicted (Sick+Sicker)
plotrix::plotCI(x    = SickSicker_targets$PropSick$Time, 
                y    = SickSicker_targets$PropSick$value, 
                ui   = SickSicker_targets$PropSick$ub,
                li   = SickSicker_targets$PropSick$lb,
                ylim = c(0, 1), 
                xlab = "Time", ylab = "Pr(Sicker | Sick+Sicker)")

#### 03.3 Run calibration algorithms ####
# Check that it works
v_params_calib <- c(r_DieMets = 0.042, # 0.03870286, 
                    r_RecurCDX2pos = (-log(0.80)/60),#0.003328773, 
                    hr_RecurCDX2neg = 2.73, # 3.601069078,
                    p_Mets  = 0.9)#0.980840626)
calibration_out(v_params_calib = v_params_calib, 
                l_params_all = l_params_all_calib)

#### 03.3.1 Specify calibration parameters ####
### Specify seed (for reproducible sequence of random numbers)
set.seed(072218)

### Number of random samples to obtain from the posterior distribution 
n_resamp <- 1000

### Names and number of input parameters to be calibrated
v_param_names  <- names(v_params_calib)
n_param        <- length(v_param_names)

### Vector with range on input search space
# Lower bound
v_lb <- c(r_DieMets       = 0.037, # O'Connell 2004 JNCI Stg IV Fig1 & Fig2;
          r_RecurCDX2pos  = 0.001,
          hr_RecurCDX2neg = 1.58, 
          p_Mets          = 0.9)  
# Upper bound
v_ub <- c(r_DieMets       = -log(1-(1-0.03))/60, # Rutter 2013 JNCI Table 4 5yr RS Colon cancer Stage IV 80+ Lower bound
          r_RecurCDX2pos  = 0.03,
          hr_RecurCDX2neg = 4.72, 
          p_Mets          = 0.99)  

### Number of calibration targets
v_target_names <- c("DFS CDX2neg",
                    "DFS CDX2pos",
                    "OS CDX2neg",
                    "OS CDX2pos",
                    "DSS CDX2neg",
                    "DSS CDX2pos")
n_target       <- length(v_target_names)

#### 03.3.2 Run IMIS algorithm ####
l_fit_imis <- IMIS::IMIS(B        =  1000,      # incremental sample size at each iteration of IMIS
                         B.re     =  n_resamp,  # desired posterior sample size
                         number_k =  10,        # maximum number of iterations in IMIS
                         D        =  0)
### Obtain posterior
m_calib_post <- l_fit_imis$resample

#### 03.4 Exploring posterior distribution ####
#### 03.4.1 Summary statistics of posterior distribution ####
### Compute posterior mean
v_calib_post_mean <- colMeans(m_calib_post)

### Compute posterior median and 95% credible interval
m_calib_post_95cr <- matrixStats::colQuantiles(m_calib_post, 
                                               probs = c(0.025, 0.5, 0.975))

### Compute posterior values for draw
v_calib_post      <- exp(log_post(m_calib_post))

### Compute maximum-a-posteriori (MAP) as the mode of the sampled values
v_calib_post_map  <- m_calib_post[which.max(v_calib_post), ]

# Summary statistics
df_posterior_summ <- data.frame(
  Parameter = v_param_names,
  Mean      = v_calib_post_mean,
  m_calib_post_95cr,
  MAP       = v_calib_post_map,
  check.names = FALSE)
df_posterior_summ

### Save summary statistics of posterior distribution
## As .RData
save(df_posterior_summ, 
     file = "data/03_summary_posterior.RData")
## As .csv
write.csv(df_posterior_summ, 
          file = "tables/03_summary_posterior.csv", 
          row.names = FALSE)

#### 03.4.2 Visualization of posterior distribution ####
### Rescale posterior to plot density of plots
v_calib_alpha <- scales::rescale(v_calib_post)

### Plot the 1000 draws from the posterior
png("figs/03_posterior_distribution_joint.png", 
    width = 8, height = 6, units = 'in', res = 300)
  s3d <- scatterplot3d::scatterplot3d(x = m_calib_post[, 1],
                       y = m_calib_post[, 2],
                       z = m_calib_post[, 3],
                       color = scales::alpha("black", v_calib_alpha),
                       xlim = c(v_lb[1], v_ub[1]), 
                       ylim = c(v_lb[2], v_ub[2]), 
                       zlim = c(v_lb[3], v_ub[3]),
                       xlab = v_param_names[1], 
                       ylab = v_param_names[2], 
                       zlab = v_param_names[3])
  ## Add center of Gaussian components
  s3d$points3d(l_fit_imis$center, col = "red", pch = 8)
  ## Add legend
  legend(s3d$xyz.convert(0.05, 1.0, 5), 
         col = c("black", "red"), 
         bg = "white", pch = c(1, 8), yjust = 0, 
         legend = c("Posterior sample", "Center of Gaussian components"), 
         cex = 1.1)
dev.off()

### Plot the 1000 draws from the posterior with marginal histograms
png("figs/03_posterior_distribution_marginal.png", 
    width = 8, height = 6, units = 'in', res = 300)
  psych::pairs.panels(m_calib_post)
dev.off()

#### 03.5 Store posterior and MAP from IMIS calibration ####
save(m_calib_post,
     v_calib_post_map,
     file = "output/03_imis_output.RData")