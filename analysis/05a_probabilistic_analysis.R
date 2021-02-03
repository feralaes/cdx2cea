################################################################################ 
# This script conducts the probabilistic sensitivity analysis (robPSA) of the  #
# cost-effectiveness of testing Stage II colon cancer patients for the absence #
# of CDX2 biomarker followed by adjuvant chemotherapy                          #
#                                                                              # 
# Authors:                                                                     #
#     - Fernando Alarid-Escudero, PhD, <fernando.alarid@cide.edu>              #
#     - Deb Schrag, MD, MPH                                                    #
#     - Karen M. Kuntz, ScD                                                    #
################################################################################
# The structure of this code is according to the DARTH framework               #
# https://github.com/DARTH-git/Decision-Modeling-Framework                     #
################################################################################

rm(list = ls()) # to clean the workspace

#### 05a.1 Load packages and functions ####
#### 05a.1.1 Load packages ####
# PSA functionality
library(dampack)    # decision-analytic modeling visualization tool
library(doParallel)

#### 05a.1.2 Load inputs ####
l_params_all <- load_all_params() # function in cdx2cea

#### 05a.1.3 Load functions ####
# no required functions

#### 05a.2 Cost-effectiveness analysis parameters ####
### Strategy names
v_names_str <- l_params_all$v_names_str
### Number of strategies
n_str <- length(v_names_str)

#### 05a.3 Setup probabilistic analysis ####
### Number of simulations
n_sim <- 1000

### Generate PSA input dataset
df_psa_input <- generate_psa_params(l_params_all = l_params_all)

### Initialize matrices for PSA output 
## Matrix of costs
df_c <- as.data.frame(matrix(0, 
                             nrow = n_sim,
                             ncol = n_str))
colnames(df_c) <- v_names_str
## Matrix of effectiveness
df_e <- as.data.frame(matrix(0, 
                             nrow = n_sim,
                             ncol = n_str))
colnames(df_e) <- v_names_str

#### 05a.4 Conduct probabilistic sensitivity analysis ####
### Run decision model on each parameter set of PSA input dataset to produce
### PSA outputs for cost and effects
n_time_init_psa_series <- Sys.time()
for(i in 1:n_sim){ # i <- 1
  l_psa_input <- update_param_list(l_params_all, df_psa_input[i,])
  df_out_temp <- calculate_ce_out(l_psa_input)
  df_c[i, ] <- df_out_temp$Cost
  df_e[i, ] <- df_out_temp$Effect
  # Display simulation progress
  if(i/(n_sim/100) == round(i/(n_sim/100),0)) {
    cat('\r', paste(i/n_sim * 100, "% done", sep = " "))
  }
}
n_time_end_psa_series <- Sys.time()
n_time_total_psa_series <- n_time_end_psa_series - n_time_init_psa_series
print(paste0("PSA with ", scales::comma(n_sim), " simulations run in series in ", 
             round(n_time_total_psa_series, 2), " ", 
             units(n_time_total_psa_series)))

l_out_psa <- run_probsa(df_psa_input = df_psa_input, 
                        n_str = n_str, 
                        parallel = TRUE)
df_c <- l_out_psa$Costs
df_e <- l_out_psa$Effects

### Create PSA object for dampack
l_psa <- make_psa_obj(cost = df_c, 
                      effectiveness = df_e, 
                      parameters = df_psa_input, 
                      strategies = v_names_str)
# Format object for plotting
l_psa$strategies              <- c("No CDX2 testing and no FOLFOX", "CDX2 testing and FOLFOX if CDX2-negative")
colnames(l_psa$effectiveness) <- c("No CDX2 testing and no FOLFOX", "CDX2 testing and FOLFOX if CDX2-negative")
colnames(l_psa$cost)          <- c("No CDX2 testing and no FOLFOX", "CDX2 testing and FOLFOX if CDX2-negative")

#### 05a.5 Save PSA objects ####
save(df_psa_input, df_c, df_e, v_names_str, n_str,
     l_psa,
     file = "output/05a_psa_dataset.RData")

#### 05a.6 Create probabilistic analysis graphs ####
data("l_psa") # stored as data object in 'cdx2cea'

### Vector with willingness-to-pay (WTP) thresholds
v_wtp <- seq(0, 150000, by = 5000)

#### 05a.6.1 Cost-effectiveness scatter plot ####
plot(l_psa) + 
  theme(legend.position = "bottom")
ggsave("figs/05a_cea_plane_scatter.png", width = 8, height = 6)

#### 05a.6.2 Conduct CEA with probabilistic output ####
### Compute expected costs and effects for each strategy from the PSA
df_out_ce_psa <- summary(l_psa)
### Calculate incremental cost-effectiveness ratios (ICERs)
df_cea_psa <- calculate_icers(cost = df_out_ce_psa$meanCost, 
                              effect = df_out_ce_psa$meanEffect,
                              strategies = df_out_ce_psa$Strategy)
df_cea_psa
### Save CEA table with ICERs
## As .RData
# save(df_cea_psa, 
#      file = "tables/05a_probabilistic_cea_results.RData")
save(df_cea_psa, 
     file = "data/05a_probabilistic_cea_results.RData")
## As .csv
write.csv(df_cea_psa, 
          file = "tables/05a_probabilistic_cea_results.csv")

#### 05a.6.3 Plot cost-effectiveness frontier ####
plot(df_cea_psa)
ggsave("figs/05a_cea_frontier_psa.png", width = 8, height = 6)

#### 05a.6.4 Cost-effectiveness acceptability curves (CEACs) and frontier (CEAF) ####
ceac_obj <- ceac(wtp = v_wtp, psa = l_psa)
### Regions of highest probability of cost-effectiveness for each strategy
summary(ceac_obj)
### CEAC & CEAF plot
gg_ceac_ceaf <- plot(ceac_obj, txtsize = 16, col = "bw") +
  xlab("Cost-effectiveness threshold (Thousand $/QALY)") +
  guides(color = guide_legend(nrow = 2),
         shape = guide_legend(nrow = 2)) +
  theme(legend.position = "bottom")
gg_ceac_ceaf
ggsave(plot = gg_ceac_ceaf, "figs/05a_ceac_ceaf.png", width = 8, height = 6)
ggsave(plot = gg_ceac_ceaf, "figs/manuscript/fig05a_ceac_ceaf.png", width = 8, height = 6)
ggsave(plot = gg_ceac_ceaf, "figs/manuscript/fig05a_ceac_ceaf.pdf", width = 8, height = 6)
ggsave(plot = gg_ceac_ceaf, "figs/manuscript/fig05a_ceac_ceaf.tiff", width = 8, height = 6)

#### 05a.6.3 Expected Loss Curves (ELCs) ####
## Population affected by the decision
# Source: Siegel RL, Miller KD, Jemal A: Cancer statistics, 2019. CA Cancer J Clin 69:7–34, 2019
pop_evi <- sum(101000*0.33*(1/(1+0.03)^(1:10))) # Annual cases * proportion that are Stage II * ten years

elc_obj <- calc_exp_loss(wtp = v_wtp, psa = l_psa)
elc_obj
gg_elc <- plot(elc_obj, log_y = TRUE, #col = "bw",
               txtsize = 16, n_x_ticks = 8, n_y_ticks = 6) +
  xlab("Cost-effectiveness threshold (Thousand $/QALY)") +
  scale_y_continuous(name = "Expected Loss (Million $)",
                     breaks = c(0, 50,100, 200, 500, 1000, 2000, 5000), 
                     trans = "log", labels = function(x){formatC((x*pop_evi)/1000000, digits = 0, format = "f")}) + 
  guides(color = guide_legend(nrow = 2),
         shape = guide_legend(nrow = 2)) +
  theme(legend.position = "bottom")
gg_elc
ggsave(plot = gg_elc, "figs/05a_elc.png", width = 8, height = 6)
ggsave(plot = gg_elc, "figs/manuscript/fig05b_elc.png", width = 8, height = 6)
ggsave(plot = gg_elc, "figs/manuscript/fig05b_elc.pdf", width = 8, height = 6)
ggsave(plot = gg_elc, "figs/manuscript/fig05b_elc.tiff", width = 8, height = 6)

#### 05a.7 Create linear regression metamodeling sensitivity analysis graphs ####
#### 05a.7.1 One-way sensitivity analysis (OWSA) ####
owsa_lrm_nmb <- owsa(l_psa, params = c("c_Test", "p_HS1", "u_S1", "u_Trt"),
                     outcome = "nmb", wtp = 150000)
plot(owsa_lrm_nmb, txtsize = 16, n_x_ticks = 5, 
     facet_scales = "free") +
  theme(legend.position = "bottom")
ggsave("figs/05a_owsa_lrm_nmb.png", width = 10, height = 6)  

#### 05a.7.2 Optimal strategy with OWSA ####
owsa_opt_strat(owsa = owsa_lrm_nmb)
ggsave("figs/05a_optimal_owsa_lrm_nmb.png", width = 8, height = 6)

#### 05a.7.3 Tornado plot ####
owsa_tornado(owsa = owsa_lrm_nmb)
ggsave("figs/05a_tornado_lrm_Treatment_nmb.png", width = 8, height = 6)

#### 05a.7.4 Two-way sensitivity analysis (TWSA) ####
twsa_lrm_nmb <- twsa(l_psa, param1 = "u_S1", param2 = "u_Trt",
                     outcome = "nmb", wtp = 150000)
plot(twsa_lrm_nmb)
ggsave("figs/05a_twsa_lrm_uS1_uTrt_nmb.png", width = 8, height = 6)  
