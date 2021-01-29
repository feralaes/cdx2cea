################################################################################ 
# This script runs the cost-effectiveness analysis of a hypothetical treatment #
# for the simulated cohort of the Sick-Sicker state-transition model (STM)     #
#                                                                              #                                                                          # 
# Authors:                                                                     #
#     - Fernando Alarid-Escudero, PhD, <fernando.alarid@cide.edu>              # 
#     - Eline Krijkamp, MS                                                     #
#     - Petros Pechlivanoglou, PhD                                             #
#     - Hawre Jalal, MD, PhD                                                   #
#     - Eva A. Enns, PhD                                                       # 
################################################################################
# The structure of this code is according to the DARTH framework               #
# https://github.com/DARTH-git/Decision-Modeling-Framework                     #
################################################################################

# rm(list = ls()) # to clean the workspace

#### 05b.1 Load packages and functions ####
#### 05b.1.1 Load packages ####
# devtools::install_github("DARTH-git/dampack") # Uncomment if dampack not installed
library(dampack) 
library(dplyr)
library(reshape2)

#### 05b.1.2 Load inputs ####
l_params_all <- load_all_params() # function in cdx2cea

#### 05b.1.3 Load functions ####
# no required functions

#### 05b.1.4 Load base-case parameters ####
### Load MAP as best calibration parameter set for deterministic analysis
data("v_calib_post_map")

### Parameters for base-case
## Update base-case parameters with calibrated values at MAP
# l_params_basecase <- update_param_list(l_params_all, v_calib_post_map)


v_calib_params_test <- c(r_DieMets = 0.03870286,
                         r_RecurCDX2pos = 0.003328773,
                         hr_RecurCDX2neg = 3.601069,
                         p_Mets = 0.9808406)
l_params_test <- update_param_list(l_params_all, v_calib_params_test)
p_CDX2neg_init <- 0.07174888
list2env(l_params_test, envir = .GlobalEnv)
l_params_basecase <- l_params_test


#### 05b.2 Time spent in each state ####
#### 05b.2.1 CDX2-negative patients ####
### Run model with treatment
l_cdx2neg_trt <- ce_model(l_params_all = l_params_basecase, 
                          p_CDX2neg_init = 1, Trt = TRUE)
### Run model without treatment
l_cdx2neg_notrt <- ce_model(l_params_all = l_params_basecase, 
                            p_CDX2neg_init = 1, Trt = FALSE)
## Extract life years
ly_cdx2neg_trt   <- l_cdx2neg_trt$tot_ly_und
ly_cdx2neg_notrt <- l_cdx2neg_notrt$tot_ly_und
## Extract cohort trace
m_M_cdx2neg_trt   <- l_cdx2neg_trt$m_M
m_M_cdx2neg_notrt <- l_cdx2neg_notrt$m_M

### Life years gained under treatment 
ly_gain_cdx2neg <- (ly_cdx2neg_trt - ly_cdx2neg_notrt)/ly_cdx2neg_notrt ### RESULT!

### Time spent without and with recurrence
df_ly_rec_cdx2neg <- data.frame(`Health State` = ordered(c("Without recurrence", "With recurrence"),
                                                         c("With recurrence", "Without recurrence")),
                            `No FOLFOX` = c(sum(m_M_cdx2neg_notrt[, 1:2]), 
                                            sum(m_M_cdx2neg_notrt[, 4]))/12,
                            `FOLFOX`    = c(sum(m_M_cdx2neg_trt[, 1:2]), 
                                            sum(m_M_cdx2neg_trt[, 4]))/12,
                            check.names = FALSE
)
df_ly_rec_cdx2neg_lng <- melt(df_ly_rec_cdx2neg, 
                              id.vars = "Health State", 
                              variable.name = "Strategy", 
                              value.name = "Years")
df_charts_ly_rec_cdx2neg <- df_ly_rec_cdx2neg_lng %>% 
  group_by(Strategy) %>%
  mutate(pos = cumsum(Years) - (0.5 * Years))
df_charts_ly_rec_cdx2neg

#### 05b.2.2 CDX2-positive patients ####
### Run model with treatment
l_cdx2pos_trt <- ce_model(l_params_all = l_params_basecase, 
                          p_CDX2neg_init = 0, Trt = TRUE)
### Run model without treatment
l_cdx2pos_notrt <- ce_model(l_params_all = l_params_basecase, 
                            p_CDX2neg_init = 0, Trt = FALSE)
## Extract life years
ly_cdx2pos_trt   <- l_cdx2pos_trt$tot_ly_und
ly_cdx2pos_notrt <- l_cdx2pos_notrt$tot_ly_und
## Extract cohort trace
m_M_cdx2pos_trt   <- l_cdx2pos_trt$m_M
m_M_cdx2pos_notrt <- l_cdx2pos_notrt$m_M

### Life years gained under treatment 
ly_gain_cdx2pos <- (ly_cdx2pos_trt - ly_cdx2pos_notrt)/ly_cdx2pos_notrt ### RESULT!

### Time spent without and with recurrence
df_ly_rec_cdx2pos <- data.frame(`Health State` = ordered(c("Without recurrence", "With recurrence"),
                                                         c("With recurrence", "Without recurrence")),
                                `No FOLFOX` = c(sum(m_M_cdx2pos_notrt[, 1:2]), 
                                                sum(m_M_cdx2pos_notrt[, 4]))/12,
                                `FOLFOX`    = c(sum(m_M_cdx2pos_trt[, 1:2]), 
                                                sum(m_M_cdx2pos_trt[, 4]))/12,
                                check.names = FALSE
)
df_ly_rec_cdx2pos_lng <- melt(df_ly_rec_cdx2pos, 
                              id.vars = "Health State", 
                              variable.name = "Strategy", 
                              value.name = "Years")
df_charts_ly_rec_cdx2pos <- df_ly_rec_cdx2pos_lng %>% 
  group_by(Strategy) %>%
  mutate(pos = cumsum(Years) - (0.5 * Years))
df_charts_ly_rec_cdx2pos

#### 05b.2.3 For both types of patients ####
df_ly_rec_all <- rbind(
  cbind(Group = paste0("CDX2-negative patients (", 
                       scales::percent(round(l_params_basecase$p_CDX2neg, 3), accuracy = 0.1), ")"), 
        df_ly_rec_cdx2neg_lng),
  cbind(Group = paste0("CDX2-positive patients (", 
                       scales::percent(round(1-l_params_basecase$p_CDX2neg, 3), accuracy = 0.1), ")"), 
        df_ly_rec_cdx2pos_lng))
df_charts_ly_rec_all <- rbind(
  cbind(Group = paste0("CDX2-negative patients (", 
                       scales::percent(round(l_params_basecase$p_CDX2neg, 3), accuracy = 0.1), ")"), 
        data.frame(df_charts_ly_rec_cdx2neg, check.names = F)),
  cbind(Group = paste0("CDX2-positive patients (", 
                       scales::percent(round(1-l_params_basecase$p_CDX2neg, 3), accuracy = 0.1), ")"), 
        data.frame(df_charts_ly_rec_cdx2pos, check.names = F)))

ggplot(df_ly_rec_all, 
       aes(x = Strategy, y = Years, fill = `Health State`)) +
  facet_wrap(~ Group, ncol = 1) +
  geom_bar(stat="identity", position = "stack") +
  geom_text(data = df_charts_ly_rec_all, 
            aes(x = Strategy, y = pos, 
                label = format(round(Years, 2), nsmall = 2), group = `Health State`),
            size=5) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
  # scale_fill_grey(start=0.8, end=0.5) +
  xlab("") + 
  scale_y_continuous("Life expectancy (years)", 
                     breaks = dampack::number_ticks(8)) +
  coord_flip() +
  guides(fill = guide_legend(title = "", reverse=T)) +
  theme_bw(base_size = 17) + 
  theme(legend.position = "bottom") +
  NULL
ggsave(filename = "figs/05b_time-states-all-patients.pdf", 
       width = 10, height = 7)
ggsave(filename = "figs/05b_time-states-all-patients.png", 
       width = 10, height = 7)

#### 05b.3 Cost-effectiveness analysis parameters ####
### Strategy names
v_names_str <- l_params_all$v_names_str
### Number of strategies
n_str <- length(v_names_str)

#### 05b.3 Compute cost-effectiveness outcomes ####
df_out_ce <- calculate_ce_out(l_params_all = l_params_basecase, 
                              n_wtp = 150000)
df_out_ce

#### 05b.4 Conduct CEA with deterministic output ####
### Calculate incremental cost-effectiveness ratios (ICERs)
df_cea_det <- dampack::calculate_icers(cost       = df_out_ce$Cost, 
                                       effect     = df_out_ce$Effect, 
                                       strategies = v_names_str)
df_cea_det
### Save CEA table with ICERs
## As .RData
save(df_cea_det,
     file = "tables/05b_deterministic_cea_results.RData")
## As .csv
write.csv(df_cea_det, 
          file = "tables/05b_deterministic_cea_results.csv")

#### 05b.5 Plot cost-effectiveness frontier ####
plot(df_cea_det)
ggsave("figs/05b_cea_frontier.png", width = 8, height = 6)

#### 05b.6 One-way sensitivity analysis (OWSA) ####
df_owsa_input <- data.frame(pars = c("p_CDX2neg", "hr_Recurr_CDXneg_Rx",
                                     "c_Test", "u_Stg4", "hr_RecurCDX2neg"),
                            min = c(0.015, 0.620, 94,  0.2, 1.69),
                            max = c(0.150, 0.970, 179, 0.7, 4.38))
df_owsa_det <- run_owsa_det(params_range = df_owsa_input,
                            params_basecase = l_params_basecase,
                            FUN = calculate_ce_out, # Function to compute outputs
                            outcomes = "ICER",      # Output to do the OWSA on
                            strategies = v_names_str
                            )
df_owsa_icer$parameter <- ordered(df_owsa_icer$parameter, 
                                  # unique(df_owsa_icer$parameter),
                                  labels = c("Cost of test ($)",
                                             "Increased recurrence in CDX2-negative patients",
                                             "Effectiveness of FOLFOX in CDX2−negative patients as a HR",
                                             "Proportion of CDX2−negative patients",
                                             "Utility of metastatic recurrence"))
gg_owsa <- plot(df_owsa_icer, 
                txtsize = 16, n_x_ticks = 5, 
                facet_ncol = 2,
                facet_scales = "free") +
  ylab("$/QALY") +
  theme(legend.position = "")
gg_owsa
ggsave(plot = gg_owsa, 
       filename = "figs/05b_owsa_icer.png", 
       width = 10, height = 8)
ggsave(plot = gg_owsa, 
       filename = "figs/05b_owsa_icer.pdf", 
       width = 10, height = 8)

#### 05b.7 Two-way sensitivity analysis (TWSA) ####

#### 05b.7.1 Proportion of CDX2-negative vs Effectiveness of FOLFOX in CDX2-negative patients (HR) ####
### $50K/QALY
df_twsa_icer_pCDX2_vs_hrCDX2negtrt_50k <- twsa_det(parm1 = "p_CDX2neg",  # parameter 1 name
                     parm2 = "hr_Recurr_CDXneg_Rx", # parameter 2 name
                     ranges = list("p_CDX2neg" = c(0.015, 0.150),
                                   "hr_Recurr_CDXneg_Rx" = c(0.620, 0.970)),
                     nsamps = 40, # number of values  
                     FUN = calculate_ce_out, # Function to compute outputs 
                     params_basecase = l_params_basecase, # Vector with base-case parameters
                     outcome = "NMB",      # Output to do the OWSA on
                     strategies = v_names_str, # Names of strategies
                     n_wtp = 50000        # Extra argument to pass to FUN
)
levels(df_twsa_icer_pCDX2_vs_hrCDX2negtrt_50k$strategy) <- c("(1) No CDX2 testing and no FOLFOX",
                                                              "(2) CDX2 testing and FOLFOX if CDX2-negative")
colnames(df_twsa_icer_pCDX2_vs_hrCDX2negtrt_50k)[1:2] <- c("Proportion of CDX2-neative patients",
                                                            "Effectiveness of FOLFOX in CDX2-negative patients (HR)")
gg_twsa_icer_pCDX2_vs_hrCDX2negtrt_50k <- plot(df_twsa_icer_pCDX2_vs_hrCDX2negtrt_50k,
                                               col = c("bw")) +
  theme(legend.position = "bottom")
gg_twsa_icer_pCDX2_vs_hrCDX2negtrt_50k 
ggsave(plot = gg_twsa_icer_pCDX2_vs_hrCDX2negtrt_50k,
       filename = "figs/05b_twsa_p_CDX2neg_hr_RecurrCDXnegRx_nmb_50k.png", 
       width = 8, height = 6)
ggsave(plot = gg_twsa_icer_pCDX2_vs_hrCDX2negtrt_50k,
       filename = "figs/05b_twsa_p_CDX2neg_hr_RecurrCDXnegRx_nmb_50k.pdf", 
       width = 8, height = 6)

### $100K/QALY
df_twsa_icer_pCDX2_vs_hrCDX2negtrt_100k <- twsa_det(parm1 = "p_CDX2neg",  # parameter 1 name
                                                    parm2 = "hr_Recurr_CDXneg_Rx", # parameter 2 name
                                                    ranges = list("p_CDX2neg" = c(0.015, 0.150),
                                                                  "hr_Recurr_CDXneg_Rx" = c(0.620, 0.970)),
                                                    nsamps = 40, # number of values  
                                                    FUN = calculate_ce_out, # Function to compute outputs 
                                                    params_basecase = l_params_basecase, # Vector with base-case parameters
                                                    outcome = "NMB",      # Output to do the OWSA on
                                                    strategies = v_names_str, # Names of strategies
                                                    n_wtp = 100000        # Extra argument to pass to FUN
)
levels(df_twsa_icer_pCDX2_vs_hrCDX2negtrt_100k$strategy) <- c("(1) No CDX2 testing and no FOLFOX",
                                                              "(2) CDX2 testing and FOLFOX if CDX2-negative")
colnames(df_twsa_icer_pCDX2_vs_hrCDX2negtrt_100k)[1:2] <- c("Proportion of CDX2-neative patients",
                                                            "Effectiveness of FOLFOX in CDX2-negative patients (HR)")
gg_twsa_icer_pCDX2_vs_hrCDX2negtrt_100k <- plot(df_twsa_icer_pCDX2_vs_hrCDX2negtrt_100k,
                                                col = c("bw")) +
  theme(legend.position = "bottom")
gg_twsa_icer_pCDX2_vs_hrCDX2negtrt_100k 
ggsave(plot = gg_twsa_icer_pCDX2_vs_hrCDX2negtrt_100k,
       filename = "figs/05b_twsa_p_CDX2neg_hr_RecurrCDXnegRx_nmb_100k.png", 
       width = 8, height = 6)
ggsave(plot = gg_twsa_icer_pCDX2_vs_hrCDX2negtrt_100k,
       filename = "figs/05b_twsa_p_CDX2neg_hr_RecurrCDXnegRx_nmb_100k.pdf", 
       width = 8, height = 6)

#### 05b.7.2 Increased recurrence in CDX2-negative patients vs Effectiveness of FOLFOX in CDX2-negative patients (HR) ####
### $50K/QALY
df_twsa_icer_hrRecurCDX2neg_vs_hrCDX2negtrt_50k <- twsa_det(parm1 = "hr_RecurCDX2neg",  # parameter 1 name
                                                   parm2 = "hr_Recurr_CDXneg_Rx", # parameter 2 name
                                                   ranges = list("hr_RecurCDX2neg" = c(1.69, 4.38),
                                                                 "hr_Recurr_CDXneg_Rx" = c(0.620, 0.970)),
                                                   nsamps = 40, # number of values  
                                                   FUN = calculate_ce_out, # Function to compute outputs 
                                                   params_basecase = l_params_basecase, # Vector with base-case parameters
                                                   outcome = "NMB",      # Output to do the OWSA on
                                                   strategies = v_names_str, # Names of strategies
                                                   n_wtp = 50000        # Extra argument to pass to FUN
)
levels(df_twsa_icer_hrRecurCDX2neg_vs_hrCDX2negtrt_50k$strategy) <- c("(1) No CDX2 testing and no FOLFOX",
                                                                      "(2) CDX2 testing and FOLFOX if CDX2-negative")
colnames(df_twsa_icer_hrRecurCDX2neg_vs_hrCDX2negtrt_50k)[1:2] <- c("Increased recurrence in CDX2-negative patients",
                                                                    "Effectiveness of FOLFOX in CDX2-negative patients (HR)")
gg_twsa_icer_hrRecurCDX2neg_vs_hrCDX2negtrt_50k <- plot(df_twsa_icer_hrRecurCDX2neg_vs_hrCDX2negtrt_50k,
                                                        col = c("bw")) +
  theme(legend.position = "bottom")
gg_twsa_icer_hrRecurCDX2neg_vs_hrCDX2negtrt_50k 
ggsave(plot = gg_twsa_icer_hrRecurCDX2neg_vs_hrCDX2negtrt_50k,
       filename = "figs/05b_twsa_hr_Recurr_CDXneg_Rx_RecurrCDXnegRx_nmb_50k.png", 
       width = 8, height = 6)
ggsave(plot = gg_twsa_icer_hrRecurCDX2neg_vs_hrCDX2negtrt_50k,
       filename = "figs/05b_twsa_hr_Recurr_CDXneg_Rx_RecurrCDXnegRx_nmb_50k.pdf", 
       width = 8, height = 6)

### $100K/QALY
df_twsa_icer_hrRecurCDX2neg_vs_hrCDX2negtrt_100k <- twsa_det(parm1 = "hr_RecurCDX2neg",  # parameter 1 name
                                                            parm2 = "hr_Recurr_CDXneg_Rx", # parameter 2 name
                                                            ranges = list("hr_RecurCDX2neg" = c(1.69, 4.38),
                                                                          "hr_Recurr_CDXneg_Rx" = c(0.620, 0.970)),
                                                            nsamps = 40, # number of values  
                                                            FUN = calculate_ce_out, # Function to compute outputs 
                                                            params_basecase = l_params_basecase, # Vector with base-case parameters
                                                            outcome = "NMB",      # Output to do the OWSA on
                                                            strategies = v_names_str, # Names of strategies
                                                            n_wtp = 100000        # Extra argument to pass to FUN
)
levels(df_twsa_icer_hrRecurCDX2neg_vs_hrCDX2negtrt_100k$strategy) <- c("(1) No CDX2 testing and no FOLFOX",
                                                                       "(2) CDX2 testing and FOLFOX if CDX2-negative")
colnames(df_twsa_icer_hrRecurCDX2neg_vs_hrCDX2negtrt_100k)[1:2] <- c("Increased recurrence in CDX2-negative patients",
                                                                     "Effectiveness of FOLFOX in CDX2-negative patients (HR)")
gg_twsa_icer_hrRecurCDX2neg_vs_hrCDX2negtrt_100k <- plot(df_twsa_icer_hrRecurCDX2neg_vs_hrCDX2negtrt_100k,
                                                         col = c("bw")) +
  theme(legend.position = "bottom")
gg_twsa_icer_hrRecurCDX2neg_vs_hrCDX2negtrt_100k 
ggsave(plot = gg_twsa_icer_hrRecurCDX2neg_vs_hrCDX2negtrt_100k,
       filename = "figs/05b_twsa_hr_Recurr_CDXneg_Rx_RecurrCDXnegRx_nmb_100k.png", 
       width = 8, height = 6)
ggsave(plot = gg_twsa_icer_hrRecurCDX2neg_vs_hrCDX2negtrt_100k,
       filename = "figs/05b_twsa_hr_Recurr_CDXneg_Rx_RecurrCDXnegRx_nmb_100k.pdf", 
       width = 8, height = 6)
