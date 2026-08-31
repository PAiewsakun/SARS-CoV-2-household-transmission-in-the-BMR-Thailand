# ###############################
# Load global variables
# ###############################
source(file.path("miscellaneous.R"))

# ###############################
# Load required libraries
# ###############################
library(lme4) # glmer, glmerControl
library(buildmer)
library(performance) #icc
library(emmeans)

library(tidyverse)
library(gtsummary)

select <- dplyr::select #Sometimes, if you run multiple scripts in one session, you can get "select" from "MASS", which can cause this script to fail. This line is to make sure that "select" is from "dplyr".

# ###############################
# Prepare participant metadata
# ###############################
# Load data, format and select variables
# -------------------------------
participant_metadata <- read.delim(participant_metadata_filepath, header = TRUE, sep = "\t") %>%
  mutate(
    age_cat = factor(age_cat, c("Child", "Adult")),
    family_id = factor(family_id), # for joining with the HH_metadata to add first_HH_case_age_cat to the df
    
    overall_res_end = factor(overall_res_end, levels = c("Neg", "Pos")),
    overall_res_end_with_symp_severity = factor(overall_res_end_with_symp_severity, levels = names(overall_res_end_with_symp_severity_cols)),
    
    sex = factor(sex, c("Female", "Male")),
    BMI = as.numeric(BMI),
    
    pre_existing_med_cond_ht = recode_factor(factor(pre_existing_med_cond_ht), "0" = "No", "1" = "Yes"),
    pre_existing_med_cond_dm = recode_factor(factor(pre_existing_med_cond_dm), "0" = "No", "1" = "Yes"),
    pre_existing_med_cond_al = recode_factor(factor(pre_existing_med_cond_al), "0" = "No", "1" = "Yes"),
    pre_existing_med_cond_hl = recode_factor(factor(pre_existing_med_cond_hl), "0" = "No", "1" = "Yes"),
    
    num_prev_COVID = as.integer(num_prev_COVID), 
    time_to_last_prev_COVID_cat = factor(time_to_last_prev_COVID_cat, levels = c("Never", ">1 yr", "0.5-1 yr", "<0.5 yr")),
    
    num_vac = as.integer(num_vac),
    vac_viral_vector = as.integer(vac_viral_vector),
    vac_mRNA = as.integer(vac_mRNA),
    vac_inact_virus = as.integer(vac_inact_virus),
    vac_other = as.integer(vac_other),
    time_to_last_vac_cat = factor(time_to_last_vac_cat, levels = c("Never", ">1 yr", "0.5-1 yr", "<0.5 yr")),
    
    major_lineage = ifelse(major_lineage %in% c("XBB", "XBL"), "XBB or XBL", major_lineage),
    major_lineage = factor(major_lineage, levels = names(major_lineage_cols)), 
    
    wave = case_when(
      (enrolment_epiyear == 2022) | (enrolment_epiyear == 2023 & enrolment_epiweek <= 10) ~ "BA.5 and BA.2.75 wave",
      enrolment_epiyear == 2023 & enrolment_epiweek > 10 ~ "XBB and XBL wave",
      enrolment_epiyear == 2024 ~ "BA.2.86 wave",
      TRUE ~ NA
    ),
    wave = factor(wave, levels = c("BA.5 and BA.2.75 wave", "XBB and XBL wave", "BA.2.86 wave")),
  ) %>%
  select(
    family_id, # for joining first_HH_case_age_cat from HH_metadata df
    
    age_cat, 
    
    overall_res_end,
    overall_res_end_with_symp_severity,
    
    sex, BMI,
    
    pre_existing_med_cond_ht, pre_existing_med_cond_dm, pre_existing_med_cond_al, pre_existing_med_cond_hl,
    
    num_prev_COVID,
    time_to_last_prev_COVID_cat,
    
    num_vac, 
    vac_viral_vector, vac_mRNA, vac_inact_virus, vac_other,
    time_to_last_vac_cat,
    
    pcr_orf1ab_ct, pcr_n_ct, pcr_e_ct,
    
    wave,
    major_lineage,
  )

# Add first_HH_case_age_cat to the df
# -------------------------------
first_HH_case_age_cat_dat <- read.delim(HH_metadata_filepath, header = TRUE, sep = "\t") %>%
  select(family_id, first_HH_case_age_cat) %>% unique %>%
  mutate(
    first_HH_case_age_cat = ifelse(first_HH_case_age_cat == "Adult & child", "Adult", first_HH_case_age_cat),
    first_HH_case_age_cat = recode_factor(factor(first_HH_case_age_cat), "Child" = "Child-first household", "Adult" = "Adult-first household"),
    first_HH_case_age_cat = factor(first_HH_case_age_cat, c("Child-first household", "Adult-first household"))
  )

participant_metadata <- left_join(
  participant_metadata,
  first_HH_case_age_cat_dat,
  by = "family_id"
) %>% 
  relocate(first_HH_case_age_cat)

# ###############################
# Identify factors associated with endpoint SARS-CoV-2 positivity 
# using multivariable mixed-effects logistic regression analysis 
# accounting for household-level clustering, with backward BIC-based model selection.
# ###############################
# The initial model included additive fixed effects of all variables + an interaction between first_HH_case_age_cat * age_cat.
# The model doesn't include the total number of vaccine doses received,...
# as it is the sum of the numbers of doses received from individual vaccine platforms, and they were all included in the initial model. 
# -------------------------------
# Prepare data for model fitting (N = 210)
endpoint_case_rate_mdl_full_dat <- participant_metadata %>% 
  select(
    overall_res_end,
    
    first_HH_case_age_cat, age_cat,
    
    sex,
    BMI,
    pre_existing_med_cond_ht, pre_existing_med_cond_dm, pre_existing_med_cond_al, pre_existing_med_cond_hl,
    
    num_prev_COVID,
    time_to_last_prev_COVID_cat,
    
    # num_vac,
    vac_viral_vector, vac_mRNA, vac_inact_virus, vac_other, 
    time_to_last_vac_cat,
    
    wave,
    
    family_id
  ) %>% na.omit()

nrow(endpoint_case_rate_mdl_full_dat)
# [1] 210

# Model fitting + stepwise BIC-based model selection
endpoint_case_rate_mdl_reduced <- buildmer(
  overall_res_end ~
    first_HH_case_age_cat * age_cat +
    
    sex +
    BMI +
    pre_existing_med_cond_ht + pre_existing_med_cond_dm + pre_existing_med_cond_al + pre_existing_med_cond_hl +
    
    num_prev_COVID +
    time_to_last_prev_COVID_cat +
    
    vac_viral_vector + vac_mRNA + vac_inact_virus + vac_other +
    time_to_last_vac_cat +
    
    wave +
    
    (1 | family_id),
  data = endpoint_case_rate_mdl_full_dat,
  family = binomial,
  buildmerControl = buildmerControl(
    include=~ (1 | family_id),
    direction = "backward",
    crit = "BIC",
    args = list(
      control = glmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 2e5)
      ) # end glmerControl
    ) # end args list
  ) # end buildmerControl
) # end buildmer

# Examine the reduced model
summary(endpoint_case_rate_mdl_reduced@model)
# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
#  Family: binomial  ( logit )
# Formula: overall_res_end ~ 1 + first_HH_case_age_cat + age_cat + (1 |      family_id)
#    Data: endpoint_case_rate_mdl_full_dat
# Control: structure(list(optimizer = c("bobyqa", "bobyqa"), restart_edge = FALSE,  
#     boundary.tol = 1e-05, calc.derivs = NULL, use.last.params = FALSE,  
#     checkControl = list(autoscale = NULL, check.nobs.vs.rankZ = "ignore",          check.nobs.vs.nlev = "stop", check.nlev.gtreq.5 = "ignore",  
#         check.nlev.gtr.1 = "stop", check.nobs.vs.nRE = "stop",          check.rankX = "message+drop.cols", check.scaleX = "warning",  
#         check.formula.LHS = "stop", check.response.not.const = "stop"),  
#     checkConv = list(check.conv.nobsmax = 10000, check.conv.nparmax = 20,          check.conv.grad = list(action = "warning", tol = 0.002,  
#             relTol = NULL), check.conv.singular = list(action = "message",  
#             tol = 1e-04), check.conv.hess = list(action = "warning",  
#             tol = 1e-06)), optCtrl = list(maxfun = 2e+05), tolPwrss = 1e-07,  
#     compDev = TRUE, nAGQ0initStep = TRUE), class = c("glmerControl",  "merControl"))
# 
#       AIC       BIC    logLik -2*log(L)  df.resid 
#     211.9     225.2    -101.9     203.9       206 
# 
# Scaled residuals: 
#     Min      1Q  Median      3Q     Max 
# -3.4695 -0.4460  0.1736  0.4668  2.2352 
# 
# Random effects:
#  Groups    Name        Variance Std.Dev.
#  family_id (Intercept) 0.7119   0.8437  
# Number of obs: 210, groups:  family_id, 60
# 
# Fixed effects:
#                                            Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                                  1.6326     0.4322   3.778 0.000158 ***
# first_HH_case_age_catAdult-first household   1.5568     0.4535   3.433 0.000597 ***
# age_catAdult                                -3.1385     0.5228  -6.004 1.93e-09 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#             (Intr) f_HH_h
# fr_HH___A-h -0.215       
# age_catAdlt -0.747 -0.277

icc(endpoint_case_rate_mdl_reduced@model)
# # Intraclass Correlation Coefficient
# 
#     Adjusted ICC: 0.178
#   Unadjusted ICC: 0.099

car::Anova(endpoint_case_rate_mdl_reduced@model)
# Analysis of Deviance Table (Type II Wald chisquare tests)
# 
# Response: overall_res_end
#                        Chisq Df Pr(>Chisq)    
# first_HH_case_age_cat 11.786  1  0.0005967 ***
# age_cat               36.045  1  1.928e-09 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Re-evaluate the reduced model using a larger dataset
# -------------------------------
# Prepare data for model fitting (N = 265)
endpoint_case_rate_mdl_full_dat <- participant_metadata %>% 
  select(
    overall_res_end,
    
    first_HH_case_age_cat, age_cat,
    
    family_id
    ) %>% na.omit()

nrow(endpoint_case_rate_mdl_full_dat)
# [1] 265

# Model fitting + stepwise BIC-based model selection
endpoint_case_rate_mdl_reduced <- buildmer(
  overall_res_end ~ first_HH_case_age_cat + age_cat + (1 | family_id),
  # overall_res_end ~ first_HH_case_age_cat * age_cat + (1 | family_id), #give the same result as above
  data = endpoint_case_rate_mdl_full_dat,
  family = binomial,
  buildmerControl = buildmerControl(
    include=~ (1 | family_id),
    direction = "backward",
    crit = "BIC",
    args = list(
      control = glmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 2e5)
      ) # end glmerControl
    ) # end args list
  ) # end buildmerControl
) # end buildmer

# Examine the reduced model
summary(endpoint_case_rate_mdl_reduced@model)
# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
#  Family: binomial  ( logit )
# Formula: overall_res_end ~ 1 + first_HH_case_age_cat + age_cat + (1 |      family_id)
#    Data: endpoint_case_rate_mdl_full_dat
# Control: structure(list(optimizer = c("bobyqa", "bobyqa"), restart_edge = FALSE,  
#     boundary.tol = 1e-05, calc.derivs = NULL, use.last.params = FALSE,  
#     checkControl = list(autoscale = NULL, check.nobs.vs.rankZ = "ignore",          check.nobs.vs.nlev = "stop", check.nlev.gtreq.5 = "ignore",  
#         check.nlev.gtr.1 = "stop", check.nobs.vs.nRE = "stop",          check.rankX = "message+drop.cols", check.scaleX = "warning",  
#         check.formula.LHS = "stop", check.response.not.const = "stop"),  
#     checkConv = list(check.conv.nobsmax = 10000, check.conv.nparmax = 20,          check.conv.grad = list(action = "warning", tol = 0.002,  
#             relTol = NULL), check.conv.singular = list(action = "message",  
#             tol = 1e-04), check.conv.hess = list(action = "warning",  
#             tol = 1e-06)), optCtrl = list(maxfun = 2e+05), tolPwrss = 1e-07,  
#     compDev = TRUE, nAGQ0initStep = TRUE), class = c("glmerControl",  "merControl"))
# 
#       AIC       BIC    logLik -2*log(L)  df.resid 
#     274.3     288.6    -133.2     266.3       261 
# 
# Scaled residuals: 
#     Min      1Q  Median      3Q     Max 
# -4.9088 -0.5433  0.1931  0.5543  1.7715 
# 
# Random effects:
#  Groups    Name        Variance Std.Dev.
#  family_id (Intercept) 0.4643   0.6814  
# Number of obs: 265, groups:  family_id, 60
# 
# Fixed effects:
#                                            Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                                  1.4012     0.3780   3.707  0.00021 ***
# first_HH_case_age_catAdult-first household   2.0847     0.4016   5.192 2.08e-07 ***
# age_catAdult                                -2.4121     0.4381  -5.505 3.69e-08 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#             (Intr) f_HH_h
# fr_HH___A-h -0.123       
# age_catAdlt -0.762 -0.364

icc(endpoint_case_rate_mdl_reduced@model)
# # Intraclass Correlation Coefficient
# 
#     Adjusted ICC: 0.124
#   Unadjusted ICC: 0.076

car::Anova(endpoint_case_rate_mdl_reduced@model)
# Analysis of Deviance Table (Type II Wald chisquare tests)
# 
# Response: overall_res_end
#                        Chisq Df Pr(>Chisq)    
# first_HH_case_age_cat 26.953  1  2.085e-07 ***
# age_cat               30.306  1  3.689e-08 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Fit the minimal model: Endpoint case rate ~ household type + age group + (1|family_id)
# The reduced model is ready the minimal model, so this is it
# -------------------------------
endpoint_case_rate_mdl_minimal <- endpoint_case_rate_mdl_reduced@model 

# ###############################
# Estimate adjusted odds ratios between various participant groups based on the min model
# ###############################
# Between child-first vs adult-first households
# -------------------------------
emmeans(
  endpoint_case_rate_mdl_minimal, 
  pairwise ~ first_HH_case_age_cat,
  type = "response",
) %>% summary(infer = c(TRUE, TRUE)) %>%
  pluck("contrasts") %>% as.data.frame()
#  contrast                                          odds.ratio         SE  df  asymp.LCL asymp.UCL null z.ratio p.value
#  (Child-first household) / (Adult-first household)   0.124343 0.04993039 Inf 0.05660028 0.2731646    1  -5.192 <0.0001
# 
# Results are averaged over the levels of: age_cat 
# Confidence level used: 0.95 
# Intervals are back-transformed from the log odds ratio scale 
# Tests are performed on the log odds ratio scale 

# Between children vs adults
# -------------------------------
emmeans(
  endpoint_case_rate_mdl_minimal, 
  pairwise ~ age_cat,
  type = "response", 
) %>% summary(infer = c(TRUE, TRUE)) %>%
  pluck("contrasts") %>% as.data.frame()
#  contrast      odds.ratio       SE  df asymp.LCL asymp.UCL null z.ratio p.value
#  Child / Adult   11.15686 4.888366 Inf  4.727036  26.33267    1   5.505 <0.0001
# 
# Results are averaged over the levels of: first_HH_case_age_cat 
# Confidence level used: 0.95 
# Intervals are back-transformed from the log odds ratio scale 
# Tests are performed on the log odds ratio scale

# ###############################
# Estimate marginalised endpoint case rates for various participant groups based on the min model
# ###############################
# For first_HH_case_age_cat * age_cat
# -------------------------------
emmeans(
  endpoint_case_rate_mdl_minimal, 
  ~ first_HH_case_age_cat * age_cat, 
  type = "response", weights = "cells"
) %>% summary(infer = c(TRUE, TRUE)) %>% as.data.frame()
#  first_HH_case_age_cat age_cat      prob         SE  df asymp.LCL asymp.UCL null z.ratio p.value
#  Child-first household Child   0.8023677 0.05994387 Inf 0.6593178 0.8949247  0.5   3.707  0.0002
#  Adult-first household Child   0.9702830 0.01489748 Inf 0.9222432 0.9889969  0.5   6.747 <0.0001
#  Child-first household Adult   0.2668046 0.05618145 Inf 0.1716766 0.3898364  0.5  -3.520  0.0004
#  Adult-first household Adult   0.7453218 0.05311787 Inf 0.6283968 0.8351101  0.5   3.837  0.0001
# 
# Confidence level used: 0.95 
# Intervals are back-transformed from the logit scale 
# Tests are performed on the logit scale 

# For child-first and adult-first households
# -------------------------------
emmeans(
  endpoint_case_rate_mdl_minimal, 
  ~ first_HH_case_age_cat, 
  type = "response", weights = "cells"
) %>% summary(infer = c(TRUE, TRUE)) %>% as.data.frame()
#  first_HH_case_age_cat      prob         SE  df asymp.LCL asymp.UCL null z.ratio p.value
#  Child-first household 0.4657442 0.06095825 Inf 0.3503740 0.5848996  0.5  -0.560  0.5753
#  Adult-first household 0.8673621 0.03611135 Inf 0.7794773 0.9236527  0.5   5.982 <0.0001
# 
# Results are averaged over the levels of: age_cat 
# Confidence level used: 0.95 
# Intervals are back-transformed from the logit scale 
# Tests are performed on the logit scale 

# For children and adults
# -------------------------------
emmeans(
  endpoint_case_rate_mdl_minimal, 
  ~ age_cat, 
  type = "response", weights = "cells"
) %>% summary(infer = c(TRUE, TRUE)) %>% as.data.frame()
#  age_cat      prob         SE  df asymp.LCL asymp.UCL null z.ratio p.value
#  Child   0.9200859 0.02983181 Inf 0.8386632 0.9622654  0.5   6.023 <0.0001
#  Adult   0.5244148 0.04995207 Inf 0.4268252 0.6201747  0.5   0.488  0.6256
# 
# Results are averaged over the levels of: first_HH_case_age_cat 
# Confidence level used: 0.95 
# Intervals are back-transformed from the logit scale 
# Tests are performed on the logit scale 

# ###############################
# Create Table 2 to summarise the results
# ###############################
# Set labels for variables in the gtsummary tables
# -------------------------------
label_list <- list(
  first_HH_case_age_cat = "Household type",
  age_cat = "Age group",
  
  overall_res_end = "Endpoint COVID-19 status",
  overall_res_end_with_symp_severity = "Endpoint COVID-19 status with symptom severity",
  
  sex = "Sex",
  BMI = "BMI",
  
  pre_existing_med_cond_ht = "Hypertension",
  pre_existing_med_cond_dm = "Diabetes",
  pre_existing_med_cond_al = "Allergic rhinitis",
  pre_existing_med_cond_hl = "Hyperlipidaemia",
  
  num_prev_COVID = "Number of previous infections",
  time_to_last_prev_COVID_cat = "Time since most recent infection",
  
  num_vac = "Total number of vaccine doses received",
  vac_viral_vector = "Number of viral-vector vaccine doses received",
  vac_mRNA = "Number of mRNA vaccine doses received",
  vac_inact_virus = "Number of inactivated-virus vaccine doses received",
  vac_other = "Number of other vaccine doses received",
  time_to_last_vac_cat = "Time since most recent vaccination",
  
  pcr_orf1ab_ct = "RT-PCR Ct value of ORF1ab gene",
  pcr_n_ct = "RT-PCR Ct value of N gene",
  pcr_e_ct = "RT-PCR Ct value of E gene",
  
  wave = "Virus wave",
  major_lineage = "Virus lineage"
)

# Descriptive statistics
# -------------------------------
Table_2A <- tbl_summary(
  data = participant_metadata,
  by = overall_res_end,
  include = c(
    first_HH_case_age_cat,
    age_cat,
    sex,
    BMI,
    pre_existing_med_cond_ht,
    pre_existing_med_cond_dm,
    pre_existing_med_cond_al,
    pre_existing_med_cond_hl,
    num_prev_COVID,
    time_to_last_prev_COVID_cat,
    num_vac,
    vac_viral_vector, vac_mRNA, vac_inact_virus, vac_other, 
    time_to_last_vac_cat,
    wave
  ),
  label = label_list,
  missing = "ifany", missing_text = "# of participants not included",
  type = list(where(is.numeric) ~ "continuous"),
  statistic = list(
    all_continuous() ~ "{mean} ({sd}; {min}-{max})" 
  ),
  digits = all_continuous() ~ 2,
  percent = "row"
) %>% 
  add_n() %>%
  bold_labels() %>%
  add_p(
    pvalue_fun = label_style_pvalue(digits = 3),
    include = everything()
  ) %>%
  add_q(method = "fdr") %>%
  modify_header(
    all_stat_cols() ~ "**{ifelse(level == 'Neg', 'Endpoint-negative','Endpoint-positive')}**\n(N = {n})",
    label ~ "**Variable**",
    p.value ~ "**p-value**",
    q.value ~ "**FDR-adjusted p-value**"
  )
# 1 missing row in the "overall_res_end" column has been removed. 266 >> 265

# Univariable mixed-effects logistic regression models
# -------------------------------
Table_2B <- tbl_uvregression(
  data = participant_metadata,
  method = glmer,
  y = overall_res_end,
  include = c(
    first_HH_case_age_cat,
    age_cat,
    sex,
    BMI,
    pre_existing_med_cond_ht,
    pre_existing_med_cond_dm,
    pre_existing_med_cond_al,
    pre_existing_med_cond_hl,
    num_prev_COVID,
    time_to_last_prev_COVID_cat,
    num_vac,
    vac_viral_vector, vac_mRNA, vac_inact_virus, vac_other, 
    time_to_last_vac_cat,
    wave
  ),
  label = label_list,
  method.args = list(
    family = binomial,
    control = glmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 2e5)
    )
  ),
  formula = "{y} ~ {x} + (1|family_id)",
  exponentiate = TRUE,
  add_estimate_to_reference_rows = TRUE,
  conf.int = TRUE,
  pvalue_fun = label_style_pvalue(digits = 3),
  hide_n = FALSE,
  show_single_row = c("pre_existing_med_cond_ht", "pre_existing_med_cond_dm", "pre_existing_med_cond_al", "pre_existing_med_cond_hl")
) %>%
  bold_labels() %>%
  add_q(method = "fdr") %>%
  modify_header(
    label ~ "**Variable**",
    estimate ~ "**Crude OR**",
    p.value ~ "**p-value**",
    q.value ~ "**FDR-adjusted p-value**"
  )
# There was a warning constructing the model for variable "vac_other". See message below.
# ! unable to evaluate scaled gradient and Model failed to converge: degenerate Hessian with 1 negative eigenvalues See ?lme4::convergence and ?lme4::troubleshooting.
# There was a warning running `tbl_regression()` for variable "vac_other". See message below.
# ! variance-covariance matrix computed from finite-difference Hessian is not positive definite or contains NA values: falling back to var-cov estimated from RX, variance-covariance matrix computed
#   from finite-difference Hessian is not positive definite or contains NA values: falling back to var-cov estimated from RX, variance-covariance matrix computed from finite-difference Hessian is
#   not positive definite or contains NA values: falling back to var-cov estimated from RX, and variance-covariance matrix computed from finite-difference Hessian is not positive definite or
#   contains NA values: falling back to var-cov estimated from RX

# Most parsimonious mixed-effects logistic regression model
# -------------------------------
Table_2C <- tbl_regression(
  endpoint_case_rate_mdl_minimal,
  exponentiate = TRUE,
  conf.int = TRUE,
  label = label_list,
  show_single_row = c(),
  add_estimate_to_reference_rows = TRUE
) %>%
  add_n() %>% 
  modify_header(
    estimate = "**Adjusted OR**",
    p.value = "**p-value**"
  ) %>%
  bold_labels()

# Combine all tables together
# -------------------------------
Table_2 <- tbl_merge(
  tbls = list(Table_2A, Table_2B, Table_2C),
  tab_spanner = c(
    "**Descriptive statistics**",
    "**Univariable mixed-effects logistic regression**",
    "**Most parsimonious mixed-effects logistic regression model**"
  )
)

Table_2 %>% as.data.frame()
#                                              **Variable** **N** **Endpoint-negative**\n(N = 95) **Endpoint-positive**\n(N = 170) **p-value** **FDR-adjusted p-value** **N** **Crude OR** **95% CI** **p-value** **FDR-adjusted p-value** **N** **Adjusted OR** **95% CI** **p-value**
# 1                                      __Household type__   265                            <NA>                             <NA>      <0.001                   <0.001   265         <NA>       <NA>        <NA>                     <NA>   265            <NA>       <NA>        <NA>
# 2                                   Child-first household  <NA>                        68 (54%)                         59 (46%)        <NA>                     <NA>  <NA>         1.00       <NA>        <NA>                     <NA>  <NA>            1.00       <NA>        <NA>
# 3                                   Adult-first household  <NA>                        27 (20%)                        111 (80%)        <NA>                     <NA>  <NA>         4.94 2.71, 8.99      <0.001                   <0.001  <NA>            8.04 3.66, 17.7      <0.001
# 4                                           __Age group__   265                            <NA>                             <NA>      <0.001                   <0.001   265         <NA>       <NA>        <NA>                     <NA>   265            <NA>       <NA>        <NA>
# 5                                                   Child  <NA>                        12 (13%)                         80 (87%)        <NA>                     <NA>  <NA>         1.00       <NA>        <NA>                     <NA>  <NA>            1.00       <NA>        <NA>
# 6                                                   Adult  <NA>                        83 (48%)                         90 (52%)        <NA>                     <NA>  <NA>         0.09 0.04, 0.22      <0.001                   <0.001  <NA>            0.09 0.04, 0.21      <0.001
# 7                                                 __Sex__   265                            <NA>                             <NA>       0.378                    0.494   265         <NA>       <NA>        <NA>                     <NA>  <NA>            <NA>       <NA>        <NA>
# 8                                                  Female  <NA>                        50 (34%)                         99 (66%)        <NA>                     <NA>  <NA>         1.00       <NA>        <NA>                     <NA>  <NA>            <NA>       <NA>        <NA>
# 9                                                    Male  <NA>                        45 (39%)                         71 (61%)        <NA>                     <NA>  <NA>         0.69 0.39, 1.22       0.200                    0.276  <NA>            <NA>       <NA>        <NA>
# 10                                                __BMI__   214       24.02 (6.39; 11.04-51.94)        20.48 (6.22; 10.18-44.58)      <0.001                   <0.001   214         0.90 0.85, 0.95      <0.001                   <0.001  <NA>            <NA>       <NA>        <NA>
# 11                         # of participants not included  <NA>                               0                               51        <NA>                     <NA>  <NA>         <NA>       <NA>        <NA>                     <NA>  <NA>            <NA>       <NA>        <NA>
# 12                                       __Hypertension__   220                        18 (72%)                          7 (28%)       0.002                    0.005   220         0.22 0.08, 0.62       0.004                    0.010  <NA>            <NA>       <NA>        <NA>
# 13                         # of participants not included  <NA>                               0                               45        <NA>                     <NA>  <NA>         <NA>       <NA>        <NA>                     <NA>  <NA>            <NA>       <NA>        <NA>
# 14                                           __Diabetes__   220                         7 (64%)                          4 (36%)       0.214                    0.303   220         0.38 0.10, 1.50       0.166                    0.257  <NA>            <NA>       <NA>        <NA>
# 15                         # of participants not included  <NA>                               0                               45        <NA>                     <NA>  <NA>         <NA>       <NA>        <NA>                     <NA>  <NA>            <NA>       <NA>        <NA>
# 16                                  __Allergic rhinitis__   220                        11 (37%)                         19 (63%)       0.438                    0.521   220         1.38 0.58, 3.29       0.461                    0.564  <NA>            <NA>       <NA>        <NA>
# 17                         # of participants not included  <NA>                               0                               45        <NA>                     <NA>  <NA>         <NA>       <NA>        <NA>                     <NA>  <NA>            <NA>       <NA>        <NA>
# 18                                    __Hyperlipidaemia__   220                        15 (75%)                          5 (25%)       0.003                    0.005   220         0.20 0.06, 0.63       0.006                    0.012  <NA>            <NA>       <NA>        <NA>
# 19                         # of participants not included  <NA>                               0                               45        <NA>                     <NA>  <NA>         <NA>       <NA>        <NA>                     <NA>  <NA>            <NA>       <NA>        <NA>
# 20                      __Number of previous infections__   265          0.84 (0.59; 0.00-2.00)           0.49 (0.59; 0.00-2.00)      <0.001                   <0.001   265         0.34 0.20, 0.57      <0.001                   <0.001  <NA>            <NA>       <NA>        <NA>
# 21                   __Time since most recent infection__   262                            <NA>                             <NA>      <0.001                   <0.001   262         <NA>       <NA>        <NA>                     <NA>  <NA>            <NA>       <NA>        <NA>
# 22                                                  Never  <NA>                        25 (21%)                         95 (79%)        <NA>                     <NA>  <NA>         1.00       <NA>        <NA>                     <NA>  <NA>            <NA>       <NA>        <NA>
# 23                                                  >1 yr  <NA>                        33 (40%)                         49 (60%)        <NA>                     <NA>  <NA>         0.38 0.18, 0.81       0.012                    0.021  <NA>            <NA>       <NA>        <NA>
# 24                                               0.5-1 yr  <NA>                        27 (54%)                         23 (46%)        <NA>                     <NA>  <NA>         0.19 0.09, 0.44      <0.001                   <0.001  <NA>            <NA>       <NA>        <NA>
# 25                                                <0.5 yr  <NA>                         8 (80%)                          2 (20%)        <NA>                     <NA>  <NA>         0.05 0.01, 0.29      <0.001                    0.003  <NA>            <NA>       <NA>        <NA>
# 26                         # of participants not included  <NA>                               2                                1        <NA>                     <NA>  <NA>         <NA>       <NA>        <NA>                     <NA>  <NA>            <NA>       <NA>        <NA>
# 27             __Total number of vaccine doses received__   265          3.29 (1.13; 0.00-6.00)           2.61 (1.37; 0.00-5.00)      <0.001                   <0.001   265         0.58 0.45, 0.77      <0.001                   <0.001  <NA>            <NA>       <NA>        <NA>
# 28      __Number of viral-vector vaccine doses received__   265          1.08 (0.92; 0.00-3.00)           0.74 (0.93; 0.00-2.00)       0.003                    0.006   265         0.56 0.39, 0.78      <0.001                    0.002  <NA>            <NA>       <NA>        <NA>
# 29              __Number of mRNA vaccine doses received__   265          1.51 (1.00; 0.00-4.00)           1.39 (0.96; 0.00-4.00)       0.512                    0.544   265         0.97 0.71, 1.33       0.846                    0.931  <NA>            <NA>       <NA>        <NA>
# 30 __Number of inactivated-virus vaccine doses received__   265          0.71 (0.93; 0.00-2.00)           0.46 (0.84; 0.00-2.00)       0.025                    0.038   265         0.67 0.48, 0.93       0.017                    0.029  <NA>            <NA>       <NA>        <NA>
# 31             __Number of other vaccine doses received__   265          0.00 (0.00; 0.00-0.00)           0.01 (0.15; 0.00-2.00)       0.459                    0.521   265        1,950  0.00, Inf       0.995                    0.995  <NA>            <NA>       <NA>        <NA>
# 32                 __Time since most recent vaccination__   263                            <NA>                             <NA>       0.006                    0.010   263         <NA>       <NA>        <NA>                     <NA>  <NA>            <NA>       <NA>        <NA>
# 33                                                  Never  <NA>                         4 (15%)                         22 (85%)        <NA>                     <NA>  <NA>         1.00       <NA>        <NA>                     <NA>  <NA>            <NA>       <NA>        <NA>
# 34                                                  >1 yr  <NA>                        59 (46%)                         70 (54%)        <NA>                     <NA>  <NA>         0.20 0.06, 0.68       0.010                    0.019  <NA>            <NA>       <NA>        <NA>
# 35                                               0.5-1 yr  <NA>                        21 (28%)                         55 (72%)        <NA>                     <NA>  <NA>         0.50 0.14, 1.77       0.279                    0.361  <NA>            <NA>       <NA>        <NA>
# 36                                                <0.5 yr  <NA>                        11 (34%)                         21 (66%)        <NA>                     <NA>  <NA>         0.37 0.09, 1.55       0.175                    0.257  <NA>            <NA>       <NA>        <NA>
# 37                         # of participants not included  <NA>                               0                                2        <NA>                     <NA>  <NA>         <NA>       <NA>        <NA>                     <NA>  <NA>            <NA>       <NA>        <NA>
# 38                                         __Virus wave__   265                            <NA>                             <NA>       0.936                    0.936   265         <NA>       <NA>        <NA>                     <NA>  <NA>            <NA>       <NA>        <NA>
# 39                                  BA.5 and BA.2.75 wave  <NA>                        26 (34%)                         50 (66%)        <NA>                     <NA>  <NA>         1.00       <NA>        <NA>                     <NA>  <NA>            <NA>       <NA>        <NA>
# 40                                       XBB and XBL wave  <NA>                        45 (36%)                         79 (64%)        <NA>                     <NA>  <NA>         0.91 0.40, 2.08       0.824                    0.931  <NA>            <NA>       <NA>        <NA>
# 41                                           BA.2.86 wave  <NA>                        24 (37%)                         41 (63%)        <NA>                     <NA>  <NA>         0.95 0.37, 2.46       0.914                    0.957  <NA>            <NA>       <NA>        <NA>

# Save the table to file
# -------------------------------
gt::gtsave(data = Table_2 %>% as_gt, filename = Table_2_file)
