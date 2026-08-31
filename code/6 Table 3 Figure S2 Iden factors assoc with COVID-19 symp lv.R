# ###############################
# Load global variables
# ###############################
source(file.path("miscellaneous.R"))

# ###############################
# Load required libraries
# ###############################
# library(nnet)
library(ordinal) #clm, clmm
library(gofcat) # brant.test
library(buildmer)
library(performance) #icc
library(emmeans)

library(tidyverse)
library(gtsummary)

library(cowplot) # plot_grid

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
# Identify factors associated with endpoint COVID-19 symptom level
# using multivariable mixed-effects ordinal logistic regression analysis (Pos: minimally/(a)symptomatic < Pos: moderately symptomatic < Pos: highly symptomatic)
# accounting for household-level clustering, with backward BIC-based model selection.
# ###############################
# The initial model included additive fixed effects of all variables + an interaction between first_HH_case_age_cat * age_cat.
# The model doesn't include the total number of vaccine doses received,...
# as it is the sum of the numbers of doses received from individual vaccine platforms, and they were all included in the initial model. 
# The model also doesn't include viral lineage, ...
# due to substantial missing data (32/123 missing observations); otherwise, ...
# the model cannot be fitted reliably, yielding large coeff estimates.
# This was likely due to the small sample size, and the complete separation issue in turn.
# "virus wave" is included however, and it can be thought of as a proxy of "virus major linage".
# -------------------------------
# Prepare data for model fitting (N = 103)
severity_prob_mdl_full_dat <- participant_metadata %>% 
  filter(overall_res_end_with_symp_severity %in% c("Pos: minimally/(a)symptomatic", "Pos: highly symptomatic", "Pos: moderately symptomatic")) %>% droplevels() %>%
  mutate(overall_res_end_with_symp_severity = factor(overall_res_end_with_symp_severity, order = TRUE, levels = c("Pos: minimally/(a)symptomatic", "Pos: moderately symptomatic", "Pos: highly symptomatic"))) %>%
  select(
    overall_res_end_with_symp_severity,
    
    first_HH_case_age_cat, age_cat,
    
    sex,
    BMI,
    pre_existing_med_cond_ht, pre_existing_med_cond_dm, pre_existing_med_cond_al, pre_existing_med_cond_hl,
    
    num_prev_COVID,
    time_to_last_prev_COVID_cat,
    
    # num_vac,
    vac_viral_vector, vac_mRNA, vac_inact_virus, vac_other, 
    time_to_last_vac_cat,
    
    pcr_orf1ab_ct, pcr_n_ct, pcr_e_ct,
    
    # major_lineage,
    wave, # Can be thought of as a proxy of "major_lineage", which, if included, would make the number of data points dropped to 87 from 103!!!
    
    family_id,
  ) %>% na.omit() 

nrow(severity_prob_mdl_full_dat)
# [1] 103

# Model fitting + stepwise BIC-based model selection
severity_prob_mdl_reduced <- buildclmm(
  overall_res_end_with_symp_severity ~
    
    first_HH_case_age_cat * age_cat +
    
    sex +
    BMI +
    
    pre_existing_med_cond_ht + pre_existing_med_cond_dm + pre_existing_med_cond_al + pre_existing_med_cond_hl +
    
    num_prev_COVID +
    time_to_last_prev_COVID_cat +
    
    vac_viral_vector + vac_mRNA + vac_inact_virus + vac_other +
    time_to_last_vac_cat +
    
    pcr_orf1ab_ct + pcr_n_ct + pcr_e_ct +
    
    wave +
    
    (1 | family_id),
  data = severity_prob_mdl_full_dat,
  buildmerControl = buildmerControl(
    include = ~ (1 | family_id),
    direction = "backward",
    crit = "BIC"
  )
)

# Check the proportional odds assumption, evaluated using a fixed-effect-only model
# and all look good :)
brant.test(clm(overall_res_end_with_symp_severity ~ pcr_orf1ab_ct, data = severity_prob_mdl_full_dat))
# 
# Brant Test:
#                  chi-sq   df   pr(>chi)
# Omnibus           0.232    1       0.63
# pcr_orf1ab_ct     0.232    1       0.63
# 
# H0: Proportional odds assumption holds

# Examine the reduced model
summary(severity_prob_mdl_reduced@model)
# Cumulative Link Mixed Model fitted with the Laplace approximation
# 
# formula: overall_res_end_with_symp_severity ~ 1 + pcr_orf1ab_ct + (1 |      family_id)
# data:    severity_prob_mdl_full_dat
# 
#  link  threshold nobs logLik AIC    niter    max.grad cond.H 
#  logit flexible  103  -89.16 186.33 145(313) 6.16e-05 4.6e+04
# 
# Random effects:
#  Groups    Name        Variance Std.Dev.
#  family_id (Intercept) 0.215    0.4637  
# Number of groups:  family_id 56 
# 
# Coefficients:
#               Estimate Std. Error z value Pr(>|z|)   
# pcr_orf1ab_ct -0.11721    0.03685  -3.181  0.00147 **
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Threshold coefficients:
#                                                           Estimate Std. Error z value
# Pos: minimally/(a)symptomatic|Pos: moderately symptomatic   -5.557      1.192  -4.662
# Pos: moderately symptomatic|Pos: highly symptomatic         -3.230      1.004  -3.217

icc(severity_prob_mdl_reduced@model)
# # Intraclass Correlation Coefficient
# 
#     Adjusted ICC: 0.061
#   Unadjusted ICC: 0.053

car::Anova(severity_prob_mdl_reduced@model)
# Analysis of Deviance Table (Type II tests)
# 
# Response: overall_res_end_with_symp_severity
#               Df  Chisq Pr(>Chisq)   
# pcr_orf1ab_ct  1 10.118   0.001468 **
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Repeat the analysis, but with a larger dataset of participants with complete data for
# symptom level, ORF1ab RT-qPCR Ct value, and household identifier (N = 110)
# -------------------------------
# Prepare data for model fitting (N = 110)
severity_prob_mdl_full_dat <- participant_metadata %>% 
  filter(overall_res_end_with_symp_severity %in% c("Pos: highly symptomatic", "Pos: moderately symptomatic", "Pos: minimally/(a)symptomatic")) %>% droplevels() %>%
  mutate(overall_res_end_with_symp_severity = factor(overall_res_end_with_symp_severity, ordered = TRUE, levels = c( "Pos: minimally/(a)symptomatic", "Pos: moderately symptomatic", "Pos: highly symptomatic"))) %>%
  select(overall_res_end_with_symp_severity, pcr_orf1ab_ct, family_id) %>% na.omit()

nrow(severity_prob_mdl_full_dat)
# [1] 110

# Model fitting + stepwise BIC-based model selection
severity_prob_mdl_reduced <- buildclmm(
  overall_res_end_with_symp_severity ~ pcr_orf1ab_ct + (1 | family_id),
  data = severity_prob_mdl_full_dat,
  
  buildmerControl = buildmerControl(
    direction = "backward",
    crit = "BIC",
    include = ~ (1 | family_id)
  )
)

# Check the proportional odds assumption, evaluated using a fixed-effect-only model
# and all look good :)
brant.test(clm(overall_res_end_with_symp_severity ~ pcr_orf1ab_ct, data = severity_prob_mdl_full_dat))
# 
# Brant Test:
#                  chi-sq   df   pr(>chi)
# Omnibus           0.215    1       0.64
# pcr_orf1ab_ct     0.215    1       0.64
# 
# H0: Proportional odds assumption holds

# Examine the reduced model
summary(severity_prob_mdl_reduced@model)
# Cumulative Link Mixed Model fitted with the Laplace approximation
# 
# formula: overall_res_end_with_symp_severity ~ 1 + pcr_orf1ab_ct + (1 |      family_id)
# data:    severity_prob_mdl_full_dat
# 
#  link  threshold nobs logLik AIC    niter    max.grad cond.H 
#  logit flexible  110  -97.58 203.15 133(294) 2.13e-04 3.8e+04
# 
# Random effects:
#  Groups    Name        Variance Std.Dev.
#  family_id (Intercept) 0.449    0.6701  
# Number of groups:  family_id 57 
# 
# Coefficients:
#               Estimate Std. Error z value Pr(>|z|)   
# pcr_orf1ab_ct -0.11062    0.03515  -3.147  0.00165 **
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Threshold coefficients:
#                                                           Estimate Std. Error z value
# Pos: minimally/(a)symptomatic|Pos: moderately symptomatic  -5.3765     1.1606  -4.632
# Pos: moderately symptomatic|Pos: highly symptomatic        -3.0973     0.9809  -3.158

icc(severity_prob_mdl_reduced@model)
# # Intraclass Correlation Coefficient
# 
#     Adjusted ICC: 0.120
#   Unadjusted ICC: 0.105

car::Anova(severity_prob_mdl_reduced@model)
# Analysis of Deviance Table (Type II tests)
# 
# Response: overall_res_end_with_symp_severity
#               Df  Chisq Pr(>Chisq)   
# pcr_orf1ab_ct  1 9.9049   0.001648 **
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Fit the minimal model: Symptom severity ~ orf1ab RT-PCR Ct value with clmm (N = 110)
# The reduced model is already the minimal model, so this is it
# -------------------------------
severity_prob_mdl_minimal <- severity_prob_mdl_reduced@model

confint(severity_prob_mdl_minimal)
#                                                                2.5 %      97.5 %
# Pos: minimally/(a)symptomatic|Pos: moderately symptomatic -7.6512857 -3.10173385
# Pos: moderately symptomatic|Pos: highly symptomatic       -5.0198221 -1.17482552
# pcr_orf1ab_ct                                             -0.1795039 -0.04172856

# ###############################
# Plot overall_res_end_with_symp_severity ~ pcr_orf1ab_ct
# ###############################
# Estimate probabilities at varying pcr_orf1ab_ct values
# -------------------------------
severity_prob_estimates <- emmeans(
  severity_prob_mdl_minimal,
  ~ overall_res_end_with_symp_severity | pcr_orf1ab_ct,
  at = list(pcr_orf1ab_ct = seq(15, 35, 1)),
  mode = "prob",
) %>% data.frame() %>%
  mutate_at(vars(prob, SE, asymp.LCL, asymp.UCL), ~ round((.)*100, 2)) %>%
  mutate(asymp.LCL = ifelse(asymp.LCL<0,0,asymp.LCL)) %>%
  mutate(rate_lab = sprintf("%s\n(%s-%s)", prob, asymp.LCL, asymp.UCL)) %>%
  mutate(
    overall_res_end_with_symp_severity = recode_factor(overall_res_end_with_symp_severity, 
                                                       "1" = "Pos: minimally/(a)symptomatic",
                                                       "2" = "Pos: moderately symptomatic",
                                                       "3" = "Pos: highly symptomatic"
    )
  )

severity_prob_estimates %>% 
  filter(pcr_orf1ab_ct %in% c(20, 25, 30)) %>%
  arrange(overall_res_end_with_symp_severity)
#   overall_res_end_with_symp_severity pcr_orf1ab_ct  prob   SE  df asymp.LCL asymp.UCL             rate_lab
# 1      Pos: minimally/(a)symptomatic            20  4.05 2.15 Inf      0.00      8.26       4.05\n(0-8.26)
# 2      Pos: minimally/(a)symptomatic            25  6.84 2.84 Inf      1.28     12.41   6.84\n(1.28-12.41)
# 3      Pos: minimally/(a)symptomatic            30 11.32 3.94 Inf      3.60     19.05   11.32\n(3.6-19.05)
# 4        Pos: moderately symptomatic            20 25.16 5.76 Inf     13.88     36.44 25.16\n(13.88-36.44)
# 5        Pos: moderately symptomatic            25 34.93 5.23 Inf     24.68     45.19 34.93\n(24.68-45.19)
# 6        Pos: moderately symptomatic            30 44.18 6.23 Inf     31.98     56.39 44.18\n(31.98-56.39)
# 7            Pos: highly symptomatic            20 70.79 7.09 Inf     56.89     84.68 70.79\n(56.89-84.68)
# 8            Pos: highly symptomatic            25 58.22 6.00 Inf     46.46     69.98 58.22\n(46.46-69.98)
# 9            Pos: highly symptomatic            30 44.49 6.34 Inf     32.07     56.92 44.49\n(32.07-56.92)

# Plot the estimated probabilities
# -------------------------------
Figure_S2A <- ggplot(
  severity_prob_estimates,
  mapping = aes(
    x = pcr_orf1ab_ct,
    y = prob,
    ymin = asymp.LCL, ymax = asymp.UCL,
    color = overall_res_end_with_symp_severity,
    fill = overall_res_end_with_symp_severity
  )
) +
  geom_line(linewidth = 0.1) +
  geom_ribbon(alpha = 0.2, color = NA) +
  
  geom_point(data = severity_prob_estimates %>% filter(pcr_orf1ab_ct %in% c(20, 25, 30)), size = 3) +
  geom_pointrange(data = severity_prob_estimates %>% filter(pcr_orf1ab_ct %in% c(20, 25,30)), size = 0.2) +
  geom_text(data = severity_prob_estimates %>% filter(pcr_orf1ab_ct %in% c(20, 25,30)), aes(label = rate_lab), col = "black",
            angle = 90, hjust = 0.5, vjust = -0.5, size = 6*0.35) +
  
  scale_y_continuous(limits = c(0,100)) +
  scale_x_continuous() +
  
  labs(
    x = "ORF1ab RT-qPCR Ct value",
    y = "Predicted probability",
    title = "Predicted probability of symptom levels\nby ORF1ab RT-qPCR Ct value"
  ) +
  scale_colour_manual(
    values = overall_res_end_with_symp_severity_cols,
    breaks = names(overall_res_end_with_symp_severity_cols),
    labels = c(
      "Case without\nsymptom data",
      "Highly\nsymptomatic",
      "Moderately\nsymptomatic",
      "Minimally/\n(a)symptomatic",
      "Uninfected contact",
      "NA"
    ),
    name = "Symptom level"
  ) +
  scale_fill_manual(
    values = overall_res_end_with_symp_severity_cols,
    breaks = names(overall_res_end_with_symp_severity_cols),
    labels = c(
      "Case without\nsymptom data",
      "Highly\nsymptomatic",
      "Moderately\nsymptomatic",
      "Minimally/\n(a)symptomatic",
      "Uninfected contact",
      "NA"
    ),
    name = "Symptom level"
  ) +
  theme_bw() + my_theme +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 6, face = "bold", margin = margin(t = 1, b = 1)), #element_blank(),
  )

# ###############################
# Create Table 3 to summarise the results
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
  
  pcr_orf1ab_ct = "RT-qPCR Ct value of ORF1ab gene",
  pcr_n_ct = "RT-qPCR Ct value of N gene",
  pcr_e_ct = "RT-qPCR Ct value of E gene",
  
  wave = "Virus wave",
  major_lineage = "Virus lineage"
)

# Descriptive statistics
# -------------------------------
Table_3A <- tbl_summary(
  data = participant_metadata %>%
    filter(overall_res_end_with_symp_severity %in% c("Pos: minimally/(a)symptomatic", "Pos: moderately symptomatic", "Pos: highly symptomatic")) %>% droplevels() %>% 
    mutate(overall_res_end_with_symp_severity = factor(overall_res_end_with_symp_severity, ordered = TRUE, levels = c("Pos: minimally/(a)symptomatic", "Pos: moderately symptomatic", "Pos: highly symptomatic"))),
  by = overall_res_end_with_symp_severity,
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
    pcr_orf1ab_ct, pcr_n_ct, pcr_e_ct,
    wave,
    major_lineage
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
  add_n()  %>%
  bold_labels() %>%
  add_p(
    pvalue_fun = label_style_pvalue(digits = 3),
    include = everything()
  ) %>%
  add_q(method = "fdr") %>%
  modify_header(
    all_stat_cols() ~ "**{ifelse(level == 'Pos: minimally/(a)symptomatic', 'Minimally/(a)symptomatic cases', ifelse(level == 'Pos: moderately symptomatic', 'Moderately symptomatic cases', 'Highly symptomatic cases'))}**\n(N = {n})",
    label ~ "**Variable**",
    p.value ~ "**p-value**",
    q.value ~ "**FDR-adjusted p-value**"
  )

# Univariable mixed-effects ordinal logistic regression models
# -------------------------------
Table_3B <- tbl_uvregression(
  data = participant_metadata %>%
    filter(overall_res_end_with_symp_severity %in% c("Pos: minimally/(a)symptomatic", "Pos: highly symptomatic", "Pos: moderately symptomatic")) %>% droplevels() %>%
    mutate(overall_res_end_with_symp_severity = factor(overall_res_end_with_symp_severity, ordered = TRUE, levels = c("Pos: minimally/(a)symptomatic", "Pos: moderately symptomatic", "Pos: highly symptomatic"))),
  method = clmm,
  y = overall_res_end_with_symp_severity,
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
    vac_viral_vector,
    vac_mRNA,
    vac_inact_virus,
    vac_other,
    time_to_last_vac_cat,
    pcr_orf1ab_ct, pcr_n_ct, pcr_e_ct,
    wave,
    major_lineage
  ),
  label = label_list,
  method.args = list(
    link = "logit",
    Hess = TRUE
  ),
  formula = "{y} ~ {x} + (1|family_id)",
  exponentiate = TRUE,
  add_estimate_to_reference_rows = TRUE,
  conf.int = TRUE,
  pvalue_fun = label_style_pvalue(digits = 3),
  hide_n = FALSE,
  show_single_row = c("pre_existing_med_cond_ht", "pre_existing_med_cond_dm", "pre_existing_med_cond_al", "pre_existing_med_cond_hl")
) %>%
  add_q(method = "fdr") %>%
  modify_header(
    estimate ~ "**Crude proportional OR**",
    p.value ~ "**p-value**",
    q.value ~ "**FDR-adjusted p-value**"
  ) %>%
  bold_labels()
# There was a warning running `tbl_regression()` for variable "major_lineage". See message below.
# ! Variance-covariance matrix of the parameters is not defined and Variance-covariance matrix of the parameters is not defined

# Multivariable mixed-effects ordinal logistic regression models
# -------------------------------
Table_3C <- tbl_regression(
  severity_prob_mdl_minimal,
  exponentiate = TRUE,
  conf.int = TRUE,
  label = label_list,
  show_single_row = c(),
  add_estimate_to_reference_rows = TRUE
) %>%
  add_n() %>% 
  modify_header(
    estimate ~ "**Adjusted proportional OR**",
    p.value ~ "**p-value**"
  ) %>%
  bold_labels()

# Combine all tables together
# -------------------------------
Table_3 <- tbl_merge(
  tbls = list(Table_3A, Table_3B, Table_3C),
  tab_spanner = c(
    "**Descriptive statistics**",
    "**Univariable mixed-effects ordinal logistic regression**",
    "**Most parsimonious mixed-effects ordinal logistic regression model**"
  )
)

Table_3 %>% as.data.frame()
#                                              **Variable** **N** **Minimally/(a)symptomatic cases**\n(N = 16) **Moderately symptomatic cases**\n(N = 42) **Highly symptomatic cases**\n(N = 65) **p-value** **FDR-adjusted p-value** **N** **Crude proportional OR** **95% CI** **p-value** **FDR-adjusted p-value** **N** **Adjusted proportional OR** **95% CI** **p-value**
# 1                                      __Household type__   123                                         <NA>                                       <NA>                                   <NA>       0.094                    0.393   123                      <NA>       <NA>        <NA>                     <NA>  <NA>                         <NA>       <NA>        <NA>
# 2                                   Child-first household  <NA>                                     5 (9.8%)                                   23 (45%)                               23 (45%)        <NA>                     <NA>  <NA>                      1.00       <NA>        <NA>                     <NA>  <NA>                         <NA>       <NA>        <NA>
# 3                                   Adult-first household  <NA>                                     11 (15%)                                   19 (26%)                               42 (58%)        <NA>                     <NA>  <NA>                      1.38 0.66, 2.89       0.394                    0.823  <NA>                         <NA>       <NA>        <NA>
# 4                                           __Age group__   123                                         <NA>                                       <NA>                                   <NA>       0.689                    0.906   123                      <NA>       <NA>        <NA>                     <NA>  <NA>                         <NA>       <NA>        <NA>
# 5                                                   Child  <NA>                                     11 (14%)                                   29 (36%)                               40 (50%)        <NA>                     <NA>  <NA>                      1.00       <NA>        <NA>                     <NA>  <NA>                         <NA>       <NA>        <NA>
# 6                                                   Adult  <NA>                                      5 (12%)                                   13 (30%)                               25 (58%)        <NA>                     <NA>  <NA>                      1.39 0.65, 2.97       0.395                    0.823  <NA>                         <NA>       <NA>        <NA>
# 7                                                 __Sex__   123                                         <NA>                                       <NA>                                   <NA>       0.779                    0.906   123                      <NA>       <NA>        <NA>                     <NA>  <NA>                         <NA>       <NA>        <NA>
# 8                                                  Female  <NA>                                     10 (15%)                                   22 (33%)                               35 (52%)        <NA>                     <NA>  <NA>                      1.00       <NA>        <NA>                     <NA>  <NA>                         <NA>       <NA>        <NA>
# 9                                                    Male  <NA>                                      6 (11%)                                   20 (36%)                               30 (54%)        <NA>                     <NA>  <NA>                      1.16 0.56, 2.39       0.694                    0.851  <NA>                         <NA>       <NA>        <NA>
# 10                                                __BMI__   119                    19.54 (4.49; 13.77-29.22)                  20.35 (6.79; 10.18-41.41)              20.75 (6.19; 11.15-44.58)       0.803                    0.906   119                      1.02 0.96, 1.08       0.574                    0.851  <NA>                         <NA>       <NA>        <NA>
# 11                         # of participants not included  <NA>                                            3                                          0                                      1        <NA>                     <NA>  <NA>                      <NA>       <NA>        <NA>                     <NA>  <NA>                         <NA>       <NA>        <NA>
# 12                                       __Hypertension__   123                                      2 (29%)                                    2 (29%)                                3 (43%)       0.399                    0.762   123                      0.51 0.11, 2.37       0.388                    0.823  <NA>                         <NA>       <NA>        <NA>
# 13                                           __Diabetes__   123                                      1 (25%)                                    1 (25%)                                2 (50%)       0.601                    0.906   123                      0.67 0.09, 5.18       0.703                    0.851  <NA>                         <NA>       <NA>        <NA>
# 14                                  __Allergic rhinitis__   123                                      2 (11%)                                    8 (42%)                                9 (47%)       0.778                    0.906   123                      0.81 0.30, 2.17       0.678                    0.851  <NA>                         <NA>       <NA>        <NA>
# 15                                    __Hyperlipidaemia__   123                                      1 (20%)                                    3 (60%)                                1 (20%)       0.315                    0.736   123                      0.33 0.06, 1.74       0.189                    0.602  <NA>                         <NA>       <NA>        <NA>
# 16                      __Number of previous infections__   123                       0.56 (0.51; 0.00-1.00)                     0.62 (0.66; 0.00-2.00)                 0.37 (0.55; 0.00-2.00)       0.084                    0.393   123                      0.56 0.31, 1.02       0.058                    0.323  <NA>                         <NA>       <NA>        <NA>
# 17                   __Time since most recent infection__   123                                         <NA>                                       <NA>                                   <NA>       0.142                    0.495   123                      <NA>       <NA>        <NA>                     <NA>  <NA>                         <NA>       <NA>        <NA>
# 18                                                  Never  <NA>                                      7 (10%)                                   20 (29%)                               43 (61%)        <NA>                     <NA>  <NA>                      1.00       <NA>        <NA>                     <NA>  <NA>                         <NA>       <NA>        <NA>
# 19                                                  >1 yr  <NA>                                      5 (14%)                                   15 (42%)                               16 (44%)        <NA>                     <NA>  <NA>                      0.54 0.23, 1.23       0.143                    0.596  <NA>                         <NA>       <NA>        <NA>
# 20                                               0.5-1 yr  <NA>                                      4 (27%)                                    5 (33%)                                6 (40%)        <NA>                     <NA>  <NA>                      0.33 0.10, 1.07       0.065                    0.323  <NA>                         <NA>       <NA>        <NA>
# 21                                                <0.5 yr  <NA>                                       0 (0%)                                   2 (100%)                                 0 (0%)        <NA>                     <NA>  <NA>                      0.22 0.02, 2.44       0.217                    0.602  <NA>                         <NA>       <NA>        <NA>
# 22             __Total number of vaccine doses received__   123                       2.31 (1.35; 0.00-5.00)                     2.19 (1.35; 0.00-5.00)                 2.34 (1.44; 0.00-5.00)       0.876                    0.906   123                      1.05 0.81, 1.36       0.705                    0.851  <NA>                         <NA>       <NA>        <NA>
# 23      __Number of viral-vector vaccine doses received__   123                       0.50 (0.89; 0.00-2.00)                     0.26 (0.63; 0.00-2.00)                 0.57 (0.88; 0.00-2.00)       0.206                    0.617   123                      1.36 0.85, 2.18       0.197                    0.602  <NA>                         <NA>       <NA>        <NA>
# 24              __Number of mRNA vaccine doses received__   123                       1.56 (0.81; 0.00-3.00)                     1.40 (1.11; 0.00-3.00)                 1.49 (1.05; 0.00-4.00)       0.906                    0.906   123                      1.02 0.72, 1.45       0.895                    0.895  <NA>                         <NA>       <NA>        <NA>
# 25 __Number of inactivated-virus vaccine doses received__   123                       0.25 (0.68; 0.00-2.00)                     0.48 (0.83; 0.00-2.00)                 0.28 (0.70; 0.00-2.00)       0.262                    0.687   123                      0.86 0.54, 1.37       0.526                    0.851  <NA>                         <NA>       <NA>        <NA>
# 26             __Number of other vaccine doses received__   123                       0.00 (0.00; 0.00-0.00)                     0.05 (0.31; 0.00-2.00)                 0.00 (0.00; 0.00-0.00)       0.381                    0.762   123                      0.55 0.10, 2.92       0.482                    0.851  <NA>                         <NA>       <NA>        <NA>
# 27                 __Time since most recent vaccination__   121                                         <NA>                                       <NA>                                   <NA>       0.851                    0.906   121                      <NA>       <NA>        <NA>                     <NA>  <NA>                         <NA>       <NA>        <NA>
# 28                                                  Never  <NA>                                     2 (9.1%)                                    8 (36%)                               12 (55%)        <NA>                     <NA>  <NA>                      1.00       <NA>        <NA>                     <NA>  <NA>                         <NA>       <NA>        <NA>
# 29                                                  >1 yr  <NA>                                      9 (19%)                                   16 (33%)                               23 (48%)        <NA>                     <NA>  <NA>                      0.68 0.25, 1.86       0.452                    0.851  <NA>                         <NA>       <NA>        <NA>
# 30                                               0.5-1 yr  <NA>                                     3 (8.3%)                                   12 (33%)                               21 (58%)        <NA>                     <NA>  <NA>                      1.15 0.39, 3.35       0.800                    0.870  <NA>                         <NA>       <NA>        <NA>
# 31                                                <0.5 yr  <NA>                                     1 (6.7%)                                    5 (33%)                                9 (60%)        <NA>                     <NA>  <NA>                      1.29 0.33, 4.94       0.714                    0.851  <NA>                         <NA>       <NA>        <NA>
# 32                         # of participants not included  <NA>                                            1                                          1                                      0        <NA>                     <NA>  <NA>                      <NA>       <NA>        <NA>                     <NA>  <NA>                         <NA>       <NA>        <NA>
# 33                    __RT-qPCR Ct value of ORF1ab gene__   110                    31.03 (5.00; 21.21-38.85)                  27.62 (6.42; 17.21-39.85)              24.66 (6.26; 15.31-40.04)       0.002                    0.011   110                      0.90 0.84, 0.96       0.002                    0.019   110                         0.90 0.84, 0.96       0.002
# 34                         # of participants not included  <NA>                                            4                                          3                                      6        <NA>                     <NA>  <NA>                      <NA>       <NA>        <NA>                     <NA>  <NA>                         <NA>       <NA>        <NA>
# 35                         __RT-qPCR Ct value of N gene__   110                    28.71 (4.90; 19.94-36.46)                  25.29 (6.21; 14.33-37.87)              22.46 (6.55; 13.52-38.77)       0.001                    0.011   110                      0.90 0.84, 0.96       0.002                    0.019  <NA>                         <NA>       <NA>        <NA>
# 36                         # of participants not included  <NA>                                            4                                          3                                      6        <NA>                     <NA>  <NA>                      <NA>       <NA>        <NA>                     <NA>  <NA>                         <NA>       <NA>        <NA>
# 37                         __RT-qPCR Ct value of E gene__   107                    29.76 (4.77; 20.67-38.11)                  26.05 (6.34; 16.32-38.45)              22.66 (5.91; 14.28-39.21)      <0.001                    0.004   107                      0.88 0.81, 0.95      <0.001                    0.016  <NA>                         <NA>       <NA>        <NA>
# 38                         # of participants not included  <NA>                                            4                                          4                                      8        <NA>                     <NA>  <NA>                      <NA>       <NA>        <NA>                     <NA>  <NA>                         <NA>       <NA>        <NA>
# 39                                         __Virus wave__   123                                         <NA>                                       <NA>                                   <NA>       0.820                    0.906   123                      <NA>       <NA>        <NA>                     <NA>  <NA>                         <NA>       <NA>        <NA>
# 40                                  BA.5 and BA.2.75 wave  <NA>                                      5 (13%)                                   13 (33%)                               22 (55%)        <NA>                     <NA>  <NA>                      1.00       <NA>        <NA>                     <NA>  <NA>                         <NA>       <NA>        <NA>
# 41                                       XBB and XBL wave  <NA>                                      8 (16%)                                   15 (31%)                               26 (53%)        <NA>                     <NA>  <NA>                      0.89 0.36, 2.16       0.791                    0.870  <NA>                         <NA>       <NA>        <NA>
# 42                                           BA.2.86 wave  <NA>                                     3 (8.8%)                                   14 (41%)                               17 (50%)        <NA>                     <NA>  <NA>                      0.93 0.36, 2.44       0.886                    0.895  <NA>                         <NA>       <NA>        <NA>
# 43                                      __Virus lineage__    91                                         <NA>                                       <NA>                                   <NA>       0.861                    0.906    91                      <NA>       <NA>        <NA>                     <NA>  <NA>                         <NA>       <NA>        <NA>
# 44                                                   BA.5  <NA>                                       0 (0%)                                    1 (25%)                                3 (75%)        <NA>                     <NA>  <NA>                      1.00       <NA>        <NA>                     <NA>  <NA>                         <NA>       <NA>        <NA>
# 45                                                BA.2.75  <NA>                                     2 (8.7%)                                    6 (26%)                               15 (65%)        <NA>                     <NA>  <NA>                      0.57       <NA>        <NA>                     <NA>  <NA>                         <NA>       <NA>        <NA>
# 46                                             XBB or XBL  <NA>                                      4 (11%)                                   10 (28%)                               22 (61%)        <NA>                     <NA>  <NA>                      0.47       <NA>        <NA>                     <NA>  <NA>                         <NA>       <NA>        <NA>
# 47                                                BA.2.86  <NA>                                     2 (7.1%)                                   12 (43%)                               14 (50%)        <NA>                     <NA>  <NA>                      0.34       <NA>        <NA>                     <NA>  <NA>                         <NA>       <NA>        <NA>
# 48                         # of participants not included  <NA>                                            8                                         13                                     11        <NA>                     <NA>  <NA>                      <NA>       <NA>        <NA>                     <NA>  <NA>                         <NA>       <NA>        <NA>

# Save the table to file
# -------------------------------
gt::gtsave(data = Table_3 %>% as_gt, filename = Table_3_file)

# ###############################
# Further characterise the relationship between symptom level and ORF1ab RT-qPCR Ct value
# with reciprocal mixed-effects linear regression analysis of factors associated with ORF1ab Ct value 
# ###############################
# Prepare data for model fitting (N = 103)
# Basically, the same one used in the full full analysis of symptom level
# -------------------------------
pcr_orf1ab_ct_mdl_full_dat <- participant_metadata %>% 
  filter(overall_res_end_with_symp_severity %in% c("Pos: minimally/(a)symptomatic", "Pos: moderately symptomatic", "Pos: highly symptomatic")) %>% droplevels() %>%
  mutate(overall_res_end_with_symp_severity = factor(overall_res_end_with_symp_severity, order = TRUE, levels = c("Pos: minimally/(a)symptomatic", "Pos: moderately symptomatic", "Pos: highly symptomatic"))) %>%
  select(
    overall_res_end_with_symp_severity,
    
    first_HH_case_age_cat, age_cat,
    
    sex,
    BMI,
    pre_existing_med_cond_ht, pre_existing_med_cond_dm, pre_existing_med_cond_al, pre_existing_med_cond_hl,
    
    num_prev_COVID,
    time_to_last_prev_COVID_cat,
    
    vac_viral_vector, vac_mRNA, vac_inact_virus, vac_other, 
    time_to_last_vac_cat,
    
    pcr_orf1ab_ct, pcr_n_ct, pcr_e_ct,
    
    wave,
    
    family_id,
  ) %>% na.omit() 

nrow(pcr_orf1ab_ct_mdl_full_dat)
#[1] 103

# Model fitting + stepwise BIC-based model selection
pcr_orf1ab_ct_mdl_reduced <- buildmer(
  pcr_orf1ab_ct ~ 
    overall_res_end_with_symp_severity + 
    
    first_HH_case_age_cat * age_cat +
    
    sex +
    BMI +
    pre_existing_med_cond_ht + pre_existing_med_cond_dm + pre_existing_med_cond_al + pre_existing_med_cond_hl +
    
    num_prev_COVID +
    time_to_last_prev_COVID_cat +
    
    vac_viral_vector + vac_mRNA + vac_inact_virus + vac_other + 
    time_to_last_vac_cat +
    
    wave +
    (1|family_id), 
  data = pcr_orf1ab_ct_mdl_full_dat,
  
  buildmerControl = buildmerControl(
    include=~ (1 | family_id),
    direction = "backward",
    crit = "BIC"
  )
)

# Examine the reduced model
summary(pcr_orf1ab_ct_mdl_reduced@model)
# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: pcr_orf1ab_ct ~ 1 + overall_res_end_with_symp_severity + pre_existing_med_cond_al +      (1 | family_id)
#    Data: pcr_orf1ab_ct_mdl_full_dat
# 
# REML criterion at convergence: 641.1
# 
# Scaled residuals: 
#     Min      1Q  Median      3Q     Max 
# -1.6922 -0.6695 -0.2325  0.6133  2.7691 
# 
# Random effects:
#  Groups    Name        Variance Std.Dev.
#  family_id (Intercept)  2.033   1.426   
#  Residual              31.195   5.585   
# Number of obs: 103, groups:  family_id, 56
# 
# Fixed effects:
#                                      Estimate Std. Error      df t value Pr(>|t|)    
# (Intercept)                           26.7876     0.8053 49.3331  33.263  < 2e-16 ***
# overall_res_end_with_symp_severity.L  -4.5676     1.4025 98.8032  -3.257  0.00154 ** 
# overall_res_end_with_symp_severity.Q   0.3237     1.1174 98.9366   0.290  0.77270    
# pre_existing_med_cond_alYes            3.8472     1.5300 98.8206   2.515  0.01353 *  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#             (Intr) o_____.L o_____.Q
# ovrl_____.L -0.597                  
# ovrl_____.Q  0.323 -0.504           
# pr_xstn___Y -0.350  0.056    0.017  

icc(pcr_orf1ab_ct_mdl_reduced@model)
# # Intraclass Correlation Coefficient
# 
#     Adjusted ICC: 0.061
#   Unadjusted ICC: 0.051

car::Anova(pcr_orf1ab_ct_mdl_reduced@model)
# Analysis of Deviance Table (Type II Wald chisquare tests)
# 
# Response: pcr_orf1ab_ct
#                                      Chisq Df Pr(>Chisq)   
# overall_res_end_with_symp_severity 13.0505  2   0.001466 **
# pre_existing_med_cond_al            6.3231  1   0.011918 * 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Repeat the analysis, but with a larger dataset of participants with complete data for
# ORF1ab RT-qPCR Ct value, symptom level, allergic rhinitis status, and household identifier (N = 110)
# to find that, indeed, both symptom level, and allergic rhinitis status are sig assoc with orf1ab RT-PCR Ct value
# -------------------------------
# Prepare data for model fitting (N = 110)
pcr_orf1ab_ct_mdl_reduced_dat <- participant_metadata %>% 
  filter(overall_res_end_with_symp_severity %in% c("Pos: minimally/(a)symptomatic", "Pos: highly symptomatic", "Pos: moderately symptomatic")) %>% droplevels() %>%
  mutate(overall_res_end_with_symp_severity = factor(overall_res_end_with_symp_severity, ordered = TRUE, levels = c( "Pos: minimally/(a)symptomatic", "Pos: moderately symptomatic", "Pos: highly symptomatic"))) %>%
  select(pcr_orf1ab_ct, overall_res_end_with_symp_severity, pre_existing_med_cond_al, family_id) %>% na.omit()

nrow(pcr_orf1ab_ct_mdl_reduced_dat)
# [1] 110

# Model fitting + stepwise BIC-based model selection
# NOTE:"Singular fit" occurred, so singular.ok = T is needed to retain (1|family_id) !!!
# Otherwise (1|family_id) gets dropped...
pcr_orf1ab_ct_mdl_reduced <- buildmer(
  pcr_orf1ab_ct ~ overall_res_end_with_symp_severity + pre_existing_med_cond_al + (1|family_id), 
  # pcr_orf1ab_ct ~ overall_res_end_with_symp_severity * pre_existing_med_cond_al + (1|family_id), # give the same result as above
  data = pcr_orf1ab_ct_mdl_reduced_dat,
  buildmerControl = buildmerControl(
    include = ~ (1 | family_id),
    direction = "backward",
    crit = "BIC",
    singular.ok = T
  )
)

# Examine the reduced model
summary(pcr_orf1ab_ct_mdl_reduced@model)
# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: pcr_orf1ab_ct ~ 1 + overall_res_end_with_symp_severity + pre_existing_med_cond_al +      (1 | family_id)
#    Data: pcr_orf1ab_ct_mdl_reduced_dat
# 
# REML criterion at convergence: 696
# 
# Scaled residuals: 
#     Min      1Q  Median      3Q     Max 
# -1.7641 -0.7278 -0.2447  0.5518  2.6267 
# 
# Random effects:
#  Groups    Name        Variance Std.Dev.
#  family_id (Intercept)  0.09095 0.3016  
#  Residual              36.35331 6.0294  
# Number of obs: 110, groups:  family_id, 57
# 
# Fixed effects:
#                                      Estimate Std. Error       df t value Pr(>|t|)    
# (Intercept)                           27.0914     0.7608  48.7813  35.609  < 2e-16 ***
# overall_res_end_with_symp_severity.L  -4.4681     1.3522 102.7852  -3.304  0.00131 ** 
# overall_res_end_with_symp_severity.Q   0.2460     1.1104 105.0137   0.222  0.82511    
# pre_existing_med_cond_alYes            4.1005     1.5568 105.9816   2.634  0.00970 ** 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#             (Intr) o_____.L o_____.Q
# ovrl_____.L -0.559                  
# ovrl_____.Q  0.281 -0.465           
# pr_xstn___Y -0.341  0.012    0.023  

confint(pcr_orf1ab_ct_mdl_reduced@model)
# Computing profile confidence intervals ...
#                                          2.5 %    97.5 %
# .sig01                                0.000000  3.038534
# .sigma                                5.222658  6.805847
# (Intercept)                          25.615141 28.563435
# overall_res_end_with_symp_severity.L -7.090090 -1.842396
# overall_res_end_with_symp_severity.Q -1.909950  2.399908
# pre_existing_med_cond_alYes           1.074010  7.117060
# Warning messages:
# 1: In nextpar(mat, cc, i, delta, lowcut, upcut) :
#   Last two rows have identical or NA .zeta values: using minstep
# 2: In FUN(X[[i]], ...) : non-monotonic profile for .sig01
# 3: In nextpar(mat, cc, i, delta, lowcut, upcut) :
#   unexpected decrease in profile: using minstep
# 4: In FUN(X[[i]], ...) : non-monotonic profile for .sigma
# 5: In confint.thpr(pp, level = level, zeta = zeta) :
#   bad spline fit for .sig01: falling back to linear interpolation
# 6: In regularize.values(x, y, ties, missing(ties), na.rm = na.rm) :
#   collapsing to unique 'x' values
# 7: In confint.thpr(pp, level = level, zeta = zeta) :
#   bad spline fit for .sigma: falling back to linear interpolation

icc(pcr_orf1ab_ct_mdl_reduced@model)
# # Intraclass Correlation Coefficient
# 
#     Adjusted ICC: 0.002
#   Unadjusted ICC: 0.002

car::Anova(pcr_orf1ab_ct_mdl_reduced@model)
# Analysis of Deviance Table (Type II Wald chisquare tests)
# 
# Response: pcr_orf1ab_ct
#                                      Chisq Df Pr(>Chisq)   
# overall_res_end_with_symp_severity 13.1256  2   0.001412 **
# pre_existing_med_cond_al            6.9377  1   0.008440 **
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Fit the minimal model: overall_res_end_with_symp_severity + pre_existing_med_cond_al + (1 | family_id)
# -------------------------------
# However, since all variables are kept in the model selection,
# the final parsimonious model obtained from the model selection is thereby 
# the minimal model fitted to the largest available data also already
pcr_orf1ab_ct_mdl_minimal <- pcr_orf1ab_ct_mdl_reduced@model

# ###############################
# Estimate marginalised Ct values for various participant groups based on the minimal model
# ###############################
# For overall_res_end_with_symp_severity * pre_existing_med_cond_al
# -------------------------------
emmeans(
  pcr_orf1ab_ct_mdl_minimal, 
  ~ overall_res_end_with_symp_severity * pre_existing_med_cond_al, 
  weights = "proportional"
) %>% as.data.frame() %>% 
  mutate(label = sprintf("%.2f\n(%.2f-%.2f)", emmean, lower.CL, upper.CL))
#   overall_res_end_with_symp_severity pre_existing_med_cond_al   emmean        SE        df lower.CL upper.CL                label
# 1      Pos: minimally/(a)symptomatic                       No 30.35122 1.8154484  87.14328 26.74290 33.95953 30.35\n(26.74-33.96)
# 2        Pos: moderately symptomatic                       No 26.89054 1.0248229  89.51884 24.85440 28.92668 26.89\n(24.85-28.93)
# 3            Pos: highly symptomatic                       No 24.03244 0.8381784  71.07770 22.36119 25.70368 24.03\n(22.36-25.70)
# 4      Pos: minimally/(a)symptomatic                      Yes 34.45169 2.2196904 105.41971 30.05066 38.85272 34.45\n(30.05-38.85)
# 5        Pos: moderately symptomatic                      Yes 30.99102 1.6262574 105.95243 27.76679 34.21525 30.99\n(27.77-34.22)
# 6            Pos: highly symptomatic                      Yes 28.13291 1.5609243 102.17609 25.03689 31.22893 28.13\n(25.04-31.23)

# For symp categories
# -------------------------------
emmeans(
  pcr_orf1ab_ct_mdl_minimal, 
  ~ overall_res_end_with_symp_severity, 
  weights = "proportional"
) %>% summary(infer = c(TRUE, TRUE)) %>% as.data.frame() 
#  overall_res_end_with_symp_severity   emmean        SE    df lower.CL upper.CL t.ratio p.value
#  Pos: minimally/(a)symptomatic      31.02220 1.7943856 90.42 27.45757 34.58684  17.288 <0.0001
#  Pos: moderately symptomatic        27.56152 0.9839816 91.79 25.60719 29.51586  28.010 <0.0001
#  Pos: highly symptomatic            24.70342 0.8020765 67.66 23.10276 26.30409  30.799 <0.0001
# 
# Results are averaged over the levels of: pre_existing_med_cond_al 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 

# For with and without allergic rhinitis
# -------------------------------
emmeans(
  pcr_orf1ab_ct_mdl_minimal, 
  ~ pre_existing_med_cond_al, 
  weights = "proportional"
) %>% summary(infer = c(TRUE, TRUE)) %>% as.data.frame() 
#  pre_existing_med_cond_al   emmean        SE     df lower.CL upper.CL t.ratio p.value
#  No                       25.73508 0.6449335  42.40 24.43391 27.03625  39.903 <0.0001
#  Yes                      29.83556 1.4454113 104.73 26.96949 32.70163  20.642 <0.0001
# 
# Results are averaged over the levels of: overall_res_end_with_symp_severity 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 

# ###############################
# Estimate adjusted mean difference in Ct value between various participant groups based on the min model
# ###############################
# Between symp categories
# -------------------------------
emmeans(
  pcr_orf1ab_ct_mdl_minimal, 
  pairwise ~ overall_res_end_with_symp_severity, 
  weights = "proportional"
) %>% summary(infer = c(TRUE, TRUE)) %>% 
  pluck("contrasts") %>% as.data.frame() 
#  contrast                                                      estimate       SE     df   lower.CL  upper.CL t.ratio p.value
#  (Pos: minimally/(a)symptomatic) - Pos: moderately symptomatic 3.460678 2.042745 103.64 -1.3966722  8.318028   1.694  0.2123
#  (Pos: minimally/(a)symptomatic) - Pos: highly symptomatic     6.318781 1.959992 102.70  1.6575892 10.979972   3.224  0.0048
#  Pos: moderately symptomatic - Pos: highly symptomatic         2.858103 1.269795 105.64 -0.1604537  5.876659   2.251  0.0675
# 
# Results are averaged over the levels of: pre_existing_med_cond_al 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 
# Conf-level adjustment: tukey method for comparing a family of 3 estimates 
# P value adjustment: tukey method for comparing a family of 3 estimates 

# Between those with VS without allergic rhinitis
# -------------------------------
emmeans(
  pcr_orf1ab_ct_mdl_minimal, 
  pairwise ~ pre_existing_med_cond_al, 
  weights = "proportional"
) %>% summary(infer = c(TRUE, TRUE)) %>% 
  pluck("contrasts") %>% as.data.frame() 
#  contrast  estimate       SE     df  lower.CL   upper.CL t.ratio p.value
#  No - Yes -4.100477 1.582995 105.98 -7.238925 -0.9620303  -2.590  0.0109
# 
# Results are averaged over the levels of: overall_res_end_with_symp_severity 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 

# ###############################
# Plot the distribution of PCR orf1ab Ct value by symptom level and allergic rhinitis conditions
# ###############################
obs_pcr_orf1ab_ct <- pcr_orf1ab_ct_mdl_reduced_dat %>%
  group_by(overall_res_end_with_symp_severity, pre_existing_med_cond_al) %>%
  summarise(
    med_pcr_orf1ab_ct = median(pcr_orf1ab_ct),
    SD_pcr_orf1ab_ct = sd(pcr_orf1ab_ct),
    .groups = "drop"
  ) %>%
  mutate(label = sprintf("%.2f\n(SD = %.2f)", med_pcr_orf1ab_ct, SD_pcr_orf1ab_ct))

obs_pcr_orf1ab_ct
#   overall_res_end_with_symp_severity pre_existing_med_cond_al med_pcr_orf1ab_ct SD_pcr_orf1ab_ct label               
#   <ord>                              <fct>                                <dbl>            <dbl> <chr>               
# 1 Pos: minimally/(a)symptomatic      No                                    30.5             5.34 "30.47\n(SD = 5.34)"
# 2 Pos: minimally/(a)symptomatic      Yes                                   29.4             3.48 "29.44\n(SD = 3.48)"
# 3 Pos: moderately symptomatic        No                                    24.8             6.37 "24.80\n(SD = 6.37)"
# 4 Pos: moderately symptomatic        Yes                                   30.1             5.61 "30.13\n(SD = 5.61)"
# 5 Pos: highly symptomatic            No                                    22.5             5.62 "22.49\n(SD = 5.62)"
# 6 Pos: highly symptomatic            Yes                                   29.1             8.14 "29.13\n(SD = 8.14)"

Figure_S2B <- ggplot(
  data = pcr_orf1ab_ct_mdl_reduced_dat, 
  mapping = aes(
    y = pcr_orf1ab_ct, 
    x = overall_res_end_with_symp_severity, 
    fill = pre_existing_med_cond_al
  ) 
) +
  geom_boxplot(outliers  = F) + 
  geom_point(position = position_jitterdodge(), shape = 21, col = "black") +
  geom_text(
    data = obs_pcr_orf1ab_ct %>% filter(pre_existing_med_cond_al == "No"), 
    mapping = aes(
      x = overall_res_end_with_symp_severity, 
      y = med_pcr_orf1ab_ct,
      label = label
    ),
    angle = 90, vjust = -1.75, hjust = 0.5,
    size = 6/.pt
  ) + 
  geom_text(
    data = obs_pcr_orf1ab_ct %>% filter(pre_existing_med_cond_al == "Yes"), 
    mapping = aes(
      x = overall_res_end_with_symp_severity, 
      y = med_pcr_orf1ab_ct,
      label = label
    ),
    angle = -90, vjust = -1.75, hjust = 0.5,
    size = 6/.pt
  ) + 
  scale_x_discrete(
    labels = c(
      "Pos: minimally/(a)symptomatic" = "Minimally/\n(a)symptomatic",
      "Pos: moderately symptomatic" = "Moderately\nsymptomatic",
      "Pos: highly symptomatic" = "Highly\nsymptomatic"
    )
  ) +
  labs(
    x = "Symptom level",
    y = "ORF1ab Ct value",
    title = "ORF1ab RT-qPCR Ct value by\nsymptom level and allergic rhinitis condition",
    fill = "Allergic rhinitis condition"
  ) +
  theme_bw() + my_theme +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 6, face = "bold", margin = margin(t = 1, b = 1)), 
  )

# ###############################
# Combine Figure S2A and S2B, and save it to file
# ###############################
# Combine Figure S2A and S2B
# -------------------------------
Figure_S2 <- plot_grid(
  Figure_S2A,
  Figure_S2B,
  labels = c("A)", "B)"), label_size = 8,
  label_x = 0, label_y = 1, hjust = 0, vjust = 1,
  nrow = 1, align = "h", axis = "tb", rel_widths = c(1, 1)
)

# Save it to file
# -------------------------------
ggsave(filename = Figure_S2_file_png, device = "png",
       plot = Figure_S2,
       width = 16, height = 8, units = "cm", dpi = 300, bg = "white")

ggsave(filename = Figure_S2_file_svg, device = "svg",
       plot = Figure_S2,
       width = 16, height = 8, units = "cm", dpi = 300, bg = "white")
