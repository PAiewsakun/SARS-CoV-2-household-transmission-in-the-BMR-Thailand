# ###############################
# Load global variables
# ###############################
source(file.path("miscellaneous.R"))

# ###############################
# Load required libraries
# ###############################
library(tidyverse)
library(emmeans)
library(gtsummary)

library(ggh4x)

# ###############################
# Load and prepare metadata
# ###############################
# Participant metadata
# -------------------------------
participant_metadata <- read.delim(participant_metadata_filepath, header = TRUE, sep = "\t")

# Add first_HH_case_age_cat to participant_metadata
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

# Add cluster_id to participant_metadata
# -------------------------------
cluster_metadata_hard_pruned_tree <- read.delim(Table_S9_file, header = TRUE, sep = "\t")
seq_members_list <- strsplit(cluster_metadata_hard_pruned_tree$seq_members, ";\\s*")
cluster_assignments <- setNames(
  rep(cluster_metadata_hard_pruned_tree$cluster_id, lengths(seq_members_list)),
  unlist(seq_members_list)
)

participant_metadata <- participant_metadata %>%
  mutate(
    cluster_id = cluster_assignments[seq_id], .after = subject_id,
    major_lineage = ifelse(major_lineage %in% c("XBB", "XBL"), "XBB or XBL", as.character(major_lineage))
  )

# Select only households with complete sequencing data
# --------------------
HH_with_complete_seq_dat <- participant_metadata %>%
  select(family_id, overall_res_end, seq_id) %>%
  filter(overall_res_end == "Pos") %>%
  group_by(family_id) %>%
  summarise(
    N_pos = n(),
    seq_avail = sum(!is.na(seq_id)) / N_pos,
    .groups = "drop"
  ) %>%
  filter(seq_avail == 1) %>%
  pull(family_id)

participant_metadata_from_HH_with_complete_seq_dat <- participant_metadata %>% 
  filter(family_id %in% HH_with_complete_seq_dat) %>%
  select(
    subject_id, family_id, age_cat, first_HH_case_age_cat, overall_res_end, group_wrt_index_case, time_to_overall_min_pos_date, major_lineage, cluster_id,
  ) %>%
  mutate_at(vars(subject_id, family_id, age_cat, first_HH_case_age_cat, overall_res_end, group_wrt_index_case, major_lineage, cluster_id), factor) %>%
  mutate(time_to_overall_min_pos_date = as.integer(time_to_overall_min_pos_date)) %>%
  arrange(family_id, overall_res_end, time_to_overall_min_pos_date)

# Add major_lineage_family_lv to participant_metadata_from_HH_with_complete_seq_dat
# -------------------------------
participant_metadata_from_HH_with_complete_seq_dat <- participant_metadata_from_HH_with_complete_seq_dat %>%
  group_by(family_id) %>%
  mutate(
    major_lineage_family_lv = unique(na.omit(major_lineage))
  ) %>% 
  ungroup()

# ###############################
# Compute some summary stats
# ###############################
# Household type count
# -------------------------------
participant_metadata_from_HH_with_complete_seq_dat %>% 
  select(family_id, first_HH_case_age_cat) %>% unique() %>% 
  pull(first_HH_case_age_cat) %>% 
  table(useNA = "ifany") %>% addmargins()
# .
# Child-first household Adult-first household                   Sum 
#                    18                    11                    29 

# Number of clusters with complete seq data
# -------------------------------
participant_metadata_from_HH_with_complete_seq_dat %>% pull(cluster_id) %>% na.omit() %>% unique() %>% length
# [1] 33

# Number of participants by overall_res_end, group_wrt_index_case, age_cat
# -------------------------------
participant_metadata_from_HH_with_complete_seq_dat %>% 
  select(overall_res_end, group_wrt_index_case, age_cat) %>% table(useNA = "ifany") %>% addmargins()
# , , age_cat = Adult
# 
#                group_wrt_index_case
# overall_res_end Co-index case Index case Post-index case / incident case Pre-index case / additional case identified at enrolment Uninfected contact Sum
#             Neg             0          0                               0                                                        0                 51  51
#             Pos             3          0                              11                                                       17                  0  31
#             Sum             3          0                              11                                                       17                 51  82
# 
# , , age_cat = Child
# 
#                group_wrt_index_case
# overall_res_end Co-index case Index case Post-index case / incident case Pre-index case / additional case identified at enrolment Uninfected contact Sum
#             Neg             0          0                               0                                                        0                  8   8
#             Pos             0         30                               3                                                        1                  0  34
#             Sum             0         30                               3                                                        1                  8  42
# 
# , , age_cat = Sum
# 
#                group_wrt_index_case
# overall_res_end Co-index case Index case Post-index case / incident case Pre-index case / additional case identified at enrolment Uninfected contact Sum
#             Neg             0          0                               0                                                        0                 59  59
#             Pos             3         30                              14                                                       18                  0  65
#             Sum             3         30                              14                                                       18                 59 124

# Number of families with more than 1 cluster
# -------------------------------
participant_metadata_from_HH_with_complete_seq_dat %>% 
  select(family_id, first_HH_case_age_cat, cluster_id) %>% na.omit() %>% unique() %>%
  group_by(family_id, first_HH_case_age_cat) %>%
  summarise(n_clusters = n_distinct(cluster_id),
            clusters = paste(unique(cluster_id), collapse = ", ")) %>%
  filter(n_clusters > 1)
#   family_id first_HH_case_age_cat n_clusters clusters
#   <fct>     <fct>                      <int> <chr>   
# 1 F023      Adult-first household          2 25, 26  
# 2 F029      Adult-first household          2 15, 14  
# 3 F031      Adult-first household          2 10, 9   
# 4 F047      Adult-first household          2 65, 61  
# 5 F054      Adult-first household          2 51, 50

# Number of clusters with more than 1 family
# -------------------------------
participant_metadata_from_HH_with_complete_seq_dat %>% 
  select(family_id, cluster_id) %>% na.omit() %>% unique() %>%
  group_by(cluster_id) %>%
  summarise(n_families = n_distinct(family_id),
            families = paste(unique(family_id), collapse = ", ")) %>%
  filter(n_families > 1)
#   cluster_id n_families families  
#   <fct>           <int> <chr>     
# 1 20                  2 F044, F045

# ###############################
# Run simulations to compute:
#  (hidden) community acquisition rates, and
#  corrected attack rates, 
# for various participant groups
# ###############################
# Setup
# -------------------------------
set.seed(100)

n_sims <- 1000

sim_hidden_comm_acq_rate_aHH_log_odds <- list()
sim_hidden_comm_acq_rate_aHH_estimates <- list()
sim_hidden_comm_acq_rate_aHH_pvals <- list()

sim_atk_rate_log_odds <- list()
sim_atk_rate_estimates <- list()
sim_atk_rate_pvals <- list()

# Run simulations
# -------------------------------
for(sim in 1:n_sims) {
  print(sim)
  
  # ###############################
  # Assign primary cases per family_id × cluster_id / per family_id
  # ###############################
  # Iden primary per family_id × cluster_id; 
  # If there are multiple cases with time_to_overall_min_pos_date == min(time_to_overall_min_pos_date), randomly pick one among them
  # ------------------------------------------------------
  primaries_cluster_x_fam <- participant_metadata_from_HH_with_complete_seq_dat %>%
    filter(overall_res_end == "Pos") %>% 
    group_by(family_id, cluster_id) %>%
    filter(time_to_overall_min_pos_date == min(time_to_overall_min_pos_date)) %>%
    slice_sample(n = 1) %>%
    ungroup() %>%
    select(subject_id)
  
  # Tag primaries_cluster_x_fam in participant_metadata_from_HH_with_complete_seq_dat
  # ------------------------------------------------------
  participant_metadata_from_HH_with_complete_seq_dat <- participant_metadata_from_HH_with_complete_seq_dat %>% 
    mutate(primary_cluster_x_fam = ifelse(subject_id %in% primaries_cluster_x_fam$subject_id, 1, 0))
  
  # Randomly assign one primary per family_id from primaries_cluster_x_fam for consistency
  # ------------------------------------------------------
  primaries_fam <- participant_metadata_from_HH_with_complete_seq_dat %>%
    filter(primary_cluster_x_fam == 1) %>% 
    group_by(family_id) %>%
    filter(time_to_overall_min_pos_date == min(time_to_overall_min_pos_date)) %>%
    slice_sample(n = 1) %>%
    ungroup() %>%
    select(subject_id)
  
  # Tag primaries_fam in participant_metadata_from_HH_with_complete_seq_dat
  # ------------------------------------------------------
  participant_metadata_from_HH_with_complete_seq_dat <- participant_metadata_from_HH_with_complete_seq_dat %>% 
    mutate(primary_fam = ifelse(subject_id %in% primaries_fam$subject_id, 1, 0)) 
  
  # ###############################
  # Compute hidden community acquisition rates using glm, focusing on adult-first households (aHH)
  # ###############################
  # Prepare data
  # ------------------------------------------------------
  hidden_comm_acq_aHH_dat <- participant_metadata_from_HH_with_complete_seq_dat %>% 
    filter(first_HH_case_age_cat == "Adult-first household") %>%
    filter(overall_res_end == "Pos") %>% 
    group_by(family_id, age_cat, major_lineage_family_lv) %>%
    summarise(
      non_primary_cluster_x_fam = sum(primary_cluster_x_fam == 0),
      primary_cluster_x_fam = sum(primary_cluster_x_fam == 1),
      
      non_primary_fam = sum(primary_fam == 0),
      primary_fam = sum(primary_fam == 1),
      
      hidden_primary = primary_cluster_x_fam - primary_fam,
      actual_non_primary = non_primary_fam - hidden_primary,
      .groups = "drop"
    ) 
  
  # Model fitting
  # ------------------------------------------------------
  # full mdl: major_lineage_family_lv + age_cat
  hidden_comm_acq_rate_aHH_full_mdl <- glm(
    cbind(hidden_primary, actual_non_primary) ~ major_lineage_family_lv + age_cat, # actual_non_primary = non_primary_fam - hidden_primary
    data = hidden_comm_acq_aHH_dat,
    family = binomial
  )
  
  # age mdl: age_cat
  hidden_comm_acq_rate_aHH_age_mdl <- glm(
    cbind(hidden_primary, actual_non_primary) ~ age_cat,
    data = hidden_comm_acq_aHH_dat,
    family = binomial
  )
  
  # Estimate rates
  # ------------------------------------------------------
  # full mdl: major_lineage_family_lv + age_cat
  hidden_comm_acq_rate_estimates_aHH_full_mdl <- emmeans(
    hidden_comm_acq_rate_aHH_full_mdl, 
    ~ major_lineage_family_lv + age_cat, 
    type = "link"
  ) %>% as.data.frame() %>%
    mutate(
      sim = sim,
      model_type = as.character(formula(hidden_comm_acq_rate_aHH_full_mdl))[3]
    )
  
  # age mdl: age_cat
  hidden_comm_acq_rate_estimates_aHH_age_mdl <- emmeans(
    hidden_comm_acq_rate_aHH_age_mdl,
    ~ age_cat, 
    type = "link"
  ) %>% as.data.frame() %>%
    mutate(
      sim = sim,
      model_type = as.character(formula(hidden_comm_acq_rate_aHH_age_mdl))[3]
    )

  # Compute p-value using 'car::Anova'
  # ------------------------------------------------------
  # full mdl: major_lineage_family_lv + age_cat
  hidden_comm_acq_rate_aHH_full_mdl_anova <- car::Anova(hidden_comm_acq_rate_aHH_full_mdl, test = "LR") %>% as.data.frame()
  hidden_comm_acq_rate_aHH_full_mdl_pvals <- data.frame(
    term = rownames(hidden_comm_acq_rate_aHH_full_mdl_anova),
    p_value = hidden_comm_acq_rate_aHH_full_mdl_anova$`Pr(>Chisq)`,
    sim = sim,
    model_type = as.character(formula(hidden_comm_acq_rate_aHH_full_mdl))[3]
  )
  
  # age mdl: age_cat
  hidden_comm_acq_rate_aHH_age_mdl_anova <- car::Anova(hidden_comm_acq_rate_aHH_age_mdl, test = "LR") %>% as.data.frame()
  hidden_comm_acq_rate_aHH_age_mdl_pvals <- data.frame(
    term = rownames(hidden_comm_acq_rate_aHH_age_mdl_anova),
    p_value = hidden_comm_acq_rate_aHH_age_mdl_anova$`Pr(>Chisq)`,
    sim = sim,
    model_type = as.character(formula(hidden_comm_acq_rate_aHH_age_mdl))[3]
  )
  
  # Log results
  # ------------------------------------------------------
  sim_hidden_comm_acq_rate_aHH_log_odds[[sim]] <- bind_rows(
    summary(hidden_comm_acq_rate_aHH_full_mdl)$coefficients %>% as.data.frame %>% 
      select("Estimate",  "Std. Error") %>% 
      rename("average" = "Estimate", "SE" = "Std. Error") %>% rownames_to_column(var = "variable") %>%
      mutate(
        sim = sim,
        model_type = as.character(formula(hidden_comm_acq_rate_aHH_full_mdl))[3]
      ),
    summary(hidden_comm_acq_rate_aHH_age_mdl)$coefficients %>% as.data.frame %>% 
      select("Estimate",  "Std. Error") %>% 
      rename("average" = "Estimate", "SE" = "Std. Error") %>% rownames_to_column(var = "variable") %>%
      mutate(
        sim = sim,
        model_type = as.character(formula(hidden_comm_acq_rate_aHH_age_mdl))[3]
      )
  )
  
  sim_hidden_comm_acq_rate_aHH_estimates[[sim]] <- bind_rows(
    hidden_comm_acq_rate_estimates_aHH_full_mdl, 
    hidden_comm_acq_rate_estimates_aHH_age_mdl
  )
  
  sim_hidden_comm_acq_rate_aHH_pvals[[sim]] <- bind_rows(
    hidden_comm_acq_rate_aHH_full_mdl_pvals, 
    hidden_comm_acq_rate_aHH_age_mdl_pvals
  )
  
  # ###############################
  # Compute corrected attack rates to using glm
  # ###############################
  # Prepare data
  # ------------------------------------------------------
  within_HH_subseq_dat <- participant_metadata_from_HH_with_complete_seq_dat %>%
    group_by(major_lineage_family_lv, first_HH_case_age_cat, family_id, age_cat) %>%
    summarise(
      secondary_cases = sum(overall_res_end == "Pos" & primary_cluster_x_fam == 0, na.rm = TRUE), # number of secondary_cases, excluding all community acquired cases (primary_cluster_x_fam == 1)
      failures = sum(overall_res_end == "Neg", na.rm = TRUE),
      .groups = "drop"
    ) %>% na.omit()
  
  # Model fitting
  # ------------------------------------------------------
  # full mdl: major_lineage_family_lv + age_cat + first_HH_case_age_cat
  atk_rate_full_mdl <- glm(
    cbind(secondary_cases, failures) ~ major_lineage_family_lv + age_cat + first_HH_case_age_cat, 
    data = within_HH_subseq_dat,
    family = binomial
  )
  
  # age_HHtype mdl: age_cat + first_HH_case_age_cat
  atk_rate_age_HHtype_mdl <- glm(
    cbind(secondary_cases, failures) ~ age_cat + first_HH_case_age_cat,
    data = within_HH_subseq_dat,
    family = binomial
  )
  
  # age mdl: age_cat
  atk_rate_age_mdl <- glm(
    cbind(secondary_cases, failures) ~ age_cat,
    data = within_HH_subseq_dat,
    family = binomial
  )
  
  # HHtype mdl: first_HH_case_age_cat
  atk_rate_HHtype_mdl <- glm(
    cbind(secondary_cases, failures) ~ first_HH_case_age_cat,
    data = within_HH_subseq_dat,
    family = binomial
  )
  
  # Estimate rates
  # ------------------------------------------------------
  # full mdl: major_lineage_family_lv + age_cat + first_HH_case_age_cat
  atk_rate_estimates_full_mdl <- emmeans(atk_rate_full_mdl, ~ major_lineage_family_lv + age_cat + first_HH_case_age_cat, type = "link") %>% as.data.frame() %>%
    mutate(
      sim = sim,
      model_type = as.character(formula(atk_rate_full_mdl))[3]
    )
  
  # age_HHtype mdl: age_cat + first_HH_case_age_cat
  atk_rate_estimates_age_HHtype_mdl <- bind_rows(
    emmeans(atk_rate_age_HHtype_mdl, ~ age_cat + first_HH_case_age_cat, type = "link") %>% as.data.frame() %>%
      mutate(
        sim = sim,
        model_type = as.character(formula(atk_rate_age_HHtype_mdl))[3]
      ),
    
    emmeans(atk_rate_age_HHtype_mdl, ~ first_HH_case_age_cat, type = "link", weights = "cells") %>% as.data.frame() %>%
      mutate(
        sim = sim,
        model_type = as.character(formula(atk_rate_age_HHtype_mdl))[3]
      ),
    
    emmeans(atk_rate_age_HHtype_mdl, ~ age_cat, type = "link", weights = "cells") %>% as.data.frame() %>%
      mutate(
        sim = sim,
        model_type = as.character(formula(atk_rate_age_HHtype_mdl))[3]
      )
  )
  
  # age mdl: age_cat
  atk_rate_estimates_age_mdl <- emmeans(atk_rate_age_mdl, ~ age_cat, type = "link") %>% as.data.frame() %>% # weights = "cells" has no practical effect because there are no other factors to average over >> removed
    mutate(
      sim = sim,
      model_type = as.character(formula(atk_rate_age_mdl))[3]
    )
  
  # HHtype mdl: first_HH_case_age_cat
  atk_rate_estimates_HHtype_mdl <- emmeans(atk_rate_HHtype_mdl, ~ first_HH_case_age_cat, type = "link") %>% as.data.frame() %>%
    mutate(
      sim = sim,
      model_type = as.character(formula(atk_rate_HHtype_mdl))[3]
    )
  
  # Compute p-values using "car::Anova"
  # ------------------------------------------------------
  # full mdl: major_lineage_family_lv + age_cat + first_HH_case_age_cat
  atk_rate_full_mdl_anova <- car::Anova(atk_rate_full_mdl, test = "LR") %>% as.data.frame()
  atk_rate_full_mdl_pvals <- data.frame(
    term = rownames(atk_rate_full_mdl_anova),
    p_value = atk_rate_full_mdl_anova$`Pr(>Chisq)`,
    sim = sim,
    model_type = as.character(formula(atk_rate_full_mdl))[3]
  )
  
  # age_HHtype mdl: age_cat + first_HH_case_age_cat
  atk_rate_age_HHtype_mdl_anova <- car::Anova(atk_rate_age_HHtype_mdl, test = "LR") %>% as.data.frame()
  atk_rate_age_HHtype_mdl_pvals <- data.frame(
    term = rownames(atk_rate_age_HHtype_mdl_anova),
    p_value = atk_rate_age_HHtype_mdl_anova$`Pr(>Chisq)`,
    sim = sim,
    model_type = as.character(formula(atk_rate_age_HHtype_mdl))[3]
  )
  
  # age mdl: age_cat
  atk_rate_age_mdl_anova <- car::Anova(atk_rate_age_mdl, test = "LR") %>% as.data.frame()
  atk_rate_age_mdl_pvals <- data.frame(
    term = rownames(atk_rate_age_mdl_anova),
    p_value = atk_rate_age_mdl_anova$`Pr(>Chisq)`,
    sim = sim,
    model_type = as.character(formula(atk_rate_age_mdl))[3]
  )
  
  # HHtype mdl: first_HH_case_age_cat
  atk_rate_HHtype_mdl_anova <- car::Anova(atk_rate_HHtype_mdl, test = "LR") %>% as.data.frame()
  atk_rate_HHtype_mdl_pvals <- data.frame(
    term = rownames(atk_rate_HHtype_mdl_anova),
    p_value = atk_rate_HHtype_mdl_anova$`Pr(>Chisq)`,
    sim = sim,
    model_type = as.character(formula(atk_rate_HHtype_mdl))[3]
  )

  # Log results
  # ------------------------------------------------------
  sim_atk_rate_log_odds[[sim]] <- bind_rows(
    summary(atk_rate_full_mdl)$coefficients %>% as.data.frame %>% 
      select("Estimate",  "Std. Error") %>% 
      rename("average" = "Estimate", "SE" = "Std. Error") %>% rownames_to_column(var = "variable") %>%
      mutate(
        sim = sim,
        model_type = as.character(formula(atk_rate_full_mdl))[3]
      ),
    summary(atk_rate_age_HHtype_mdl)$coefficients %>% as.data.frame %>% 
      select("Estimate",  "Std. Error") %>% 
      rename("average" = "Estimate", "SE" = "Std. Error") %>% rownames_to_column(var = "variable") %>%
      mutate(
        sim = sim,
        model_type = as.character(formula(atk_rate_age_HHtype_mdl))[3]
      ),
    summary(atk_rate_age_mdl)$coefficients %>% as.data.frame %>% 
      select("Estimate",  "Std. Error") %>% 
      rename("average" = "Estimate", "SE" = "Std. Error") %>% rownames_to_column(var = "variable") %>%
      mutate(
        sim = sim,
        model_type = as.character(formula(atk_rate_age_mdl))[3]
      ),
    summary(atk_rate_HHtype_mdl)$coefficients %>% as.data.frame %>% 
      select("Estimate",  "Std. Error") %>% 
      rename("average" = "Estimate", "SE" = "Std. Error") %>% rownames_to_column(var = "variable") %>%
      mutate(
        sim = sim,
        model_type = as.character(formula(atk_rate_HHtype_mdl))[3]
      )
  )
  
  sim_atk_rate_estimates[[sim]] <- bind_rows(
    atk_rate_estimates_full_mdl, 
    atk_rate_estimates_age_HHtype_mdl,
    atk_rate_estimates_age_mdl,
    atk_rate_estimates_HHtype_mdl
  )
  sim_atk_rate_pvals[[sim]] <- bind_rows(
    atk_rate_full_mdl_pvals, 
    atk_rate_age_HHtype_mdl_pvals,
    atk_rate_age_mdl_pvals,
    atk_rate_HHtype_mdl_pvals
  )
  
}

# ###############################
# Summarise results
# ###############################
# Hidden community acquisition rates using glm, focusing on adult-first households (aHH)
# -------------------------------
# p value distribution
bind_rows(sim_hidden_comm_acq_rate_aHH_pvals) %>% 
  group_by(model_type, term) %>%
  summarise(
    median_p = median(p_value, na.rm = TRUE),
    p_025 = quantile(p_value, 0.025, na.rm = TRUE),
    p_975 = quantile(p_value, 0.975, na.rm = TRUE),
    prop_signif = mean(p_value < 0.05, na.rm = TRUE),
    .groups = "drop"
  ) %>% as.data.frame()
#                          model_type                    term   median_p        p_025     p_975 prop_signif
# 1                           age_cat                 age_cat 0.09421782 0.0018322484 0.6232214       0.276
# 2 major_lineage_family_lv + age_cat                 age_cat 0.03658846 0.0006413995 0.4271064       0.504 ******************
# 3 major_lineage_family_lv + age_cat major_lineage_family_lv 0.11895445 0.0504413092 0.1787449       0.000 ******************

# Rate estimates
bind_rows(sim_hidden_comm_acq_rate_aHH_estimates) %>% 
  mutate(model_type = factor(model_type)) %>%
  group_by(model_type, age_cat, major_lineage_family_lv) %>%
  summarise(
    N = n(),
    link_bar = mean(emmean),               # mean logit
    W_bar = mean(SE^2),                    # avg within-run var
    B = var(emmean),                       # between-run var
    T_var = W_bar + (1 + 1/N) * B,        # total variance
    df = (N - 1) * (1 + W_bar / ((1 + 1/N) * B))^2,
    lower95 = link_bar - qt(0.975, df) * sqrt(T_var),
    upper95 = link_bar + qt(0.975, df) * sqrt(T_var),
    .groups = "drop"
  ) %>%
  mutate(
    mean_prob = round(plogis(link_bar) * 100, 2),
    lower95 = round(plogis(lower95) * 100, 2),
    upper95 = round(plogis(upper95) * 100, 2)
  ) %>%
  select(
    age_cat, major_lineage_family_lv, model_type,
    N,
    mean_prob, lower95, upper95) %>% 
  arrange(model_type) %>% as.data.frame()
#   age_cat major_lineage_family_lv                        model_type    N mean_prob lower95 upper95
# 1   Adult                    <NA>                           age_cat 1000      0.10    0.00  100.00 ******************
# 2   Child                    <NA>                           age_cat 1000     28.20   10.02   58.09 ******************
# 3   Adult                 BA.2.75 major_lineage_family_lv + age_cat 1000      0.00    0.00  100.00
# 4   Adult                 BA.2.86 major_lineage_family_lv + age_cat 1000      0.08    0.00  100.00
# 5   Adult              XBB or XBL major_lineage_family_lv + age_cat 1000      0.11    0.00  100.00
# 6   Child                 BA.2.75 major_lineage_family_lv + age_cat 1000      0.00    0.00  100.00
# 7   Child                 BA.2.86 major_lineage_family_lv + age_cat 1000     38.17    7.93   81.57
# 8   Child              XBB or XBL major_lineage_family_lv + age_cat 1000     45.35   12.77   82.47

# log-odds diff estimates
bind_rows(sim_hidden_comm_acq_rate_aHH_log_odds) %>% 
  mutate(model_type = factor(model_type)) %>%
  group_by(model_type, variable) %>%
  summarise(
    N = n(),
    link_bar = mean(average),               # mean logit
    W_bar = mean(SE^2),                    # avg within-run var
    B = var(average),                       # between-run var
    T_var = W_bar + (1 + 1/N) * B,        # total variance
    df = (N - 1) * (1 + W_bar / ((1 + 1/N) * B))^2,
    lower95 = link_bar - qt(0.975, df) * sqrt(T_var),
    upper95 = link_bar + qt(0.975, df) * sqrt(T_var),
    .groups = "drop"
  ) %>%
  mutate(
    average = round(link_bar, 2),
    lower95 = round(lower95, 2),
    upper95 = round(upper95, 2)
  ) %>%
  select(
    variable, model_type,
    N,
    average, lower95, upper95) %>% 
  arrange(model_type) %>% as.data.frame()
#                            variable                        model_type    N average   lower95  upper95
# 1                       (Intercept)                           age_cat 1000   -6.87  -3173.30  3159.56
# 2                      age_catChild                           age_cat 1000    5.93  -3160.50  3172.37
# 3                       (Intercept) major_lineage_family_lv + age_cat 1000  -25.76 -10240.74 10189.22
# 4                      age_catChild major_lineage_family_lv + age_cat 1000    6.59  -4838.56  4851.74
# 5    major_lineage_family_lvBA.2.86 major_lineage_family_lv + age_cat 1000   18.69  -8974.10  9011.48
# 6 major_lineage_family_lvXBB or XBL major_lineage_family_lv + age_cat 1000   18.99  -8973.81  9011.78

# Secondary attack rates using glm
# -------------------------------
# p value distribution
bind_rows(sim_atk_rate_pvals) %>% 
  group_by(model_type, term) %>%
  summarise(
    median_p = median(p_value, na.rm = TRUE),
    p_025 = quantile(p_value, 0.025, na.rm = TRUE),
    p_975 = quantile(p_value, 0.975, na.rm = TRUE),
    prop_signif = mean(p_value < 0.05, na.rm = TRUE),
    .groups = "drop"
  ) %>% as.data.frame()
#                                                  model_type                    term     median_p        p_025        p_975 prop_signif
# 1                                                   age_cat                 age_cat 3.914106e-02 7.562692e-03 1.567068e-01        0.71
# 2                           age_cat + first_HH_case_age_cat                 age_cat 1.223465e-01 5.763385e-02 2.398557e-01        0.00
# 3                           age_cat + first_HH_case_age_cat   first_HH_case_age_cat 4.367341e-08 2.300829e-08 1.029996e-07        1.00
# 4                                     first_HH_case_age_cat   first_HH_case_age_cat 1.668255e-08 1.668255e-08 1.668255e-08        1.00 
# 5 major_lineage_family_lv + age_cat + first_HH_case_age_cat                 age_cat 1.400062e-01 6.116028e-02 2.686246e-01        0.00 ******************
# 6 major_lineage_family_lv + age_cat + first_HH_case_age_cat   first_HH_case_age_cat 1.845332e-07 9.632190e-08 4.413987e-07        1.00 ******************
# 7 major_lineage_family_lv + age_cat + first_HH_case_age_cat major_lineage_family_lv 5.269078e-01 4.778353e-01 5.618540e-01        0.00 ******************

# log-odds diff estimates
bind_rows(sim_atk_rate_log_odds) %>% 
  mutate(model_type = factor(model_type)) %>%
  group_by(model_type, variable) %>%
  summarise(
    N = n(),
    link_bar = mean(average),              # mean logit
    W_bar = mean(SE^2),                    # avg within-run var
    B = var(average),                      # between-run var
    T_var = W_bar + (1 + 1/N) * B,         # total variance
    df = (N - 1) * (1 + W_bar / ((1 + 1/N) * B))^2,
    lower95 = link_bar - qt(0.975, df) * sqrt(T_var),
    upper95 = link_bar + qt(0.975, df) * sqrt(T_var),
    .groups = "drop"
  ) %>%
  mutate(
    coeff_average = round(link_bar, 2),
    coeff_lower95 = round(lower95, 2),
    coeff_upper95 = round(upper95, 2),
    
    exp_coeff_average = round(exp(link_bar), 2),
    exp_coeff_lower95 = round(exp(lower95), 2),
    exp_coeff_upper95 = round(exp(upper95), 2)
  ) %>%
  select(
    variable, model_type,
    N,
    
    coeff_average, coeff_lower95, coeff_upper95,
    exp_coeff_average, exp_coeff_lower95, exp_coeff_upper95
    ) %>% 
  arrange(model_type) %>% as.data.frame()
#                                      variable                                                model_type    N coeff_average coeff_lower95 coeff_upper95 exp_coeff_average exp_coeff_lower95 exp_coeff_upper95
# 1                                 (Intercept)                                                   age_cat 1000         -0.89         -1.41         -0.37              0.41              0.24              0.69
# 2                                age_catChild                                                   age_cat 1000          1.12          0.02          2.22              3.06              1.02              9.21
# 3                                 (Intercept)                           age_cat + first_HH_case_age_cat 1000         -2.27         -3.20         -1.34              0.10              0.04              0.26
# 4                                age_catChild                           age_cat + first_HH_case_age_cat 1000          1.03         -0.30          2.37              2.81              0.74             10.65
# 5  first_HH_case_age_catAdult-first household                           age_cat + first_HH_case_age_cat 1000          2.76          1.65          3.88             15.88              5.21             48.40
# 6                                 (Intercept)                                     first_HH_case_age_cat 1000         -2.06         -2.91         -1.21              0.13              0.05              0.30 ******************
# 7  first_HH_case_age_catAdult-first household                                     first_HH_case_age_cat 1000          2.79          1.70          3.89             16.32              5.47             48.71 ******************
# 8                                 (Intercept) major_lineage_family_lv + age_cat + first_HH_case_age_cat 1000         -1.70         -3.07         -0.34              0.18              0.05              0.71
# 9                                age_catChild major_lineage_family_lv + age_cat + first_HH_case_age_cat 1000          1.06         -0.36          2.47              2.88              0.70             11.82
# 10 first_HH_case_age_catAdult-first household major_lineage_family_lv + age_cat + first_HH_case_age_cat 1000          2.68          1.55          3.82             14.64              4.72             45.43
# 11             major_lineage_family_lvBA.2.86 major_lineage_family_lv + age_cat + first_HH_case_age_cat 1000         -0.42         -1.96          1.12              0.66              0.14              3.08
# 12                major_lineage_family_lvBA.5 major_lineage_family_lv + age_cat + first_HH_case_age_cat 1000        -16.65      -4751.32       4718.02              0.00              0.00               Inf
# 13          major_lineage_family_lvXBB or XBL major_lineage_family_lv + age_cat + first_HH_case_age_cat 1000         -0.75         -2.20          0.70              0.47              0.11              2.01

# Rate estimates
bind_rows(sim_atk_rate_estimates) %>% 
  mutate(model_type = factor(model_type)) %>%
  group_by(model_type, first_HH_case_age_cat, age_cat, major_lineage_family_lv) %>%
  summarise(
    N = n(),
    link_bar = mean(emmean),               # mean logit
    W_bar = mean(SE^2),                    # avg within-run var
    B = var(emmean),                       # between-run var
    T_var = W_bar + (1 + 1/N) * B,         # total variance
    df = (N - 1) * (1 + W_bar / ((1 + 1/N) * B))^2,
    lower95 = link_bar - qt(0.975, df) * sqrt(T_var),
    upper95 = link_bar + qt(0.975, df) * sqrt(T_var),
    .groups = "drop"
  ) %>%
  mutate(
    mean_prob = round(plogis(link_bar) * 100, 2),
    lower95 = round(plogis(lower95) * 100, 2),
    upper95 = round(plogis(upper95) * 100, 2)
  ) %>%
  select(
    first_HH_case_age_cat, age_cat, major_lineage_family_lv, model_type,
    N,
    mean_prob, lower95, upper95) %>% 
  arrange(model_type) %>% as.data.frame()
#    first_HH_case_age_cat age_cat major_lineage_family_lv                                                model_type    N mean_prob lower95 upper95
# 1                   <NA>   Adult                    <NA>                                                   age_cat 1000     29.05   19.61   40.75
# 2                   <NA>   Child                    <NA>                                                   age_cat 1000     55.65   32.64   76.46
# 3  Child-first household   Adult                    <NA>                           age_cat + first_HH_case_age_cat 1000      9.34    3.90   20.70
# 4  Child-first household   Child                    <NA>                           age_cat + first_HH_case_age_cat 1000     22.44    7.30   51.51
# 5  Child-first household    <NA>                    <NA>                           age_cat + first_HH_case_age_cat 1000     10.74    4.79   22.38
# 6  Adult-first household   Adult                    <NA>                           age_cat + first_HH_case_age_cat 1000     62.05   43.59   77.58
# 7  Adult-first household   Child                    <NA>                           age_cat + first_HH_case_age_cat 1000     82.12   56.00   94.31
# 8  Adult-first household    <NA>                    <NA>                           age_cat + first_HH_case_age_cat 1000     68.48   51.64   81.56
# 9                   <NA>   Adult                    <NA>                           age_cat + first_HH_case_age_cat 1000     22.47   12.74   36.51
# 10                  <NA>   Child                    <NA>                           age_cat + first_HH_case_age_cat 1000     57.38   29.08   81.55
# 11 Child-first household    <NA>                    <NA>                                     first_HH_case_age_cat 1000     11.32    5.18   22.99 ******************
# 12 Adult-first household    <NA>                    <NA>                                     first_HH_case_age_cat 1000     67.57   51.14   80.57 ******************
# 13 Child-first household   Adult                 BA.2.75 major_lineage_family_lv + age_cat + first_HH_case_age_cat 1000     15.39    4.43   41.67
# 14 Child-first household   Adult                 BA.2.86 major_lineage_family_lv + age_cat + first_HH_case_age_cat 1000     10.69    3.46   28.54
# 15 Child-first household   Adult                    BA.5 major_lineage_family_lv + age_cat + first_HH_case_age_cat 1000      0.00    0.00  100.00
# 16 Child-first household   Adult              XBB or XBL major_lineage_family_lv + age_cat + first_HH_case_age_cat 1000      7.89    2.64   21.26
# 17 Child-first household   Child                 BA.2.75 major_lineage_family_lv + age_cat + first_HH_case_age_cat 1000     34.35    9.09   73.26
# 18 Child-first household   Child                 BA.2.86 major_lineage_family_lv + age_cat + first_HH_case_age_cat 1000     25.61    6.07   64.72
# 19 Child-first household   Child                    BA.5 major_lineage_family_lv + age_cat + first_HH_case_age_cat 1000      0.00    0.00  100.00
# 20 Child-first household   Child              XBB or XBL major_lineage_family_lv + age_cat + first_HH_case_age_cat 1000     19.76    5.17   52.65
# 21 Adult-first household   Adult                 BA.2.75 major_lineage_family_lv + age_cat + first_HH_case_age_cat 1000     72.70   40.63   91.20
# 22 Adult-first household   Adult                 BA.2.86 major_lineage_family_lv + age_cat + first_HH_case_age_cat 1000     63.67   37.54   83.63
# 23 Adult-first household   Adult                    BA.5 major_lineage_family_lv + age_cat + first_HH_case_age_cat 1000      0.00    0.00  100.00
# 24 Adult-first household   Adult              XBB or XBL major_lineage_family_lv + age_cat + first_HH_case_age_cat 1000     55.63   32.45   76.59
# 25 Adult-first household   Child                 BA.2.75 major_lineage_family_lv + age_cat + first_HH_case_age_cat 1000     88.45   58.07   97.70
# 26 Adult-first household   Child                 BA.2.86 major_lineage_family_lv + age_cat + first_HH_case_age_cat 1000     83.45   49.35   96.31
# 27 Adult-first household   Child                    BA.5 major_lineage_family_lv + age_cat + first_HH_case_age_cat 1000      0.00    0.00  100.00
# 28 Adult-first household   Child              XBB or XBL major_lineage_family_lv + age_cat + first_HH_case_age_cat 1000     78.29   45.99   93.85
