# ###############################
# Load global variables
# ###############################
source(file.path("miscellaneous.R"))

# ###############################
# Load required libraries
# ###############################
library(tidyverse)
library(gtsummary)
library(rstatix)

# ###############################
# Helper functions
# ###############################
# A function to make HH summary table
# -------------------------------
make_HH_sum_table_func <- function(HH_metadata, missing = "ifany"){
  HH_sum_table <- HH_metadata %>%
    tbl_summary(
      by = first_HH_case_age_cat,
      missing = missing, missing_text = "# of households not included",
      type = list(where(is.numeric) ~ "continuous"),
      statistic = list(all_continuous() ~ "{mean} ({sd}; {min}-{max})"),
      digits = all_continuous() ~ 2,
      label = label_list
    ) %>%
    add_overall() %>%
    add_n() %>%
    modify_header(
      all_stat_cols() ~ "**{ifelse(level == 'Child', 'Child-first', ifelse(level == 'Adult', 'Adult-first', ifelse(level == 'Adult & child', 'Adult-and-child co-first', 'Overall')))}**\n(N = {n})",
      label = "**Variable**"
      ) %>%
    bold_labels() %>%
    add_p(pvalue_fun = label_style_pvalue(digits = 3))
  return(HH_sum_table)
}

replace_pvalues <- function(tbl, posthoc_comp_res){
  tbl$table_body <- tbl$table_body %>%
    left_join(
      posthoc_comp_res %>%
        select(variable, p.adj) %>%
        mutate(row_type = "label"),
      by = c("variable", "row_type")
    ) %>%
    select(-any_of("p.value")) %>%
    rename(p.value = p.adj)
  
  tbl %>%
    modify_footnote_header(
      footnote = "Pairwise Wilcoxon rank sum test; Pairwise Fisher's exact test",
      columns = "p.value",
      replace = TRUE,
      text_interpret = c("md", "html")
    ) %>%
    modify_column_hide(columns = c(n, stat_0, stat_1, stat_2))
}

# ###############################
# Prepare household metadata
# ###############################
# Load data and format it a bit
# -------------------------------
HH_metadata <- read.delim(HH_metadata_filepath, header = TRUE, sep = "\t") %>%
  mutate(first_HH_case_age_cat = factor(first_HH_case_age_cat, c("Child", "Adult", "Adult & child"))) 

# Select relevant variables
# -------------------------------
HH_metadata <- HH_metadata %>%
  select(
    # family_id,
    # province, district,
    # enrolment_epiyear, enrolment_epiweek,
    first_HH_case_age_cat,
    
    # Overall measures
    ##############################
    HH_size,
    
    # num_children, num_adults,
    prop_children, prop_adults,
    child_to_adult_ratio,
    
    # num_males, num_females,
    prop_males, prop_females,
    male_to_female_ratio,
    
    mean_age,
    mean_BMI,
    
    # num_any_comorbidity,
    # num_pre_existing_med_cond_hd,
    # num_pre_existing_med_cond_ht,
    # num_pre_existing_med_cond_ld,
    # num_pre_existing_med_cond_lim,
    # num_pre_existing_med_cond_kd,
    # num_pre_existing_med_cond_bd,
    # num_pre_existing_med_cond_gi,
    # num_pre_existing_med_cond_im,
    # num_pre_existing_med_cond_dm,
    # num_pre_existing_med_cond_ar,
    # num_pre_existing_med_cond_tb,
    # num_pre_existing_med_cond_ml,
    # num_pre_existing_med_cond_al,
    # num_pre_existing_med_cond_hl,
    # num_pre_existing_med_cond_as,
    prop_any_comorbidity,
    prop_pre_existing_med_cond_hd,
    prop_pre_existing_med_cond_ht,
    prop_pre_existing_med_cond_ld,
    prop_pre_existing_med_cond_lim,
    prop_pre_existing_med_cond_kd,
    prop_pre_existing_med_cond_bd,
    prop_pre_existing_med_cond_gi,
    prop_pre_existing_med_cond_im,
    prop_pre_existing_med_cond_dm,
    prop_pre_existing_med_cond_ar,
    prop_pre_existing_med_cond_tb,
    prop_pre_existing_med_cond_ml,
    prop_pre_existing_med_cond_al,
    prop_pre_existing_med_cond_hl,
    prop_pre_existing_med_cond_as,
    
    # num_pos_enrol,
    # num_neg_enrol,
    prop_pos_enrol,
    # num_pos_end,
    # num_neg_end,
    # num_unconf_end,
    prop_pos_end,
    # num_pos_add,
    prop_pos_add,

    mean_most_severe_fv,
    mean_most_severe_co,
    mean_most_severe_st,
    mean_most_severe_rn,
    mean_most_severe_db,
    mean_most_severe_he,
    mean_most_severe_my,
    mean_most_severe_ch,
    mean_most_severe_fa,
    mean_most_severe_di,
    mean_most_severe_lo,
    mean_most_severe_vo,
    mean_most_severe_ho,
    max_most_severe_fv,
    max_most_severe_co,
    max_most_severe_st,
    max_most_severe_rn,
    max_most_severe_db,
    max_most_severe_he,
    max_most_severe_my,
    max_most_severe_ch,
    max_most_severe_fa,
    max_most_severe_di,
    max_most_severe_lo,
    max_most_severe_vo,
    max_most_severe_ho,

    # num_min_symp_cases,
    # num_mod_symp_cases,
    # num_highly_symp_cases,
    # num_no_symp_data_cases,
    prop_min_symp_cases,
    prop_mod_symp_cases,
    prop_highly_symp_cases,
    prop_no_symp_data_cases,

    mean_pcr_orf1ab_ct,
    mean_pcr_n_ct,
    mean_pcr_e_ct,
    min_pcr_orf1ab_ct,
    min_pcr_n_ct,
    min_pcr_e_ct,

    mean_num_prev_COVID,
    # num_prev_COVID_lt_6mo,
    # num_prev_COVID_6_12mo,
    # num_prev_COVID_gt_12mo,
    # num_prev_COVID_never,
    # num_prev_COVID_unknown,
    prop_prev_COVID_lt_6mo,
    prop_prev_COVID_6_12mo,
    prop_prev_COVID_gt_12mo,
    prop_prev_COVID_never,
    prop_prev_COVID_unknown,

    mean_num_vac,
    mean_num_vac_mRNA,
    mean_num_vac_inact,
    mean_num_vac_vector,
    mean_num_vac_other,
    # num_last_vac_lt_6mo,
    # num_last_vac_6_12mo,
    # num_last_vac_gt_12mo,
    # num_last_vac_never,
    # num_last_vac_unknown,
    prop_last_vac_lt_6mo,
    prop_last_vac_6_12mo,
    prop_last_vac_gt_12mo,
    prop_last_vac_never,
    prop_last_vac_unknown,
    
    # num_seq_avail,
    prop_seq_avail,
    major_lineage,

    # Children-specific measures
    ##############################
    prop_males_children,
    prop_females_children,
    male_to_female_ratio_children,

    mean_age_children,
    mean_BMI_children,
    
    # num_any_comorbidity_children,
    # num_pre_existing_med_cond_hd_children,
    # num_pre_existing_med_cond_ht_children,
    # num_pre_existing_med_cond_ld_children,
    # num_pre_existing_med_cond_lim_children,
    # num_pre_existing_med_cond_kd_children,
    # num_pre_existing_med_cond_bd_children,
    # num_pre_existing_med_cond_gi_children,
    # num_pre_existing_med_cond_im_children,
    # num_pre_existing_med_cond_dm_children,
    # num_pre_existing_med_cond_ar_children,
    # num_pre_existing_med_cond_tb_children,
    # num_pre_existing_med_cond_ml_children,
    # num_pre_existing_med_cond_al_children,
    # num_pre_existing_med_cond_hl_children,
    # num_pre_existing_med_cond_as_children,
    prop_any_comorbidity_children,
    prop_pre_existing_med_cond_hd_children,
    prop_pre_existing_med_cond_ht_children,
    prop_pre_existing_med_cond_ld_children,
    prop_pre_existing_med_cond_lim_children,
    prop_pre_existing_med_cond_kd_children,
    prop_pre_existing_med_cond_bd_children,
    prop_pre_existing_med_cond_gi_children,
    prop_pre_existing_med_cond_im_children,
    prop_pre_existing_med_cond_dm_children,
    prop_pre_existing_med_cond_ar_children,
    prop_pre_existing_med_cond_tb_children,
    prop_pre_existing_med_cond_ml_children,
    prop_pre_existing_med_cond_al_children,
    prop_pre_existing_med_cond_hl_children,
    prop_pre_existing_med_cond_as_children,
    
    # num_pos_enrol_children,
    # num_neg_enrol_children,
    prop_pos_enrol_children,
    # num_pos_end_children,
    # num_neg_end_children,
    # num_unconf_end_children,
    prop_pos_end_children,
    # num_pos_add_children,
    prop_pos_add_children,
    
    mean_most_severe_fv_children,
    mean_most_severe_co_children,
    mean_most_severe_st_children,
    mean_most_severe_rn_children,
    mean_most_severe_db_children,
    mean_most_severe_he_children,
    mean_most_severe_my_children,
    mean_most_severe_ch_children,
    mean_most_severe_fa_children,
    mean_most_severe_di_children,
    mean_most_severe_lo_children,
    mean_most_severe_vo_children,
    mean_most_severe_ho_children,
    max_most_severe_fv_children,
    max_most_severe_co_children,
    max_most_severe_st_children,
    max_most_severe_rn_children,
    max_most_severe_db_children,
    max_most_severe_he_children,
    max_most_severe_my_children,
    max_most_severe_ch_children,
    max_most_severe_fa_children,
    max_most_severe_di_children,
    max_most_severe_lo_children,
    max_most_severe_vo_children,
    max_most_severe_ho_children,
    
    # num_min_symp_cases_children,
    # num_mod_symp_cases_children,
    # num_highly_symp_cases_children,
    # num_no_symp_data_cases_children,
    prop_min_symp_cases_children,
    prop_mod_symp_cases_children,
    prop_highly_symp_cases_children,
    prop_no_symp_data_cases_children,
    
    mean_pcr_orf1ab_ct_children,
    mean_pcr_n_ct_children,
    mean_pcr_e_ct_children,
    mean_num_prev_COVID_children,
    min_pcr_orf1ab_ct_children,
    min_pcr_n_ct_children,
    min_pcr_e_ct_children,
    
    # num_prev_COVID_lt_6mo_children,
    # num_prev_COVID_6_12mo_children,
    # num_prev_COVID_gt_12mo_children,
    # num_prev_COVID_never_children,
    # num_prev_COVID_unknown_children,
    prop_prev_COVID_lt_6mo_children,
    prop_prev_COVID_6_12mo_children,
    prop_prev_COVID_gt_12mo_children,
    prop_prev_COVID_never_children,
    prop_prev_COVID_unknown_children,
    
    mean_num_vac_children,
    mean_num_vac_mRNA_children,
    mean_num_vac_inact_children,
    mean_num_vac_vector_children,
    mean_num_vac_other_children,
    # num_last_vac_lt_6mo_children,
    # num_last_vac_6_12mo_children,
    # num_last_vac_gt_12mo_children,
    # num_last_vac_never_children,
    # num_last_vac_unknown_children,
    prop_last_vac_lt_6mo_children,
    prop_last_vac_6_12mo_children,
    prop_last_vac_gt_12mo_children,
    prop_last_vac_never_children,
    prop_last_vac_unknown_children,
    
    # num_seq_avail_children,
    prop_seq_avail_children,
    
    # Adult-specific measures
    ##############################
    prop_males_adults,
    prop_females_adults,
    male_to_female_ratio_adults,
    
    mean_age_adults,
    mean_BMI_adults,
    
    # num_any_comorbidity_adults,
    # num_pre_existing_med_cond_hd_adults,
    # num_pre_existing_med_cond_ht_adults,
    # num_pre_existing_med_cond_ld_adults,
    # num_pre_existing_med_cond_lim_adults,
    # num_pre_existing_med_cond_kd_adults,
    # num_pre_existing_med_cond_bd_adults,
    # num_pre_existing_med_cond_gi_adults,
    # num_pre_existing_med_cond_im_adults,
    # num_pre_existing_med_cond_dm_adults,
    # num_pre_existing_med_cond_ar_adults,
    # num_pre_existing_med_cond_tb_adults,
    # num_pre_existing_med_cond_ml_adults,
    # num_pre_existing_med_cond_al_adults,
    # num_pre_existing_med_cond_hl_adults,
    # num_pre_existing_med_cond_as_adults,
    prop_any_comorbidity_adults,
    prop_pre_existing_med_cond_hd_adults,
    prop_pre_existing_med_cond_ht_adults,
    prop_pre_existing_med_cond_ld_adults,
    prop_pre_existing_med_cond_lim_adults,
    prop_pre_existing_med_cond_kd_adults,
    prop_pre_existing_med_cond_bd_adults,
    prop_pre_existing_med_cond_gi_adults,
    prop_pre_existing_med_cond_im_adults,
    prop_pre_existing_med_cond_dm_adults,
    prop_pre_existing_med_cond_ar_adults,
    prop_pre_existing_med_cond_tb_adults,
    prop_pre_existing_med_cond_ml_adults,
    prop_pre_existing_med_cond_al_adults,
    prop_pre_existing_med_cond_hl_adults,
    prop_pre_existing_med_cond_as_adults,
    
    # num_pos_enrol_adults,
    # num_neg_enrol_adults,
    prop_pos_enrol_adults,
    # num_pos_end_adults,
    # num_neg_end_adults,
    # num_unconf_end_adults,
    prop_pos_end_adults,
    # num_pos_add_adults,
    prop_pos_add_adults,

    mean_most_severe_fv_adults,
    mean_most_severe_co_adults,
    mean_most_severe_st_adults,
    mean_most_severe_rn_adults,
    mean_most_severe_db_adults,
    mean_most_severe_he_adults,
    mean_most_severe_my_adults,
    mean_most_severe_ch_adults,
    mean_most_severe_fa_adults,
    mean_most_severe_di_adults,
    mean_most_severe_lo_adults,
    mean_most_severe_vo_adults,
    mean_most_severe_ho_adults,
    max_most_severe_fv_adults,
    max_most_severe_co_adults,
    max_most_severe_st_adults,
    max_most_severe_rn_adults,
    max_most_severe_db_adults,
    max_most_severe_he_adults,
    max_most_severe_my_adults,
    max_most_severe_ch_adults,
    max_most_severe_fa_adults,
    max_most_severe_di_adults,
    max_most_severe_lo_adults,
    max_most_severe_vo_adults,
    max_most_severe_ho_adults,
    
    # num_min_symp_cases_adults,
    # num_mod_symp_cases_adults,
    # num_highly_symp_cases_adults,
    # num_no_symp_data_cases_adults,
    prop_min_symp_cases_adults,
    prop_mod_symp_cases_adults,
    prop_highly_symp_cases_adults,
    prop_no_symp_data_cases_adults,
    
    mean_pcr_orf1ab_ct_adults,
    mean_pcr_n_ct_adults,
    mean_pcr_e_ct_adults,
    min_pcr_orf1ab_ct_adults,
    min_pcr_n_ct_adults,
    min_pcr_e_ct_adults,
    
    mean_num_prev_COVID_adults,
    # num_prev_COVID_lt_6mo_adults,
    # num_prev_COVID_6_12mo_adults,
    # num_prev_COVID_gt_12mo_adults,
    # num_prev_COVID_never_adults,
    # num_prev_COVID_unknown_adults,
    prop_prev_COVID_lt_6mo_adults,
    prop_prev_COVID_6_12mo_adults,
    prop_prev_COVID_gt_12mo_adults,
    prop_prev_COVID_never_adults,
    prop_prev_COVID_unknown_adults,
    
    mean_num_vac_adults,
    mean_num_vac_mRNA_adults,
    mean_num_vac_inact_adults,
    mean_num_vac_vector_adults,
    mean_num_vac_other_adults,
    # num_last_vac_lt_6mo_adults,
    # num_last_vac_6_12mo_adults,
    # num_last_vac_gt_12mo_adults,
    # num_last_vac_never_adults,
    # num_last_vac_unknown_adults,
    prop_last_vac_lt_6mo_adults,
    prop_last_vac_6_12mo_adults,
    prop_last_vac_gt_12mo_adults,
    prop_last_vac_never_adults,
    prop_last_vac_unknown_adults,
    
    # num_seq_avail_adults,
    prop_seq_avail_adults,
  ) %>%
  mutate(across(starts_with("prop_"), ~ .x * 100)) %>%
  mutate(across(where(is.integer), as.numeric)) 

# Set labels
# -------------------------------
label_list <- list(
  # =============================
  # Overall measures
  # =============================
  # General demographic variables
  HH_size = "Number of household members",
  prop_children = "Percentage of children in the household",
  prop_adults = "Percentage of adults in the household",
  child_to_adult_ratio = "Child-to-adult ratio",
  
  prop_males = "Percentage of males in the household",
  prop_females = "Percentage of females in the household",
  male_to_female_ratio = "Male-to-female ratio",
  
  mean_age = "Mean age",
  mean_BMI = "Mean BMI†",
  
  # Comorbidities
  # † among household members with available information for that variable
  prop_any_comorbidity = "Percentage with any pre-existing medical condition†",
  prop_pre_existing_med_cond_hd = "Percentage with heart disease†",
  prop_pre_existing_med_cond_ht = "Percentage with hypertension†",
  prop_pre_existing_med_cond_ld = "Percentage with lung disease†",
  prop_pre_existing_med_cond_lim = "Percentage with low immunity†",
  prop_pre_existing_med_cond_kd = "Percentage with kidney disease†",
  prop_pre_existing_med_cond_bd = "Percentage with brain disorder†",
  prop_pre_existing_med_cond_gi = "Percentage with gastrointestinal disease†",
  prop_pre_existing_med_cond_im = "Percentage with immunodeficiency†",
  prop_pre_existing_med_cond_dm = "Percentage with diabetes†",
  prop_pre_existing_med_cond_ar = "Percentage with arthritis†",
  prop_pre_existing_med_cond_tb = "Percentage with tuberculosis†",
  prop_pre_existing_med_cond_ml = "Percentage with malaria†",
  prop_pre_existing_med_cond_al = "Percentage with allergic rhinitis†",
  prop_pre_existing_med_cond_hl = "Percentage with hyperlipidaemia†",
  prop_pre_existing_med_cond_as = "Percentage with asthma†",
  
  # Infection
  prop_pos_enrol = "Overall baseline case rate at enrolment (percentage of baseline cases at enrolment per household)",
  prop_pos_end = "Overall endpoint case rate (percentage of cases identified by the end of follow-up per household)",
  prop_pos_add = "Incident rate (crude attack rate; percentage of incident cases among participants at risk at enrolment)",
  
  # Symptom severity scores 
  # ‡ among household cases identified by the end of follow-up with available information for that variable
  mean_most_severe_fv = "Mean maximum fever severity score‡",
  mean_most_severe_co = "Mean maximum cough severity score‡",
  mean_most_severe_st = "Mean maximum sore throat severity score‡",
  mean_most_severe_rn = "Mean maximum runny nose severity score‡",
  mean_most_severe_db = "Mean maximum difficulty breathing severity score‡",
  mean_most_severe_he = "Mean maximum headache severity score‡",
  mean_most_severe_my = "Mean maximum myalgia severity score‡",
  mean_most_severe_ch = "Mean maximum chills severity score‡",
  mean_most_severe_fa = "Mean maximum fatigue severity score‡",
  mean_most_severe_di = "Mean maximum diarrhoea severity score‡",
  mean_most_severe_lo = "Mean maximum loss of smell/taste severity score‡",
  mean_most_severe_vo = "Mean maximum vomiting severity score‡",
  mean_most_severe_ho = "Mean maximum hoarseness severity score‡",
  
  max_most_severe_fv = "Highest fever severity score‡",
  max_most_severe_co = "Highest cough severity score‡",
  max_most_severe_st = "Highest sore throat severity score‡",
  max_most_severe_rn = "Highest runny nose severity score‡",
  max_most_severe_db = "Highest difficulty breathing severity score‡",
  max_most_severe_he = "Highest headache severity score‡",
  max_most_severe_my = "Highest myalgia severity score‡",
  max_most_severe_ch = "Highest chills severity score‡",
  max_most_severe_fa = "Highest fatigue severity score‡",
  max_most_severe_di = "Highest diarrhoea severity score‡",
  max_most_severe_lo = "Highest loss of smell/taste severity score‡",
  max_most_severe_vo = "Highest vomiting severity score‡",
  max_most_severe_ho = "Highest hoarseness severity score‡",
  
  # Symptom groups
  prop_min_symp_cases = "Percentage of minimally/(a)symptomatic cases‡",
  prop_mod_symp_cases = "Percentage of moderately symptomatic cases‡",
  prop_highly_symp_cases = "Percentage of highly symptomatic cases‡",
  prop_no_symp_data_cases = "Percentage of cases without symptom data",
  
  # PCR Ct values
  # ‡ among household cases identified by the end of follow-up with available information
  mean_pcr_orf1ab_ct = "Mean ORF1ab Ct value‡",
  mean_pcr_n_ct = "Mean N gene Ct value‡",
  mean_pcr_e_ct = "Mean E gene Ct value‡",
  
  min_pcr_orf1ab_ct = "Lowest ORF1ab Ct value‡",
  min_pcr_n_ct = "Lowest N gene Ct value‡",
  min_pcr_e_ct = "Lowest E gene Ct value‡",
  
  # Previous infection history
  mean_num_prev_COVID = "Mean number of previous COVID-19 infections†",
  prop_prev_COVID_lt_6mo = "Percentage with last infection within 6 months†",
  prop_prev_COVID_6_12mo = "Percentage with last infection 6–12 months earlier†",
  prop_prev_COVID_gt_12mo = "Percentage with last infection more than 12 months earlier†",
  prop_prev_COVID_never = "Percentage with no previous COVID-19 infection†",
  prop_prev_COVID_unknown = "Percentage with unknown previous COVID-19 history",
  
  # Vaccination
  mean_num_vac = "Mean number of COVID-19 vaccine doses†",
  mean_num_vac_mRNA = "Mean number of mRNA vaccine doses†",
  mean_num_vac_inact = "Mean number of inactivated-virus vaccine doses†",
  mean_num_vac_vector = "Mean number of viral vector vaccine doses†",
  mean_num_vac_other = "Mean number of other COVID-19 vaccine doses†",
  
  prop_last_vac_lt_6mo = "Percentage with last vaccination <6 months†",
  prop_last_vac_6_12mo = "Percentage with last vaccination 6–12 months†",
  prop_last_vac_gt_12mo = "Percentage with last vaccination >12 months†",
  prop_last_vac_never = "Percentage who had never received COVID-19 vaccination†",
  prop_last_vac_unknown = "Percentage with unknown COVID-19 vaccination history",
  
  # Sequencing
  prop_seq_avail = "Percentage of cases with available viral genome sequence data‡",
  major_lineage = "Household major SARS-CoV-2 lineage‡",
  
  # =============================
  # Children-specific measures
  # =============================
  prop_males_children = "Children: Percentage of males in the household",
  prop_females_children = "Children: Percentage of females in the household",
  male_to_female_ratio_children = "Children: Male-to-female ratio",
  
  mean_age_children = "Children: Mean age",
  mean_BMI_children = "Children: Mean BMI†",
  
  # Comorbidities
  prop_any_comorbidity_children = "Children: Percentage with any pre-existing medical condition†",
  prop_pre_existing_med_cond_hd_children = "Children: Percentage with heart disease†",
  prop_pre_existing_med_cond_ht_children = "Children: Percentage with hypertension†",
  prop_pre_existing_med_cond_ld_children = "Children: Percentage with lung disease†",
  prop_pre_existing_med_cond_lim_children = "Children: Percentage with low immunity†",
  prop_pre_existing_med_cond_kd_children = "Children: Percentage with kidney disease†",
  prop_pre_existing_med_cond_bd_children = "Children: Percentage with brain disorder†",
  prop_pre_existing_med_cond_gi_children = "Children: Percentage with gastrointestinal disease†",
  prop_pre_existing_med_cond_im_children = "Children: Percentage with immunodeficiency†",
  prop_pre_existing_med_cond_dm_children = "Children: Percentage with diabetes†",
  prop_pre_existing_med_cond_ar_children = "Children: Percentage with arthritis†",
  prop_pre_existing_med_cond_tb_children = "Children: Percentage with tuberculosis†",
  prop_pre_existing_med_cond_ml_children = "Children: Percentage with malaria†",
  prop_pre_existing_med_cond_al_children = "Children: Percentage with allergic rhinitis†",
  prop_pre_existing_med_cond_hl_children = "Children: Percentage with hyperlipidaemia†",
  prop_pre_existing_med_cond_as_children = "Children: Percentage with asthma†",
  
  # Infection
  prop_pos_enrol_children = "Children: Overall baseline case rate at enrolment",
  prop_pos_end_children = "Children: Overall endpoint case rate",
  prop_pos_add_children = "Children: Incident rate (crude attack rate)",
  
  # Symptom severity scores
  mean_most_severe_fv_children = "Children: Mean maximum fever severity score‡",
  mean_most_severe_co_children = "Children: Mean maximum cough severity score‡",
  mean_most_severe_st_children = "Children: Mean maximum sore throat severity score‡",
  mean_most_severe_rn_children = "Children: Mean maximum runny nose severity score‡",
  mean_most_severe_db_children = "Children: Mean maximum difficulty breathing severity score‡",
  mean_most_severe_he_children = "Children: Mean maximum headache severity score‡",
  mean_most_severe_my_children = "Children: Mean maximum myalgia severity score‡",
  mean_most_severe_ch_children = "Children: Mean maximum chills severity score‡",
  mean_most_severe_fa_children = "Children: Mean maximum fatigue severity score‡",
  mean_most_severe_di_children = "Children: Mean maximum diarrhoea severity score‡",
  mean_most_severe_lo_children = "Children: Mean maximum loss of smell/taste severity score‡",
  mean_most_severe_vo_children = "Children: Mean maximum vomiting severity score‡",
  mean_most_severe_ho_children = "Children: Mean maximum hoarseness severity score‡",
  
  max_most_severe_fv_children = "Children: Highest fever severity score‡",
  max_most_severe_co_children = "Children: Highest cough severity score‡",
  max_most_severe_st_children = "Children: Highest sore throat severity score‡",
  max_most_severe_rn_children = "Children: Highest runny nose severity score‡",
  max_most_severe_db_children = "Children: Highest difficulty breathing severity score‡",
  max_most_severe_he_children = "Children: Highest headache severity score‡",
  max_most_severe_my_children = "Children: Highest myalgia severity score‡",
  max_most_severe_ch_children = "Children: Highest chills severity score‡",
  max_most_severe_fa_children = "Children: Highest fatigue severity score‡",
  max_most_severe_di_children = "Children: Highest diarrhoea severity score‡",
  max_most_severe_lo_children = "Children: Highest loss of smell/taste severity score‡",
  max_most_severe_vo_children = "Children: Highest vomiting severity score‡",
  max_most_severe_ho_children = "Children: Highest hoarseness severity score‡",
  
  # Symptom groups
  prop_min_symp_cases_children = "Children: Percentage of minimally/(a)symptomatic cases‡",
  prop_mod_symp_cases_children = "Children: Percentage of moderately symptomatic cases‡",
  prop_highly_symp_cases_children = "Children: Percentage of highly symptomatic cases‡",
  prop_no_symp_data_cases_children = "Children: Percentage of cases without symptom data",
  
  # PCR Ct values
  mean_pcr_orf1ab_ct_children = "Children: Mean ORF1ab Ct value‡",
  mean_pcr_n_ct_children = "Children: Mean N gene Ct value‡",
  mean_pcr_e_ct_children = "Children: Mean E gene Ct value‡",
  
  min_pcr_orf1ab_ct_children = "Children: Lowest ORF1ab Ct value‡",
  min_pcr_n_ct_children = "Children: Lowest N gene Ct value‡",
  min_pcr_e_ct_children = "Children: Lowest E gene Ct value‡",
  
  # Previous COVID-19 history
  mean_num_prev_COVID_children = "Children: Mean number of previous COVID-19 infections†",
  prop_prev_COVID_lt_6mo_children = "Children: Percentage with last infection within 6 months†",
  prop_prev_COVID_6_12mo_children = "Children: Percentage with last infection 6–12 months earlier†",
  prop_prev_COVID_gt_12mo_children = "Children: Percentage with last infection more than 12 months earlier†",
  prop_prev_COVID_never_children = "Children: Percentage with no previous COVID-19 infection†",
  prop_prev_COVID_unknown_children = "Children: Percentage with unknown previous COVID-19 history",
  
  # Vaccination
  mean_num_vac_children = "Children: Mean number of COVID-19 vaccine doses†",
  mean_num_vac_mRNA_children = "Children: Mean number of mRNA vaccine doses†",
  mean_num_vac_inact_children = "Children: Mean number of inactivated-virus vaccine doses†",
  mean_num_vac_vector_children = "Children: Mean number of viral vector vaccine doses†",
  mean_num_vac_other_children = "Children: Mean number of other COVID-19 vaccine doses†",
  
  prop_last_vac_lt_6mo_children = "Children: Percentage with last vaccination <6 months†",
  prop_last_vac_6_12mo_children = "Children: Percentage with last vaccination 6–12 months†",
  prop_last_vac_gt_12mo_children = "Children: Percentage with last vaccination >12 months†",
  prop_last_vac_never_children = "Children: Percentage who had never received COVID-19 vaccination†",
  prop_last_vac_unknown_children = "Children: Percentage with unknown COVID-19 vaccination history",
  
  # Sequencing
  prop_seq_avail_children = "Children: Percentage of cases with available viral genome sequence data‡",
  
  # =============================
  # Adult-specific measures
  # =============================
  prop_males_adults = "Adults: Percentage of males in the household",
  prop_females_adults = "Adults: Percentage of females in the household",
  male_to_female_ratio_adults = "Adults: Male-to-female ratio",
  
  mean_age_adults = "Adults: Mean age",
  mean_BMI_adults = "Adults: Mean BMI†",
  
  # Comorbidities
  prop_any_comorbidity_adults = "Adults: Percentage with any pre-existing medical condition†",
  prop_pre_existing_med_cond_hd_adults = "Adults: Percentage with heart disease†",
  prop_pre_existing_med_cond_ht_adults = "Adults: Percentage with hypertension†",
  prop_pre_existing_med_cond_ld_adults = "Adults: Percentage with lung disease†",
  prop_pre_existing_med_cond_lim_adults = "Adults: Percentage with low immunity†",
  prop_pre_existing_med_cond_kd_adults = "Adults: Percentage with kidney disease†",
  prop_pre_existing_med_cond_bd_adults = "Adults: Percentage with brain disorder†",
  prop_pre_existing_med_cond_gi_adults = "Adults: Percentage with gastrointestinal disease†",
  prop_pre_existing_med_cond_im_adults = "Adults: Percentage with immunodeficiency†",
  prop_pre_existing_med_cond_dm_adults = "Adults: Percentage with diabetes†",
  prop_pre_existing_med_cond_ar_adults = "Adults: Percentage with arthritis†",
  prop_pre_existing_med_cond_tb_adults = "Adults: Percentage with tuberculosis†",
  prop_pre_existing_med_cond_ml_adults = "Adults: Percentage with malaria†",
  prop_pre_existing_med_cond_al_adults = "Adults: Percentage with allergic rhinitis†",
  prop_pre_existing_med_cond_hl_adults = "Adults: Percentage with hyperlipidaemia†",
  prop_pre_existing_med_cond_as_adults = "Adults: Percentage with asthma†",
  
  # Infection
  prop_pos_enrol_adults = "Adults: Overall baseline case rate at enrolment",
  prop_pos_end_adults = "Adults: Overall endpoint case rate",
  prop_pos_add_adults = "Adults: Incident rate (crude attack rate)",
  
  # Symptom severity scores
  mean_most_severe_fv_adults = "Adults: Mean maximum fever severity score‡",
  mean_most_severe_co_adults = "Adults: Mean maximum cough severity score‡",
  mean_most_severe_st_adults = "Adults: Mean maximum sore throat severity score‡",
  mean_most_severe_rn_adults = "Adults: Mean maximum runny nose severity score‡",
  mean_most_severe_db_adults = "Adults: Mean maximum difficulty breathing severity score‡",
  mean_most_severe_he_adults = "Adults: Mean maximum headache severity score‡",
  mean_most_severe_my_adults = "Adults: Mean maximum myalgia severity score‡",
  mean_most_severe_ch_adults = "Adults: Mean maximum chills severity score‡",
  mean_most_severe_fa_adults = "Adults: Mean maximum fatigue severity score‡",
  mean_most_severe_di_adults = "Adults: Mean maximum diarrhoea severity score‡",
  mean_most_severe_lo_adults = "Adults: Mean maximum loss of smell/taste severity score‡",
  mean_most_severe_vo_adults = "Adults: Mean maximum vomiting severity score‡",
  mean_most_severe_ho_adults = "Adults: Mean maximum hoarseness severity score‡",
  
  max_most_severe_fv_adults = "Adults: Highest fever severity score‡",
  max_most_severe_co_adults = "Adults: Highest cough severity score‡",
  max_most_severe_st_adults = "Adults: Highest sore throat severity score‡",
  max_most_severe_rn_adults = "Adults: Highest runny nose severity score‡",
  max_most_severe_db_adults = "Adults: Highest difficulty breathing severity score‡",
  max_most_severe_he_adults = "Adults: Highest headache severity score‡",
  max_most_severe_my_adults = "Adults: Highest myalgia severity score‡",
  max_most_severe_ch_adults = "Adults: Highest chills severity score‡",
  max_most_severe_fa_adults = "Adults: Highest fatigue severity score‡",
  max_most_severe_di_adults = "Adults: Highest diarrhoea severity score‡",
  max_most_severe_lo_adults = "Adults: Highest loss of smell/taste severity score‡",
  max_most_severe_vo_adults = "Adults: Highest vomiting severity score‡",
  max_most_severe_ho_adults = "Adults: Highest hoarseness severity score‡",
  
  # Symptom groups
  prop_min_symp_cases_adults = "Adults: Percentage of minimally/(a)symptomatic cases‡",
  prop_mod_symp_cases_adults = "Adults: Percentage of moderately symptomatic cases‡",
  prop_highly_symp_cases_adults = "Adults: Percentage of highly symptomatic cases‡",
  prop_no_symp_data_cases_adults = "Adults: Percentage of cases without symptom data",
  
  # PCR Ct values
  mean_pcr_orf1ab_ct_adults = "Adults: Mean ORF1ab Ct value‡",
  mean_pcr_n_ct_adults = "Adults: Mean N gene Ct value‡",
  mean_pcr_e_ct_adults = "Adults: Mean E gene Ct value‡",
  
  min_pcr_orf1ab_ct_adults = "Adults: Lowest ORF1ab Ct value‡",
  min_pcr_n_ct_adults = "Adults: Lowest N gene Ct value‡",
  min_pcr_e_ct_adults = "Adults: Lowest E gene Ct value‡",
  
  # Previous COVID-19 history
  mean_num_prev_COVID_adults = "Adults: Mean number of previous COVID-19 infections†",
  prop_prev_COVID_lt_6mo_adults = "Adults: Percentage with last infection within 6 months†",
  prop_prev_COVID_6_12mo_adults = "Adults: Percentage with last infection 6–12 months earlier†",
  prop_prev_COVID_gt_12mo_adults = "Adults: Percentage with last infection more than 12 months earlier†",
  prop_prev_COVID_never_adults = "Adults: Percentage with no previous COVID-19 infection†",
  prop_prev_COVID_unknown_adults = "Adults: Percentage with unknown previous COVID-19 history",
  
  # Vaccination
  mean_num_vac_adults = "Adults: Mean number of COVID-19 vaccine doses†",
  mean_num_vac_mRNA_adults = "Adults: Mean number of mRNA vaccine doses†",
  mean_num_vac_inact_adults = "Adults: Mean number of inactivated-virus vaccine doses†",
  mean_num_vac_vector_adults = "Adults: Mean number of viral vector vaccine doses†",
  mean_num_vac_other_adults = "Adults: Mean number of other COVID-19 vaccine doses†",
  
  prop_last_vac_lt_6mo_adults = "Adults: Percentage with last vaccination <6 months†",
  prop_last_vac_6_12mo_adults = "Adults: Percentage with last vaccination 6–12 months†",
  prop_last_vac_gt_12mo_adults = "Adults: Percentage with last vaccination >12 months†",
  prop_last_vac_never_adults = "Adults: Percentage who had never received COVID-19 vaccination†",
  prop_last_vac_unknown_adults = "Adults: Percentage with unknown COVID-19 vaccination history",
  
  # Sequencing
  prop_seq_avail_adults = "Adults: Percentage of cases with available viral genome sequence data‡"
)

# ###############################
# Create Table S5: Overview of the study cohort at the household level: 3 levels
# ###############################
# Omnibus test across all HH types
# -------------------------------
HH_sum_table_all_HHs <- make_HH_sum_table_func(HH_metadata)

sig_vars <- HH_sum_table_all_HHs$table_body %>%
  filter(!is.na(p.value), p.value < 0.05) %>%
  pull(variable)

# Pairwise post hoc tests for sig continuous variables using Wilcoxon
# -------------------------------
group_pairs <- list(
  c("Child", "Adult"),
  c("Child", "Adult & child"),
  c("Adult", "Adult & child")
)

cont_vars <- c(
  # -----------------------------
  # Overall household characteristics
  # -----------------------------
  "HH_size",
  "prop_children",
  "prop_adults",
  "child_to_adult_ratio",
  "prop_males",
  "prop_females",
  "male_to_female_ratio",
  "mean_age",
  "mean_BMI",
  
  # Comorbidities
  "prop_any_comorbidity",
  grep("^prop_pre_existing_med_cond_", names(HH_metadata), value = TRUE),

  # Infection outcomes
  "prop_pos_enrol",
  "prop_pos_end",
  "prop_pos_add",
  
  # Symptom severity scores
  grep("^mean_most_severe_", names(HH_metadata), value = TRUE),
  grep("^max_most_severe_", names(HH_metadata), value = TRUE),
  
  # Symptom profile groups
  "prop_min_symp_cases",
  "prop_mod_symp_cases",
  "prop_highly_symp_cases",
  "prop_no_symp_data_cases",
  
  # PCR Ct values
  "mean_pcr_orf1ab_ct",
  "mean_pcr_n_ct",
  "mean_pcr_e_ct",
  "min_pcr_orf1ab_ct",
  "min_pcr_n_ct",
  "min_pcr_e_ct",
  
  # Previous infection history
  "mean_num_prev_COVID",
  "prop_prev_COVID_lt_6mo",
  "prop_prev_COVID_6_12mo",
  "prop_prev_COVID_gt_12mo",
  "prop_prev_COVID_never",
  "prop_prev_COVID_unknown",
  
  # Vaccination
  "mean_num_vac",
  "mean_num_vac_mRNA",
  "mean_num_vac_inact",
  "mean_num_vac_vector",
  "mean_num_vac_other",
  "prop_last_vac_lt_6mo",
  "prop_last_vac_6_12mo",
  "prop_last_vac_gt_12mo",
  "prop_last_vac_never",
  "prop_last_vac_unknown",
  
  # Sequencing
  "prop_seq_avail",
  
  # -----------------------------
  # Child-specific measures
  # -----------------------------
  grep("_children$", names(HH_metadata), value = TRUE),
  
  # -----------------------------
  # Adult-specific measures
  # -----------------------------
  grep("_adults$", names(HH_metadata), value = TRUE)
) %>% unique()

cont_posthoc_results <- list()

for(v in sig_vars[sig_vars %in% cont_vars]) {
  var_results <- list()
  for(i in seq_along(group_pairs)) {
    pair <- group_pairs[[i]]
    tmp_dat <- HH_metadata %>% filter(first_HH_case_age_cat %in% pair)
    x <- tmp_dat %>%
      filter(first_HH_case_age_cat == pair[1]) %>%
      pull(all_of(v))
    y <- tmp_dat %>%
      filter(first_HH_case_age_cat == pair[2]) %>%
      pull(all_of(v))
    x <- x[!is.na(x)]
    y <- y[!is.na(y)]
    if(length(x) == 0 || length(y) == 0) {
      p_val <- NA
      statistic <- NA
    } else {
      test <- tryCatch(
        wilcox.test(
          x,
          y#, exact = FALSE
        ),
        error = function(e) NULL
      )
      if(is.null(test)) {
        p_val <- NA
        statistic <- NA
      } else {
        p_val <- test$p.value
        statistic <- unname(test$statistic)
      }
    }
    
    var_results[[i]] <- data.frame(
      .y. = v,
      group1 = pair[1],
      group2 = pair[2],
      n1 = length(x),
      n2 = length(y),
      statistic = statistic,
      p = p_val
    )
  }
  
  var_results <- bind_rows(var_results)
  var_results$p.adj <- p.adjust(var_results$p, method = "fdr")
  var_results$p.adj.signif <- rstatix::add_significance(var_results %>% select(p.adj))$p.adj.signif
  
  cont_posthoc_results[[v]] <- var_results
}

cont_posthoc_results <- bind_rows(cont_posthoc_results)

# Pairwise post hoc tests for sig categorical variables using Fisher tests
# -------------------------------
cat_vars <- c(
  "major_lineage"
)

cat_posthoc_results <- list()

for(v in sig_vars[sig_vars %in% cat_vars]) {
  var_results <- list()
  for(i in seq_along(group_pairs)) {
    pair <- group_pairs[[i]]
    tmp_dat <- HH_metadata %>%
      filter(first_HH_case_age_cat %in% pair) %>%
      select(first_HH_case_age_cat, all_of(v)) %>%
      filter(!is.na(.data[[v]])) %>%
      droplevels()
    
    n1 <- sum(tmp_dat$first_HH_case_age_cat == pair[1])
    n2 <- sum(tmp_dat$first_HH_case_age_cat == pair[2])
    
    if(n1 == 0 || n2 == 0) {
      p_val <- NA
    } else {
      tab <- table(tmp_dat[[v]], tmp_dat$first_HH_case_age_cat)
      if(nrow(tab) < 2 || ncol(tab) < 2) {
        p_val <- NA
      } else {
        fisher_res <- tryCatch(
          fisher.test(tab),
          error = function(e) NULL
        )
        if(is.null(fisher_res)) {
          p_val <- NA
        } else {
          p_val <- fisher_res$p.value
        }
      }
    }
    var_results[[i]] <- data.frame(
      .y. = v,
      group1 = pair[1],
      group2 = pair[2],
      n1 = n1,
      n2 = n2,
      statistic = NA,
      p = p_val
    )
  }
  var_results <- bind_rows(var_results)
  var_results$p.adj <- p.adjust(var_results$p, method = "fdr")
  var_results$p.adj.signif <- rstatix::add_significance(var_results %>% select(p.adj))$p.adj.signif
  
  cat_posthoc_results[[v]] <- var_results
}

cat_posthoc_results <- bind_rows(cat_posthoc_results)

# Combine all results
# -------------------------------
all_posthoc_results <- bind_rows(cont_posthoc_results, cat_posthoc_results) %>% rename("variable" = ".y.")

# For the sake of binding FDR-adjusted p values to the main table, 
# make placeholder gtsummary tables comparing ...
# -------------------------------
HH_sum_table_child_VS_adult_first_HHs <- make_HH_sum_table_func(HH_metadata %>% filter(first_HH_case_age_cat %in% c("Child", "Adult")) %>% droplevels()) # Child-first HHs VS Adult-first HHs ...
HH_sum_table_child_VS_adult_and_child_co_first_HHs <- make_HH_sum_table_func(HH_metadata %>% filter(first_HH_case_age_cat %in% c("Child", "Adult & child")) %>% droplevels()) # Child-first HHs VS Adult-and-child co-first HHs ...
HH_sum_table_adult_VS_adult_and_child_co_first_HHs <- make_HH_sum_table_func(HH_metadata %>% filter(first_HH_case_age_cat %in% c("Adult", "Adult & child")) %>% droplevels()) # Adult-first HHs VS Adult-and-child co-first HHs ...

# ... then replace p values with the FDR-adjusted ones
# -------------------------------
HH_sum_table_child_VS_adult_first_HHs <- replace_pvalues(
  HH_sum_table_child_VS_adult_first_HHs,
  all_posthoc_results %>% filter(group1 == "Child", group2 == "Adult")
)

HH_sum_table_child_VS_adult_and_child_co_first_HHs <- replace_pvalues(
  HH_sum_table_child_VS_adult_and_child_co_first_HHs,
  all_posthoc_results %>% filter(group1 == "Child", group2 == "Adult & child")
)

HH_sum_table_adult_VS_adult_and_child_co_first_HHs <- replace_pvalues(
  HH_sum_table_adult_VS_adult_and_child_co_first_HHs,
  all_posthoc_results %>% filter(group1 == "Adult", group2 == "Adult & child")
)

# Combine all tables together and save it to file
# -------------------------------
Table_S5 <- tbl_merge(
  tbls = list(
    HH_sum_table_all_HHs,
    HH_sum_table_child_VS_adult_first_HHs,
    HH_sum_table_child_VS_adult_and_child_co_first_HHs,
    HH_sum_table_adult_VS_adult_and_child_co_first_HHs
  ),
  tab_spanner = c(
    "Overall comparison",
    "Child-first VS adult-first households",
    "Child-first VS adult-and-child co-first households",
    "Adult-first VS adult-and-child co-first households"
  )
)

gt::gtsave(
  data = Table_S5 %>% as_gt %>%
    gt::tab_source_note(source_note = "† Comorbidity, previous infection, vaccination, and BMI variables were calculated among household members with available information.\n‡ Symptom severity scores, PCR Ct values, and lineage information were calculated among cases identified by the end of follow-up with available information."),
  filename = Table_S5_file
)

# ###############################
# Create Table 1: Overview of the study cohort at the household level: 2 levels--Child-first HHs VS Adult-first HHs
# ###############################
# Change "Adult & child" to "Adult"
# ----------------------------------
HH_metadata_2HHlv <- HH_metadata %>% mutate(first_HH_case_age_cat = factor(recode_factor(first_HH_case_age_cat, "Adult & child" = "Adult"), c("Child", "Adult")))

# Split the df into 3: Overall variables, Children variables, and Adult variables
# and rename the children/adult columns to match the overall naming scheme
# ----------------------------------
HH_metadata_2HHlv <- HH_metadata_2HHlv %>% 
  rename("prop_children_overall"="prop_children", "prop_adults_overall"="prop_adults") # We do this to so that these two overall variables will be kept in the overall table, and not child-specific measure and adult-specific measure table

HH_metadata_overall <- HH_metadata_2HHlv %>%
  select(-ends_with("_children"), -ends_with("_adults")) %>% 
  rename("prop_children"="prop_children_overall", "prop_adults"="prop_adults_overall") # Rename them back

HH_metadata_children <- HH_metadata_2HHlv %>%
  select(first_HH_case_age_cat, ends_with("_children")) %>%
  rename_with(~ gsub("_children$", "", .x))

HH_metadata_adults <- HH_metadata_2HHlv %>%
  select(first_HH_case_age_cat, ends_with("_adults")) %>%
  rename_with(~ gsub("_adults$", "", .x))

# Generate 3 tbl_summary tables using the same make_HH_sum_table_func()
# ----------------------------------
# Note: tbl_merge prioritises the structure of the first table.
# If the first table doesn't have a missing row for a particular variable, but subsequent tables have them, 
# missing rows of that particular variable of subsequent tables will be pushed to the bottom!
# So, we set missing = "always" so that the 3 tables will be merged correctly,
# and then manually remove uninformative rows afterward
tbl_overall <- make_HH_sum_table_func(HH_metadata_overall, missing = "always")
tbl_children <- make_HH_sum_table_func(HH_metadata_children, missing = "always")
tbl_adults <- make_HH_sum_table_func(HH_metadata_adults, missing = "always")

# Merge the tables horizontally
# ----------------------------------
Table_1 <- tbl_merge(
  tbls = list(
    tbl_overall,
    tbl_children,
    tbl_adults
  ),
  tab_spanner = c(
    "**Overall measures**",
    "**Child-specific measures**",
    "**Adult-specific measures**"
  )
)

# Remove uninformative missing rows
# ----------------------------------
Table_1$table_body <- Table_1$table_body %>%
  rowwise() %>%
  filter(!(row_type == "missing" && all(c_across(matches("^stat_[012]_[123]$")) %in% c("0", NA)))) %>%
  ungroup()

# Save it to file
# -------------------------------
gt::gtsave(
  data = Table_1 %>% as_gt %>% 
    gt::tab_source_note(source_note = "† Comorbidity, previous infection, vaccination, and BMI variables were calculated among household members with available information.\n‡ Symptom severity scores, PCR Ct values, and lineage information were calculated among cases identified by the end of follow-up with available information."),
  filename = Table_1_file
)