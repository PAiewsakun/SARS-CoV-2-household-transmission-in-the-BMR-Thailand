# ###############################
# Load global variables
# ###############################
source(file.path("miscellaneous.R"))

# ###############################
# Load required libraries
# ###############################
library(readxl)

library(tidyverse)
library(purrr)

library(umap)
library(cluster)  # clusGap, silhouette
library(mclust)   # ARI
library(aricode)  # NMI

library(ggplot2)
library(patchwork)
library(ComplexHeatmap)

# ###############################
# Helper functions
# ###############################
# Change cluster labels according to mean severity
# -------------------------------
change_cluster_labels_according_to_cluster_mean_severity <- function(cluster_assignments, symp_severity_dat) {
  tmp <- symp_severity_dat %>%
    mutate(cluster = cluster_assignments) %>%
    group_by(cluster) %>%
    summarise(avg_severity = mean(c_across(most_severe_fv:most_severe_ho)), .groups = "drop") %>%
    arrange(avg_severity) %>%
    mutate(label = c("Minimally/(a)symptomatic", "Moderately symptomatic", "Highly symptomatic"))
  
  cluster_labels <- setNames(tmp$label, tmp$cluster)
  cluster_labels[as.character(cluster_assignments)] %>% as.vector() %>% factor(levels = c("Minimally/(a)symptomatic", "Moderately symptomatic", "Highly symptomatic"))
}

# ###############################
# Symptom data clustering analysis to classify symptom level
# ###############################
# Load data
# -------------------------------
symp_severity_dat <- read.delim(participant_metadata_filepath, header = TRUE, sep = "\t") %>%
  select(
    subject_id, sex, age_cat, overall_res_end, 
    most_severe_fv:most_severe_ho
  ) %>%
  drop_na(most_severe_fv:most_severe_ho)

# Clustering loop, while also computing clustering scores (gap, tot_withinss, and silhouette)
# -------------------------------
set.seed(100)

n_reps <- 50
max_cluster_num <- 10
max_UMAP_dim <- 4

n_samples <- nrow(symp_severity_dat)
cluster_assignments_mat <- array(NA, dim = c(n_samples, n_reps, max_UMAP_dim, max_cluster_num))
clustering_score <- NULL
clustering_robustness <- NULL

for (rep in 1:n_reps) {
  cat(sprintf("=== Replicate %d / %d ===\n", rep, n_reps))
  
  for (UMAP_dim in 2:max_UMAP_dim) {
    cat(sprintf(" UMAP dim = %d\n", UMAP_dim))
    
    UMAP_coords <- umap(symp_severity_dat %>% select(most_severe_fv:most_severe_ho), n_components = UMAP_dim)$layout
    gap_stat <- clusGap(UMAP_coords, FUN = kmeans, K.max = max_cluster_num, B = 100)
    
    for (cluster_num in 1:max_cluster_num) {
      kmeans_clus_res <- kmeans(UMAP_coords, centers = cluster_num, nstart = 25)
      
      cluster_assignments_mat[, rep, UMAP_dim, cluster_num] <- kmeans_clus_res$cluster
      
      clustering_score = rbind(
        clustering_score,
        tibble(
          replicate = rep,
          UMAP_dim = UMAP_dim,
          cluster_num = cluster_num,
          tot_withinss = kmeans_clus_res$tot.withinss,
          silhouette = ifelse(cluster_num > 1, mean(silhouette(kmeans_clus_res$cluster, dist(UMAP_coords))[, 3]), NA),
          gap_stat = gap_stat$Tab[cluster_num, "gap"]
        )
      )
    } # end cluster_num
  } # end UMAP_dim
} # end rep

clustering_score <- clustering_score %>% mutate(UMAP_dim = factor(UMAP_dim))

clustering_score_summary <- clustering_score %>%
  group_by(UMAP_dim, cluster_num) %>%
  summarise(
    mean_tot_withinss = mean(tot_withinss, na.rm = TRUE),
    sd_tot_withinss   = sd(tot_withinss, na.rm = TRUE),
    mean_silhouette   = mean(silhouette, na.rm = TRUE),
    sd_silhouette     = sd(silhouette, na.rm = TRUE),
    mean_gap          = mean(gap_stat, na.rm = TRUE),
    sd_gap            = sd(gap_stat, na.rm = TRUE),
    .groups = "drop"
  )

write_tsv(clustering_score, "out/clustering_score.txt")

# Compute pairwise ARI & NMI distributions across replicates
# -------------------------------
for (UMAP_dim in 2:max_UMAP_dim) {
  cat(sprintf(" UMAP dim = %d\n", UMAP_dim))
  for (cluster_num in 1:max_cluster_num) {
    cat(sprintf(" cluster_num = %d/%d\n", cluster_num, max_cluster_num))
    
    mat <- cluster_assignments_mat[, , UMAP_dim, cluster_num]
    
    if (ncol(mat) < 2 || all(is.na(mat))) next
    
    combs <- combn(ncol(mat), 2)
    
    for (cc in 1:ncol(combs)) {
      rep1 <- combs[1, cc]
      rep2 <- combs[2, cc]
      
      clustering_robustness <- rbind(
        clustering_robustness,
        tibble(
          UMAP_dim = UMAP_dim,
          cluster_num = cluster_num,
          rep1 = rep1,
          rep2 = rep2,
          ARI = if (cluster_num > 1) adjustedRandIndex(mat[, rep1], mat[, rep2]) else NA,
          NMI = NMI(mat[, rep1], mat[, rep2])
        )
      )
    } # end comb rep
  } # end cluster_num
} # end UMAP_dim

clustering_robustness <- clustering_robustness %>% mutate(UMAP_dim = factor(UMAP_dim))

clustering_robustness_summary <- clustering_robustness %>%
  group_by(UMAP_dim, cluster_num) %>%
  summarise(
    mean_ARI = mean(ARI, na.rm = TRUE),
    sd_ARI = sd(ARI, na.rm = TRUE),
    mean_NMI = mean(NMI, na.rm = TRUE),
    sd_NMI = sd(NMI, na.rm = TRUE),
    .groups = "drop"
  )

write_tsv(clustering_robustness, "out/clustering_robustness.txt")

# ###############################
# Plot clustering scores, robustness scores, and severity heatmap
# ###############################
# Total within-cluster sum of squares
# -------------------------------
elbow_plot <- ggplot(data = clustering_score, aes(x = cluster_num, y = tot_withinss, color = UMAP_dim, group = interaction(UMAP_dim, replicate))) +
  geom_line(alpha = 0.25) + geom_point(alpha = 0.25) +
  geom_ribbon(data = clustering_score_summary, aes(x = cluster_num, ymin = mean_tot_withinss - sd_tot_withinss,
                                                   ymax = mean_tot_withinss + sd_tot_withinss, fill = UMAP_dim),
              alpha = 0.15, inherit.aes = FALSE) +
  
  geom_line(data = clustering_score_summary, aes(x = cluster_num, y = mean_tot_withinss, group = UMAP_dim), linewidth = 0.75*1.5, color = "black") +
  geom_line(data = clustering_score_summary, aes(x = cluster_num, y = mean_tot_withinss, color = UMAP_dim, group = UMAP_dim), linewidth = 0.75) +
  
  geom_vline(xintercept = 3, linetype = "dashed") +
  
  scale_x_continuous(name = "Number of clusters", breaks = seq(1, 10, 1), limits = c(0.4, 10.6), expand = c(0,0)) + 
  scale_color_discrete(name = "UMAP embedding\ndimensions", labels = c("2D", "3D", "4D")) +
  scale_fill_discrete(name = "UMAP embedding\ndimensions", labels = c("2D", "3D", "4D")) +
  labs(y = "Total within-cluster\nsum of squares") +
  
  theme_bw() + my_theme +
  theme(legend.title.position = "top", legend.direction = "horizontal")

# Silhouette
# -------------------------------
silhouette_plot <- ggplot(data = clustering_score, aes(x = cluster_num, y = silhouette, color = UMAP_dim, group = interaction(UMAP_dim, replicate))) +
  geom_line(alpha = 0.25) + geom_point(alpha = 0.25) +
  geom_ribbon(data = clustering_score_summary, aes(x = cluster_num, ymin = mean_silhouette - sd_silhouette,
                                                   ymax = mean_silhouette + sd_silhouette, fill = UMAP_dim),
              alpha = 0.15, inherit.aes = FALSE) +
  
  geom_line(data = clustering_score_summary, aes(x = cluster_num, y = mean_silhouette, group = UMAP_dim), linewidth = 0.75*1.5, color = "black") +
  geom_line(data = clustering_score_summary, aes(x = cluster_num, y = mean_silhouette, color = UMAP_dim, group = UMAP_dim), linewidth = 0.75) +
  
  geom_vline(xintercept = 3, linetype = "dashed") +
  
  scale_x_continuous(name = "Number of clusters", breaks = seq(1, 10, 1), limits = c(0.4, 10.6), expand = c(0,0)) + 
  scale_color_discrete(name = "UMAP embedding\ndimensions", labels = c("2D", "3D", "4D")) +
  scale_fill_discrete(name = "UMAP embedding\ndimensions", labels = c("2D", "3D", "4D")) +
  labs(y = "Mean silhouette width") +
  
  theme_bw() + my_theme +
  theme(legend.position = "none")

# Gap
# -------------------------------
gap_plot <- ggplot(data = clustering_score, aes(x = cluster_num, y = gap_stat, color = UMAP_dim, group = interaction(UMAP_dim, replicate))) +
  geom_line(alpha = 0.25) + geom_point(alpha = 0.25) +
  geom_ribbon(data = clustering_score_summary, aes(x = cluster_num, ymin = mean_gap - sd_gap,
                                                   ymax = mean_gap + sd_gap, fill = UMAP_dim),
              alpha = 0.15, inherit.aes = FALSE) +
  
  geom_line(data = clustering_score_summary, aes(x = cluster_num, y = mean_gap, group = UMAP_dim), linewidth = 0.75*1.5, color = "black") +
  geom_line(data = clustering_score_summary, aes(x = cluster_num, y = mean_gap, color = UMAP_dim, group = UMAP_dim), linewidth = 0.75) +
  
  geom_vline(xintercept = 3, linetype = "dashed") +
  
  scale_x_continuous(name = "Number of clusters", breaks = seq(1, 10, 1), limits = c(0.4, 10.6), expand = c(0,0)) + 
  scale_color_discrete(name = "UMAP embedding\ndimensions", labels = c("2D", "3D", "4D")) +
  scale_fill_discrete(name = "UMAP embedding\ndimensions", labels = c("2D", "3D", "4D")) +
  labs(y = "Gap statistic") +
  
  theme_bw() + my_theme +
  theme(legend.position = "none")

# Adjusted Rand index
# -------------------------------
ARI_plot <- ggplot() +
  geom_line(data = clustering_robustness_summary, aes(x = factor(cluster_num), y = mean_ARI, group = UMAP_dim), color = "black", linewidth = 0.75*1.5) +
  geom_line(data = clustering_robustness_summary, aes(x = factor(cluster_num), y = mean_ARI, group = UMAP_dim, color = UMAP_dim), linewidth = 0.75) +
  
  geom_boxplot(data = clustering_robustness, aes(x = factor(cluster_num), y = ARI, fill = UMAP_dim), alpha = 0.6, outlier.size = 1, position = position_dodge(width = 0.8)) +
  
  geom_vline(xintercept = 3, linetype = "dashed") +
  
  scale_color_discrete(name = "UMAP embedding\ndimensions", labels = c("2D", "3D", "4D")) +
  scale_fill_discrete(name = "UMAP embedding\ndimensions", labels = c("2D", "3D", "4D")) +
  labs(x = "Number of clusters", y = "Adjusted\nrand index") +
  
  theme_bw() + my_theme +
  theme(legend.position = "none")

# Normalized mutual information
# -------------------------------
NMI_plot <- ggplot() +
  geom_line(data = clustering_robustness_summary, aes(x = factor(cluster_num), y = mean_NMI, group = UMAP_dim), color = "black", linewidth = 0.75*1.5) +
  geom_line(data = clustering_robustness_summary, aes(x = factor(cluster_num), y = mean_NMI, group = UMAP_dim, color = UMAP_dim), linewidth = 0.75) +
  
  geom_boxplot(data = clustering_robustness, aes(x = factor(cluster_num), y = NMI, fill = UMAP_dim), alpha = 0.6, outlier.size = 1, position = position_dodge(width = 0.8)) +
  
  geom_vline(xintercept = 3, linetype = "dashed") +
  
  scale_color_discrete(name = "UMAP embedding\ndimensions", labels = c("2D", "3D", "4D")) +
  scale_fill_discrete(name = "UMAP embedding\ndimensions", labels = c("2D", "3D", "4D")) +
  labs(x = "Number of clusters", y = "Normalised\nmutual information") +
  
  theme_bw() + my_theme +
  theme(legend.position = "none")

# Symptom severity heatmap
# -------------------------------
# Re-perform the clustering using the selected best settings
set.seed(100)

best_UMAP_dim <- 3
best_cluster_num <- 3

UMAP_coords <- umap(symp_severity_dat %>% select(most_severe_fv:most_severe_ho), n_components = best_UMAP_dim)$layout
cluster_assignments <- kmeans(UMAP_coords, centers = best_cluster_num)
symp_cluster_res <- factor(change_cluster_labels_according_to_cluster_mean_severity(cluster_assignments$cluster, symp_severity_dat))

# Order rows by severity (min > mod > high), then by overall_res_end
symp_severity_dat <- symp_severity_dat %>%
  mutate(
    symp_cluster_res = factor(change_cluster_labels_according_to_cluster_mean_severity(cluster_assignments$cluster, symp_severity_dat))
  ) %>%
  arrange(symp_cluster_res, overall_res_end)

# Make heatmap matrix
symptom_severity_heatmap <- symp_severity_dat %>%
  column_to_rownames(var = "subject_id") %>% 
  select(most_severe_fv:most_severe_ho) %>%
  as.matrix()

# Plot heatmap
symptom_severity_heatmap_plot <- Heatmap(
  column_title = "Symptom level heatmap",
  column_title_gp = gpar(fontsize = 6, fontface = "bold"),
  
  symptom_severity_heatmap,
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_labels = c("Fever", "Cough", "Sore throat", "Runny nose", 
                    "Difficult\nbreathing", "Headache", "Myalgia", 
                    "Chills", "Fatigue", "Diarrhoea", 
                    "Loss of\nsmell/teste", "Vomit", "Hoarseness"),
  column_names_gp = gpar(fontsize = 6),
  
  row_split = data.frame(
    Severity = symp_severity_dat$symp_cluster_res#, "COVID-19\nstatus" = symp_severity_dat$overall_res_end # 2-level split
  ),
  border = TRUE,
  row_title_rot = 90,
  row_title_gp = gpar(fontsize = 6),
  
  left_annotation = rowAnnotation(
    "Symptom lv" = symp_severity_dat$symp_cluster_res,
    "COVID-19\nstatus" = symp_severity_dat$overall_res_end,
    annotation_name_gp = gpar(fontsize = 6, fontface = "bold"),
    simple_anno_size = unit(3, "mm"),
    
    annotation_legend_param = list(
      "Symptom lv" = list(
        title = "Symptom\nlevel",
        at = c("Minimally/(a)symptomatic", "Moderately symptomatic", "Highly symptomatic"),
        labels = c("Min", "Mod", "High"),
        title_gp = gpar(fontsize = 6, fontface = "bold"),
        labels_gp = gpar(fontsize = 6),
        grid_width = unit(2.5, "mm"), grid_height = unit(2.5, "mm")
      ),
      
      "COVID-19\nstatus" = list(
        title = "Endpoint\nCOVID-19\nstatus",
        at = c("Neg", "Pos"),
        labels = c("-ve", "+ve"),
        title_gp = gpar(fontsize = 6, fontface = "bold"),
        labels_gp = gpar(fontsize = 6),
        grid_width = unit(2.5, "mm"), grid_height = unit(2.5, "mm")
      )
    ),
    
    col = list(
      "Symptom lv" = severity_cols,
      "COVID-19\nstatus" = overall_res_end_cols
    )
  ),
  
  heatmap_legend_param = list(
    title = "Symptom\nscore",
    at = c(0, 1, 2, 3, 4),
    labels = c("0", "1", "2", "3", "4"),
    color_bar = "discrete",
    title_gp = gpar(fontsize = 6, fontface = "bold"),
    labels_gp = gpar(fontsize = 6),
    grid_width = unit(2.5, "mm"), grid_height = unit(2.5, "mm")
  ),
  col = circlize::colorRamp2(breaks = c(0, 2, 4), colors = c("white", "orange", "red"))
  
)

# Combine all the plots together into Figure S1
Figure_S1 <- cowplot::plot_grid(
  (elbow_plot / silhouette_plot / gap_plot / ARI_plot / NMI_plot ) + plot_layout(axis_titles = "collect"),
  grid.grabExpr(draw(symptom_severity_heatmap_plot)),
  nrow = 1,
  rel_widths = c(1, 1)
)

# Save Figure S1 to file
ggsave(filename = Figure_S1_file_png, device = "png",
       plot = Figure_S1, 
       width = 16, height = 20, units = "cm", dpi = 300, bg = "white")

ggsave(filename = Figure_S1_file_svg, device = "svg",
       plot = Figure_S1, 
       width = 16, height = 20, units = "cm", dpi = 300, bg = "white")

# Save clustering_eval_summary to file
clustering_eval_summary <- left_join(
  clustering_score_summary,
  clustering_robustness_summary,
  by = c("UMAP_dim" = "UMAP_dim", "cluster_num" = "cluster_num")
)

write_tsv(clustering_eval_summary, Table_S4_file)
