# -------------------------------------------------------------------------
# Script: Professional qPCR Analysis Pipeline (Livak Method)
# Author: Hossein Noorollahi (Ph.D. Candidate in Molecular Genetics)
# Date: December 2025
# Description: Automated analysis pipeline for qPCR data.
#              Includes: Data Parsing, Livak Calculation (2^-ddCt), 
#              Welch's t-test, 95% CI, and High-Impact Visualizations.
# License: MIT License
# -------------------------------------------------------------------------

# --- 1. Setup & Package Loading ---
message("--- 1. Initializing Environment ---")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, ggrepel, RColorBrewer, scales, gridExtra)
message("Libraries loaded successfully.")


# --- 2. Raw Data Input ---
# USER NOTE: Replace this block with your own CSV data if needed.
# Format: Sample, Target, Ct
message("\n--- 2. Loading Raw Data ---")

# Demo Data: Standardized for reproducibility
# Treat_A: Strong bimodal effect (Genes 1-5 UP, 6-10 DOWN)
# Treat_B: Specific targets only (Gene 2 UP, Gene 6 DOWN)
raw_data <- "Sample,Target,Ct
Control,GAPDH,18.1
Control,GAPDH,18.0
Control,GAPDH,18.2
Treat_A,GAPDH,18.1
Treat_A,GAPDH,18.2
Treat_A,GAPDH,18.0
Treat_B,GAPDH,18.3
Treat_B,GAPDH,18.1
Treat_B,GAPDH,18.2
Control,Gene_01,25.0
Control,Gene_01,25.2
Control,Gene_01,24.9
Treat_A,Gene_01,20.1
Treat_A,Gene_01,20.0
Treat_A,Gene_01,20.2
Treat_B,Gene_01,24.8
Treat_B,Gene_01,25.1
Treat_B,Gene_01,24.9
Control,Gene_02,26.0
Control,Gene_02,26.1
Control,Gene_02,25.9
Treat_A,Gene_02,21.0
Treat_A,Gene_02,21.2
Treat_A,Gene_02,20.9
Treat_B,Gene_02,21.5
Treat_B,Gene_02,21.4
Treat_B,Gene_02,21.6
Control,Gene_03,24.5
Control,Gene_03,24.6
Control,Gene_03,24.4
Treat_A,Gene_03,19.5
Treat_A,Gene_03,19.4
Treat_A,Gene_03,19.6
Treat_B,Gene_03,24.4
Treat_B,Gene_03,24.7
Treat_B,Gene_03,24.5
Control,Gene_04,25.5
Control,Gene_04,25.4
Control,Gene_04,25.6
Treat_A,Gene_04,20.5
Treat_A,Gene_04,20.6
Treat_A,Gene_04,20.4
Treat_B,Gene_04,25.3
Treat_B,Gene_04,25.5
Treat_B,Gene_04,25.4
Control,Gene_05,27.0
Control,Gene_05,27.1
Control,Gene_05,26.9
Treat_A,Gene_05,22.0
Treat_A,Gene_05,22.1
Treat_A,Gene_05,21.9
Treat_B,Gene_05,26.8
Treat_B,Gene_05,27.2
Treat_B,Gene_05,26.9
Control,Gene_06,22.0
Control,Gene_06,22.1
Control,Gene_06,21.9
Treat_A,Gene_06,27.0
Treat_A,Gene_06,27.2
Treat_A,Gene_06,26.9
Treat_B,Gene_06,26.5
Treat_B,Gene_06,26.4
Treat_B,Gene_06,26.6
Control,Gene_07,23.0
Control,Gene_07,23.2
Control,Gene_07,22.8
Treat_A,Gene_07,28.0
Treat_A,Gene_07,28.1
Treat_A,Gene_07,27.9
Treat_B,Gene_07,22.9
Treat_B,Gene_07,23.1
Treat_B,Gene_07,22.8
Control,Gene_08,21.5
Control,Gene_08,21.6
Control,Gene_08,21.4
Treat_A,Gene_08,26.5
Treat_A,Gene_08,26.6
Treat_A,Gene_08,26.4
Treat_B,Gene_08,21.4
Treat_B,Gene_08,21.7
Treat_B,Gene_08,21.5
Control,Gene_09,24.0
Control,Gene_09,24.1
Control,Gene_09,23.9
Treat_A,Gene_09,29.0
Treat_A,Gene_09,29.1
Treat_A,Gene_09,28.9
Treat_B,Gene_09,23.8
Treat_B,Gene_09,24.2
Treat_B,Gene_09,24.0
Control,Gene_10,25.0
Control,Gene_10,25.1
Control,Gene_10,24.9
Treat_A,Gene_10,30.0
Treat_A,Gene_10,30.2
Treat_A,Gene_10,29.9
Treat_B,Gene_10,25.1
Treat_B,Gene_10,24.9
Treat_B,Gene_10,25.0"

df_raw <- read_csv(raw_data, show_col_types = FALSE)
message("Data Loaded.")


# --- 3. Analysis Logic ---
message("\n--- 3. Running Analysis ---")

ref_gene <- "GAPDH"
control_group <- "Control"
treatments <- setdiff(unique(df_raw$Sample), control_group) 
target_genes <- setdiff(unique(df_raw$Target), ref_gene)

# Mean Ct Ref
df_ref_mean <- df_raw %>%
  filter(Target == ref_gene) %>%
  group_by(Sample) %>%
  summarise(Ref_Mean_Ct = mean(Ct), .groups = 'drop')

# Calculate dCt
df_analysis <- df_raw %>%
  filter(Target != ref_gene) %>%
  left_join(df_ref_mean, by = "Sample") %>%
  mutate(dCt = Ct - Ref_Mean_Ct)

results_list <- list()

for (trt in treatments) {
  for (g in target_genes) {
    dCt_ctrl <- df_analysis %>% filter(Sample == control_group, Target == g) %>% pull(dCt)
    dCt_trt  <- df_analysis %>% filter(Sample == trt, Target == g) %>% pull(dCt)
    
    if(length(dCt_ctrl) > 1 & length(dCt_trt) > 1) {
      t_test <- t.test(dCt_trt, dCt_ctrl, var.equal = FALSE)
      p_val <- t_test$p.value
      se_diff <- t_test$stderr
      df_val <- t_test$parameter
    } else {
      p_val <- NA; se_diff <- NA; df_val <- NA
    }
    
    mean_dCt_ctrl <- mean(dCt_ctrl)
    mean_dCt_trt <- mean(dCt_trt)
    ddCt <- mean_dCt_trt - mean_dCt_ctrl
    fold_change <- 2^(-ddCt)
    log2fc <- log2(fold_change)
    
    if (!is.na(df_val)) {
      t_crit <- qt(0.975, df = df_val)
      se_log2fc <- se_diff 
    } else {
      se_log2fc <- 0
    }
    
    sig_label <- "ns"
    if (!is.na(p_val)) {
      if (p_val < 0.001) sig_label <- "***"
      else if (p_val < 0.01) sig_label <- "**"
      else if (p_val < 0.05) sig_label <- "*"
    }
    
    reg_status <- "NS"
    if (!is.na(p_val) && p_val < 0.05 && abs(log2fc) > 1) {
      reg_status <- ifelse(log2fc > 0, "Upregulated", "Downregulated")
    }
    
    results_list[[length(results_list)+1]] <- data.frame(
      Treatment = trt,
      Target = g,
      Fold_Change = fold_change,
      Log2FC = log2fc,
      P_Value = p_val,
      SE = se_log2fc,
      Significance = sig_label,
      Regulation = reg_status
    )
  }
}

df_results <- bind_rows(results_list)
write_csv(df_results, "qPCR_Detailed_Analysis.csv")
message("Analysis Complete.")


# --- 5. Visualization ---
message("\n--- 5. Generating Visualizations ---")

# --- Define Professional Color Palette (Global) ---
# Treat_A = Teal/Greenish Blue (Scientific look)
# Treat_B = Burnt Orange/Red (High contrast)
my_colors <- c("Treat_A" = "#00AFBB", "Treat_B" = "#E7B800")

theme_set(theme_minimal(base_size = 14) +
            theme(
              axis.line = element_line(color = "black", linewidth = 0.6),
              panel.grid.major.x = element_blank(),
              panel.grid.minor = element_blank(),
              plot.title = element_text(face = "bold", hjust = 0.5, size=16),
              legend.position = "top",
              axis.text = element_text(color="black")
            ))

# --- A. Global Bar Plot (Improved Colors) ---
p_global <- ggplot(df_results, aes(x = Target, y = Log2FC, fill = Treatment)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7, color = "black", size=0.3) +
  geom_errorbar(aes(ymin = Log2FC - SE, ymax = Log2FC + SE), 
                position = position_dodge(0.8), width = 0.25, size=0.5) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  
  # Professional Manual Colors
  scale_fill_manual(values = my_colors) +
  
  labs(title = "Global Gene Expression Profile", 
       subtitle = "Relative to GAPDH (log2 Fold Change)",
       y = expression(log[2]("Fold Change")),
       x = "Target Gene") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Global_Expression_Profile.png", p_global, width = 12, height = 7, dpi = 300)

# --- B. Heatmap (Restored & Styled) ---
message("Generating Heatmap...")
p_heat <- ggplot(df_results, aes(x = Treatment, y = Target, fill = Log2FC)) +
  geom_tile(color = "white", size = 0.5) +
  # Red-White-Blue Gradient (Standard for expression)
  scale_fill_gradient2(low = "#313695", mid = "white", high = "#A50026", midpoint = 0) +
  geom_text(aes(label = Significance), color = "black", size = 5, vjust=0.7) +
  labs(title = "Expression Heatmap", 
       subtitle = "Color: log2FC | Label: Significance (*)",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(face="bold"))

ggsave("Heatmap_Expression.png", p_heat, width = 6, height = 8, dpi = 300)

# --- C. Volcano Plots ---
fc_cutoff <- 1
p_cutoff <- 0.05

for (trt in treatments) {
  data_volcano <- df_results %>% filter(Treatment == trt)
  
  p_vol <- ggplot(data_volcano, aes(x = Log2FC, y = -log10(P_Value), color = Regulation)) +
    geom_point(size = 4, alpha = 0.85) +
    # Consistent color scheme
    scale_color_manual(values = c("Upregulated" = "#A50026", "Downregulated" = "#313695", "NS" = "grey75")) +
    geom_text_repel(aes(label = ifelse(Regulation != "NS", as.character(Target), "")), 
                    box.padding = 0.6, max.overlaps = Inf, fontface="bold") +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", color="grey30") +
    geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color="grey30") +
    labs(title = paste("Volcano Plot:", trt),
         subtitle = "Sig: p < 0.05 & |log2FC| > 1",
         x = expression(log[2]("Fold Change")),
         y = expression(-log[10]("p-value")))
  
  ggsave(paste0("Volcano_", trt, ".png"), p_vol, width = 8, height = 6, dpi = 300)
}

# --- D. Individual Bar Plots ---
message("Generating Individual Plots...")
dir.create("Individual_Gene_Plots", showWarnings = FALSE)

for (gene in target_genes) {
  data_gene <- df_results %>% filter(Target == gene)
  max_y <- max(abs(data_gene$Log2FC)) + max(data_gene$SE) + 1
  
  p_ind <- ggplot(data_gene, aes(x = Treatment, y = Log2FC, fill = Treatment)) +
    geom_bar(stat = "identity", width = 0.6, color = "black", show.legend = FALSE) +
    geom_errorbar(aes(ymin = Log2FC - SE, ymax = Log2FC + SE), width = 0.2) +
    geom_hline(yintercept = 0, color = "black") +
    geom_text(aes(label = Significance, y = ifelse(Log2FC > 0, Log2FC + SE + 0.2, Log2FC - SE - 0.5)), 
              size = 6, fontface = "bold") +
    
    # Use same global palette
    scale_fill_manual(values = my_colors) +
    scale_y_continuous(limits = c(-max_y, max_y)) +
    labs(title = gene, subtitle = "Fold Change vs Control", y = "log2(Fold Change)", x = NULL) +
    theme_bw(base_size = 14) +
    theme(panel.grid.major.x = element_blank(), plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave(file.path("Individual_Gene_Plots", paste0("Plot_", gene, ".png")), p_ind, width = 4, height = 5, dpi = 300)
}

message("\n--- Full Pipeline Finished. All plots generated. ---")