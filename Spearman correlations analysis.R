############################################
# Immune infiltration vs. Immune Thermometer
# Scatter (with LM smooth) + Boxplots + Spearman correlations
# Author: RX-2025
# ------------------------------------------
# Input:
#   "Immune infiltration.csv"
#   Required columns:
#     - temperature          (numeric)
#     - group                (factor/character with 2+ groups)
#     - B.cells
#     - Activated.CD4.T.cell
#     - Activated.CD8.T.cell
#     - NK.cells
# Output:
#   figures/  (PNG files)
#   spearman_correlations.csv
############################################

# 0) Working directory ---------------------------------------------------------
setwd("YOUR FOLDER")  # <-- change to your folder path

# 1) Packages ------------------------------------------------------------------
required_pkgs <- c("ggplot2", "ggthemes", "scales")
to_install <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)
lapply(required_pkgs, library, character.only = TRUE)

# 2) Read data -----------------------------------------------------------------
# Keep R's default check.names=TRUE so "Activated.CD4.T.cell" style names are OK
Immune_infiltration_data <- read.csv("Immune infiltration.csv", stringsAsFactors = FALSE)

# Quick sanity check (optional)
# str(Immune_infiltration_data)
# head(Immune_infiltration_data)

# 3) Output folder -------------------------------------------------------------
if (!dir.exists("figures")) dir.create("figures")

# 4) A helper for color gradient ----------------------------------------------
grad2 <- scale_color_gradient2(
  limits    = c(-20, 20),
  low       = "#4F94CD",
  mid       = "grey",
  high      = "#8B2252",
  midpoint  = 0
)

# 5) Plot helpers --------------------------------------------------------------
plot_scatter_lm <- function(df, y_col, y_label) {
  p <- ggplot(df, aes(x = temperature, y = .data[[y_col]])) +
    geom_smooth(method = "lm", se = TRUE, color = "grey") +
    geom_point(aes(color = temperature), size = 2.5) +
    grad2 +
    labs(x = "Immune Thermometer (temperature)", y = y_label) +
    theme_bw()
  p
}

plot_box <- function(df, y_col, y_label) {
  p <- ggplot(df, aes(x = group, y = .data[[y_col]], color = group)) +
    stat_boxplot(geom = "errorbar", width = 0.4) +
    geom_boxplot(outlier.shape = NA) +
    scale_color_manual(values = c("#4F94CD", "#8B2252")) +
    labs(x = "Group", y = y_label, color = "Group") +
    theme_classic()
  p
}

# 6) Targets to analyze --------------------------------------------------------
cell_markers <- c(
  "B.cells"               = "B cells",
  "Activated.CD4.T.cell"  = "Activated CD4+ T cells",
  "Activated.CD8.T.cell"  = "Activated CD8+ T cells",
  "NK.cells"              = "NK cells"
)

# 7) Run correlations + save plots --------------------------------------------
spearman_out <- data.frame(
  Marker = character(),
  Spearman_rho = numeric(),
  P_value = numeric(),
  stringsAsFactors = FALSE
)

for (colname in names(cell_markers)) {
  y_label <- cell_markers[[colname]]
  
  # Scatter + LM smooth
  p_scatter <- plot_scatter_lm(Immune_infiltration_data, colname, y_label)
  ggsave(
    filename = file.path("figures", paste0("scatter_", colname, ".png")),
    plot = p_scatter, width = 5, height = 4, dpi = 300
  )
  
  # Boxplot by group
  p_box <- plot_box(Immune_infiltration_data, colname, y_label)
  ggsave(
    filename = file.path("figures", paste0("boxplot_", colname, ".png")),
    plot = p_box, width = 4.5, height = 4, dpi = 300
  )
  
  # Spearman correlation (exact=FALSE for large N)
  ct <- suppressWarnings(cor.test(
    Immune_infiltration_data$temperature,
    Immune_infiltration_data[[colname]],
    method = "spearman",
    exact = FALSE
  ))
  
  spearman_out <- rbind(
    spearman_out,
    data.frame(
      Marker = y_label,
      Spearman_rho = unname(ct$estimate),
      P_value = ct$p.value,
      stringsAsFactors = FALSE
    )
  )
}

# 8) Save correlation table ----------------------------------------------------
write.csv(spearman_out, "spearman_correlations.csv", row.names = FALSE)

# 9) Session info (optional, helpful for reproducibility) ---------------------
# sessionInfo()
