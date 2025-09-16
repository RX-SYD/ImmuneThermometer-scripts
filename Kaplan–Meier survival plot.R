############################################
# Kaplan–Meier survival plot (Group vs OS/PFS...)
# Author: RX2025
# ------------------------------------------
# Input:
#   surv_data.txt 
#   Required columns (header, in any order):
#     - Time   : numeric, survival time (e.g., months/Days)
#     - Status : 1 = event/death, 0 = censored
#     - Group  : factor with 2+ levels (e.g., "low", "high")
# Output:
#   km_plot.pdf, km_plot.png
############################################

# 0) Working directory ---------------------------------------------------------
setwd("YOUR FOLDER")  # <-- change to your path

# 1) Options -------------------------------------------------------------------
options(warn = -1)  # suppress warnings in batch runs (optional)

# 2) Packages (auto-install if missing) ----------------------------------------
req_pkgs <- c("ggplot2", "survival", "survminer")
to_install <- req_pkgs[!(req_pkgs %in% installed.packages()[, "Package"])]
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)
suppressPackageStartupMessages({
  library(ggplot2)
  library(survival)
  library(survminer)
})

# 3) Load data -----------------------------------------------------------------
# Expect a tab-delimited file with header. We’ll coerce column names and types.
dat <- read.table("surv_data.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Standardize column names if needed
colnames(dat) <- tolower(colnames(dat))
# Accept common variants
name_map <- list(
  time   = c("time", "os.time", "futime"),
  status = c("status", "os", "fustat", "event"),
  group  = c("group", "score_group", "risk_group")
)

get_col <- function(candidates, df_names) {
  hits <- intersect(candidates, df_names)
  if (length(hits) == 0) stop(paste("Missing required column among:", paste(candidates, collapse=", ")))
  hits[1]
}

time_col   <- get_col(name_map$time,   colnames(dat))
status_col <- get_col(name_map$status, colnames(dat))
group_col  <- get_col(name_map$group,  colnames(dat))

# Subset & coerce
data <- data.frame(
  Time   = as.numeric(dat[[time_col]]),
  Status = as.numeric(dat[[status_col]]),
  Group  = as.factor(dat[[group_col]])
)

# Basic checks
if (any(!data$Status %in% c(0, 1))) {
  stop("Status must be coded as 0 (censored) or 1 (event).")
}
if (nlevels(data$Group) < 2) stop("Group must have at least two levels.")

# Optional: set factor order to display low vs high consistently
if (all(c("low", "high") %in% tolower(levels(data$Group)))) {
  levels(data$Group) <- tolower(levels(data$Group))
  data$Group <- factor(data$Group, levels = c("low", "high"))
}

# 4) Fit KM model --------------------------------------------------------------
fit <- survfit(Surv(Time, Status == 1) ~ Group, data = data)

# 5) Plot settings -------------------------------------------------------------
palette_cols <- c("#1B9E77", "#D95F02")  # adjust if >2 groups

# 6) Build KM plot -------------------------------------------------------------
p <- ggsurvplot(
  fit,
  data = data,
  risk.table = FALSE,           # set TRUE if you want the risk table
  pval = TRUE,                  # log-rank test p-value
  conf.int = FALSE,             # CI ribbon off (set TRUE if needed)
  palette = palette_cols,
  xlab = "Time (months)",
  ggtheme = theme_classic(),
  ncensor.plot = FALSE,
  conf.int.style = "ribbon",
  surv.median.line = NULL,
)

# 7) Save outputs --------------------------------------------------------------
# PDF
ggsave("km_plot.pdf", plot = p$plot, width = 6, height = 5, device = cairo_pdf)
# PNG
ggsave("km_plot.png", plot = p$plot, width = 1600/300, height = 1300/300, dpi = 300)

# Optional: also save with risk table
# if (TRUE) {
#   ggsave("km_plot_with_risktable.pdf", plot = p, width = 7, height = 7, device = cairo_pdf)
# }

# 8) Session info (for reproducibility) ---------------------------------------
# sessionInfo()

############################################
# End of script
############################################
