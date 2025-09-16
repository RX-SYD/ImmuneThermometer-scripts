############################################
# Forest plot for Cox regression results
# Author: RX2025
# ------------------------------------------
# Input:
#   cox_data.csv 
#   Required columns (example):
#     1. CancerType (or variable name)
#     2. SampleSize (n)
#     3. HR (hazard ratio)
#     4. Lower95 (lower CI)
#     5. Upper95 (upper CI)
# Output:
#   forestplot.pdf (or forestplot.png)
############################################

# 0) Working directory ---------------------------------------------------------
setwd("YOUR FOLDER")   # <-- replace with your actual path

# 1) Packages ------------------------------------------------------------------
if (!require("forestplot")) install.packages("forestplot", dependencies = TRUE)
library(forestplot)
if (!require("grid")) install.packages("grid")  # for unit()

# 2) Load input data -----------------------------------------------------------
sample <- read.csv("cox_data.csv", stringsAsFactors = FALSE)

# 3) Prepare table text --------------------------------------------------------
tabletext1 <- as.character(sample[, 1])  # variable names
tabletext2 <- round(as.numeric(sample[, 2]), 5)  # sample size
tabletext3 <- paste(round(as.numeric(sample[, 3]), 3),
                    round(as.numeric(sample[, 4]), 2), sep = "(")
tabletext4 <- paste(tabletext3, round(as.numeric(sample[, 5]), 2), sep = "-")
tabletext5 <- paste0(tabletext4, ")")
tabletext  <- cbind(tabletext1, tabletext2, tabletext5)

# 4) Draw forest plot ----------------------------------------------------------
pdf("forestplot.pdf", width = 8, height = 6)   # save as PDF
forestplot(
  labeltext = tabletext,  # text info
  mean      = round(sample[, 3], 3),  # HR
  lower     = round(sample[, 4], 2),  # CI lower
  upper     = round(sample[, 5], 2),  # CI upper
  boxsize   = 0.8,                    # box size
  graph.pos = 3,                      # graph column position
  graphwidth= unit(0.4, "npc"),       # graph width
  fn.ci_norm= "fpDrawDiamondCI",      # diamond shape
  col       = fpColors(box = "steelblue", lines = "black", zero = "black"),
  lwd.ci    = 1,
  ci.vertices.height = 0.1,
  ci.vertices = TRUE,                 # CI line caps
  zero      = 1,                      # reference line at HR=1
  lwd.zero  = 1,
  grid      = TRUE,
  lwd.xaxis = 1,
  title     = "Hazard Ratio",
  xlab      = "",                     # X-axis label
  clip      = c(-Inf, 3),             # axis limits
  colgap    = unit(0.5, "cm"),
  new_page  = TRUE
)
dev.off()

# 5) Optionally also save PNG --------------------------------------------------
# png("forestplot.png", width = 1200, height = 900, res = 150)
# <repeat forestplot() code here>
# dev.off()

############################################
# End of script
############################################
