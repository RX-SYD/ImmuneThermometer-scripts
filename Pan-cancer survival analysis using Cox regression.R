############################################
# Pan-cancer survival analysis using Cox regression
# Author: RX-2025
# ------------------------------------------
# Input:
#   pansurv.csv   (clinical survival data, rows = samples, cols = OS, OS.time, type)
#   panexpr.csv   (expression/signature data, rows = genes/signatures, cols = samples)
# Output:
#   hr.txt        (hazard ratios)
#   pvalue.txt    (p-values)
#   hrlow.txt     (95% CI lower bound)
#   hrhigh.txt    (95% CI upper bound)
############################################

# 0) Working directory ---------------------------------------------------------
setwd("YOUR FOLDER")   # <-- replace with your folder path

# 1) Install & load required packages ------------------------------------------
required_pkgs <- c("data.table", "survival", "pheatmap")
to_install <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)
lapply(required_pkgs, library, character.only = TRUE)

# 2) Load input data -----------------------------------------------------------
pansurv <- read.csv("pansurv.csv", row.names = 1, check.names = FALSE)
panexpr <- read.csv("panexpr.csv", row.names = 1, check.names = FALSE)

# 3) Match samples between expression and survival -----------------------------
comsam <- intersect(colnames(panexpr), rownames(pansurv))
pansurv <- pansurv[comsam, ]
panexpr <- panexpr[, comsam]

# 4) Define target signature/gene ----------------------------------------------
gene_group <- data.frame(Symbol = "temperature", stringsAsFactors = FALSE)
tumors <- unique(pansurv$type)

# 5) Initialize result matrices ------------------------------------------------
survland.coxhr   <- matrix(NA, nrow = nrow(gene_group), ncol = length(tumors),
                           dimnames = list(gene_group$Symbol, tumors))
survland.coxp    <- matrix(NA, nrow = nrow(gene_group), ncol = length(tumors),
                           dimnames = list(gene_group$Symbol, tumors))
survland.coxplot <- matrix(0,  nrow = nrow(gene_group), ncol = length(tumors),
                           dimnames = list(gene_group$Symbol, tumors))
survland.hrlow   <- matrix(0,  nrow = nrow(gene_group), ncol = length(tumors),
                           dimnames = list(gene_group$Symbol, tumors))
survland.hrhigh  <- matrix(0,  nrow = nrow(gene_group), ncol = length(tumors),
                           dimnames = list(gene_group$Symbol, tumors))

# 6) Cox regression per tumor type ---------------------------------------------
for (g in gene_group$Symbol) {
  for (t in tumors) {
    
    # Sample IDs for current tumor
    sam <- rownames(pansurv[pansurv$type == t, ])
    
    # Expression values for target signature
    expr <- as.numeric(panexpr[g, sam])
    
    # Prepare survival dataframe
    expr.surv <- data.frame(
      futime = pansurv[sam, "OS.time"],
      fustat = pansurv[sam, "OS"],
      expr   = expr,
      stringsAsFactors = FALSE
    )
    
    # Split into high/low groups by median
    expr.surv$group <- ifelse(expr > median(expr), "high", "low")
    expr.surv$group <- factor(expr.surv$group, levels = c("low", "high"))
    
    # Cox regression
    cox <- coxph(Surv(futime, fustat) ~ group, data = expr.surv)
    coxSummary <- summary(cox)
    
    # Extract HR, p-value, 95% CI
    hr     <- as.numeric(coxSummary$coefficients[, "exp(coef)"])[1]
    pvalue <- as.numeric(coxSummary$coefficients[, "Pr(>|z|)"])[1]
    hrlow  <- as.numeric(coxSummary$conf.int[, "lower .95"])[1]
    hrhigh <- as.numeric(coxSummary$conf.int[, "upper .95"])[1]
    
    # Save results
    survland.coxhr[g, t]   <- hr
    survland.coxp[g, t]    <- pvalue
    survland.hrlow[g, t]   <- hrlow
    survland.hrhigh[g, t]  <- hrhigh
    
    # Mark significance for plotting
    if (pvalue < 0.05) {
      survland.coxplot[g, t] <- ifelse(hr > 1, 1, -1)  # 1 = risk factor, -1 = protective
    }
  }
}

# 7) Save results --------------------------------------------------------------
write.table(survland.coxhr,   file = "hr.txt",      sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
write.table(survland.coxp,    file = "pvalue.txt",  sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
write.table(survland.hrlow,   file = "hrlow.txt",   sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
write.table(survland.hrhigh,  file = "hrhigh.txt",  sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

############################################
# End of script
############################################
