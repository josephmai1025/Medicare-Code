source("Helpers.R")

ensure_package("haven")
ensure_package("ggplot2")
ensure_package("rdrobust")
ensure_package("dplyr")
ensure_package("rdlocrand")
ensure_package("fixest")

data_path <- find_data_file()
data <- read_dta(data_path)
data <- subset(data, astatflg == 1 & age >= 0 & age <= 100)
data <- data[!is.na(data$insured) & !is.na(data$age) & !is.na(data$sampweight), ]

#### BANDWIDTH SETUP ####
bw_lookup <- c(
  insured        = 5.47,
  employed       = 2.93,
  health_r       = 7.20,
  flu_shot       = 5.68,
  presc_med      = 4.72,
  wellvisit_1yr  = 5.54,
  cholcheck_1yr  = 7.28,
  bpcheck_1yr    = 5.54
)
get_bw <- function(var_name) bw_lookup[[var_name]]

#### Robustness Check #1: Bandwidths ####

#### OUTCOME DEFINITIONS ####
outcome_vars <- c("health_r", "flu_shot", "presc_med",
                  "wellvisit_1yr", "cholcheck_1yr", "bpcheck_1yr")
outcome_labels <- c(
  health_r      = "Self-Reported Health",
  flu_shot      = "Flu Shot",
  presc_med     = "Prescription Medication",
  wellvisit_1yr = "Well Visit",
  cholcheck_1yr = "Cholesterol Check",
  bpcheck_1yr   = "Blood Pressure Check"
)

bandwidths <- 5:10

#### CORE FUZZY RD FUNCTION ####
run_fuzzy_rd_bw <- function(df, y_var_name, h) {
  # Row-wise complete cases only for the variables actually used in this call
  sub <- df[!is.na(df[[y_var_name]]), ]
  
  rd <- rdrobust(
    y = sub[[y_var_name]],
    x = sub$age,
    c = 65,
    p = 1,
    h = h,
    fuzzy = sub$insured,
    masspoints = "adjust",
    weights = sub$sampweight
  )
  
  est <- rd$coef[1]
  se  <- rd$se[3]
  z <- est / se
  p_val <- 2 * (1 - pnorm(abs(z)))
  stars <- ifelse(p_val < 0.01, "***",
                  ifelse(p_val < 0.05, "**",
                         ifelse(p_val < 0.1, "*", "")))
  
  data.frame(
    Outcome   = y_var_name,
    Bandwidth = h,
    Estimate  = round(est, 3),
    SE        = round(se, 3),
    Stars     = stars,
    stringsAsFactors = FALSE
  )
}

#### FULL GRID ####
late_results <- list()
for (h in bandwidths) {
  for (v in outcome_vars) {
    res <- tryCatch(
      run_fuzzy_rd_bw(data, v, h),
      error = function(e) {
        message("Failed: ", v, " at h=", h, " -- ", conditionMessage(e))
        data.frame(Outcome = v, Bandwidth = h, Estimate = NA, SE = NA, Stars = "",
                   stringsAsFactors = FALSE)
      }
    )
    late_results[[length(late_results) + 1]] <- res
  }
}
late_df <- do.call(rbind, late_results)
late_df$Outcome <- outcome_labels[late_df$Outcome]

print(late_df)

#### BUILD LATEX TABLE ####
late_df <- late_df %>%
  mutate(
    est = paste0(Estimate, Stars),
    se  = paste0("(", SE, ")")
  )

order_vec <- unname(outcome_labels)

latex_rows <- ""
for (outcome in order_vec) {
  sub <- late_df %>% filter(Outcome == outcome) %>% arrange(Bandwidth)
  est_row <- paste(sub$est, collapse = " & ")
  se_row  <- paste(sub$se, collapse = " & ")
  latex_rows <- paste0(
    latex_rows,
    outcome, " & ", est_row, " \\\\\n",
    " & ", se_row, " \\\\\n\n"
  )
}

latex_table <- paste0(
  "\\begin{table}[htbp]\n",
  "\\centering\n",
  "\\caption{Local Average Treatment Effects at Age 65 Across Bandwidths}\n",
  "\\label{tab:late_bw}\n",
  "\\begin{threeparttable}\n",
  "\\small\n",
  "\\begin{tabular}{lcccccc}\n",
  "\\hline\\hline\n",
  " & (5) & (6) & (7) & (8) & (9) & (10) \\\\\n",
  "\\hline\n",
  latex_rows,
  "\\hline\n",
  "\\end{tabular}\n",
  "\\begin{tablenotes}\n",
  "\\footnotesize\n",
  "\\item \\textit{Notes:} Each column reports fuzzy RD estimates (LATE) at age 65 using the bandwidth indicated. Medicare coverage is instrumented by eligibility at age 65. Robust standard errors are reported in parentheses. *** $p<0.01$, ** $p<0.05$, * $p<0.1$.\n",
  "\\end{tablenotes}\n",
  "\\end{threeparttable}\n",
  "\\end{table}\n"
)
cat(latex_table)

#### Robustness Check #2: Placebos ####

run_placebo <- function(varname, cutoff, h = 10) {
  h <- get_bw(varname)
  rd <- rdrobust(
    y = data[[varname]],
    x = data$age,
    c = cutoff,
    h = h,
    masspoints = "adjust",
    weights = data$sampweight
  )
  
  est <- as.numeric(rd$coef[2])
  se  <- as.numeric(rd$se[3])
  z <- est / se
  p <- 2 * (1 - pnorm(abs(z)))
  
  return(list(
    est = est,
    se = se,
    p = p
  ))
}

fake_cutoffs <- c(44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54)

outcomes <- list(
  insured = "Insurance",
  health_r = "Self-Reported Health",
  flu_shot = "Flu Shot",
  presc_med = "Prescription Meds",
  wellvisit_1yr = "Well Visit",
  cholcheck_1yr = "Cholesterol Check",
  bpcheck_1yr = "BP Check"
)

cutoff_vec <- fake_cutoffs

row_labels <- paste0(cutoff_vec)

est_matrix <- sapply(names(outcomes), function(var) {
  sapply(fake_cutoffs, function(cut) {
    res <- run_placebo(var, cut)
    
    stars <- ifelse(res$p < 0.01, "***",
                    ifelse(res$p < 0.05, "**",
                           ifelse(res$p < 0.1, "*", "")))
    
    paste0(round(res$est, 3), stars)
  })
})

se_matrix <- sapply(names(outcomes), function(var) {
  sapply(fake_cutoffs, function(cut) {
    res <- run_placebo(var, cut)
    paste0("(", round(res$se, 3), ")")
  })
})

est_rows <- apply(est_matrix, 1, function(x) {
  paste(x, collapse = " & ")
})

se_rows <- apply(se_matrix, 1, function(x) {
  paste(x, collapse = " & ")
})

body <- paste(
  sapply(1:length(cutoff_vec), function(i) {
    paste0(
      row_labels[i], " & ", est_rows[i], " \\\\\n",
      " & ", se_rows[i], " \\\\\n"
    )
  }),
  collapse = ""
)

col_names <- paste(unlist(outcomes), collapse = " & ")

latex_table <- paste0(
  "\\begin{table}[htbp]\n",
  "\\centering\n",
  "\\caption{Placebo Tests at Alternative Cutoffs}\n",
  "\\label{tab:placebo}\n",
  "\\begin{threeparttable}\n",
  "\\small\n",
  "\\resizebox{\\textwidth}{!}{%\n",
  "\\begin{tabular}{l", paste(rep("c", length(outcomes)), collapse = ""), "}\n",
  "\\hline\\hline\n",
  "Cutoff & ", col_names, " \\\\\n",
  "\\hline\n",
  body,
  "\\hline\n",
  "\\end{tabular}%\n",
  "}\n",
  "\\begin{tablenotes}\n",
  "\\small\n",
  "\\item Notes: Each entry reports bias-corrected RD estimates at placebo cutoffs. Robust standard errors in parentheses. ",
  "Significance levels: *** p<0.01, ** p<0.05, * p<0.1.\n",
  "\\end{tablenotes}\n",
  "\\end{threeparttable}\n",
  "\\end{table}\n"
)

cat(latex_table)


#### Robustness Check #3: Local Randomization ####


ri_health <- rdrandinf(
  Y = data$health_r,
  R = data$age,
  cutoff = 65,
  wl = 64,
  wr = 65,
  reps = 1000
)

print(ri_health)

ri_insured <- rdrandinf(
  Y = data$insured,
  R = data$age,
  cutoff = 65,
  wl = 64,
  wr = 65,
  reps = 1000
)

print(ri_insured)

ri_flushot <- rdrandinf(
  Y = data$flu_shot,
  R = data$age,
  cutoff = 65,
  wl = 64,
  wr = 65,
  reps = 1000
)

print(ri_flushot)

ri_employed <- rdrandinf(
  Y = data$worked12m,
  R = data$age,
  cutoff = 65,
  wl = 64,
  wr = 65,
  reps = 1000
)

print(ri_employed)


ri_meds <- rdrandinf(
  Y = data$presc_med,
  R = data$age,
  cutoff = 65,
  wl = 64,
  wr = 65,
  reps = 1000
)

print(ri_meds)

ri_wellvisit <- rdrandinf(
  Y = data$wellvisit_1yr,
  R = data$age,
  cutoff = 65,
  wl = 64,
  wr = 65,
  reps = 1000
)

print(ri_wellvisit)

ri_cholcheck <- rdrandinf(
  Y = data$cholcheck_1yr,
  R = data$age,
  cutoff = 65,
  wl = 64,
  wr = 65,
  reps = 1000
)

print(ri_cholcheck)

ri_bpcheck <- rdrandinf(
  Y = data$bpcheck_1yr,
  R = data$age,
  cutoff = 65,
  wl = 64,
  wr = 65,
  reps = 1000
)

print(ri_bpcheck)

ri_aggregate <- list(
  ri_insured, ri_employed, ri_health, ri_flushot,
  ri_meds, ri_wellvisit, ri_cholcheck, ri_bpcheck
)

ri_table <- data.frame(
  Outcome = c(
    "Insured",
    "Employed",
    "Self-Reported Health",
    "Flu Shot",
    "Prescription Meds",
    "Well Visit",
    "Cholesterol Check",
    "BP Check"
  ),
  Estimate = sapply(ri_aggregate, function(x) x[["obs.stat"]]),
  P_value = sapply(ri_aggregate, function(x) x[["p.value"]])
)

ri_table$stars <- ifelse(ri_table$P_value < 0.01, "***",
                         ifelse(ri_table$P_value < 0.05, "**",
                                ifelse(ri_table$P_value < 0.1, "*", "")))

ri_table$Estimate <- paste0(round(ri_table$Estimate, 3), ri_table$stars)

est_lookup <- setNames(ri_table$Estimate, ri_table$Outcome)

order_vec <- ri_table$Outcome

est_row <- paste(est_lookup[order_vec], collapse = " & ")
col_names <- paste(order_vec, collapse = " & ")
col_nums <- paste(paste0("(", seq_along(order_vec), ")"), collapse = " & ")

latex_table <- paste0(
  "\\begin{table}[htbp]\n",
  "\\centering\n",
  "\\caption{Local Randomization Estimates}\n",
  "\\label{tab:local_rand}\n",
  "\\begin{threeparttable}\n",
  "\\small\n",
  "\\resizebox{\\textwidth}{!}{%\n",
  "\\begin{tabular}{l", paste(rep("c", length(order_vec)), collapse = ""), "}\n",
  "\\hline\\hline\n",
  " & ", col_names, " \\\\\n",
  " & ", col_nums, " \\\\\n",
  "\\hline\n",
  "Difference in Means & ", est_row, " \\\\\n",
  "\\hline\n",
  "\\end{tabular}%\n",
  "}\n",
  "\n",
  "\\begin{tablenotes}[flushleft]\n",
  "\\footnotesize\n",
  "\\item \\textit{Notes:} Each column reports the difference in means across the cutoff using the local randomization approach within the window [64, 65]. *** $p<0.01$, ** $p<0.05$, * $p<0.1$.\n",
  "\\end{tablenotes}\n",
  "\\end{threeparttable}\n",
  "\\end{table}\n"
)

cat(latex_table)

### Robustness Check #4: Manual Regression ####


data$above_65 <- as.numeric(data$age >= 65)
data$age_c    <- data$age - 65

run_ols <- function(y_var, var_name) {
  h <- get_bw(var_name)
  data_bw       <- data[abs(data$age_c) <= h, ]
  fml <- as.formula(paste(y_var, "~ above_65 + age_c + age_c:above_65"))
  
  fit <- feols(fml,
               data    = data_bw,
               weights = ~sampweight,
               vcov = "hetero")
  
  est   <- coef(fit)["above_65"]
  se    <- se(fit)["above_65"]
  pval  <- pvalue(fit)["above_65"]
  stars <- ifelse(pval < 0.01, "***",
                  ifelse(pval < 0.05, "**",
                         ifelse(pval < 0.1,  "*", "")))
  
  data.frame(
    Outcome  = var_name,
    Estimate = paste0(round(est, 3), stars),
    SE       = paste0("(", round(se, 3), ")")
  )
}

# Run all outcomes
outcomes <- list(
  c("insured",       "insured"),
  c("worked12m",     "employed"),
  c("health_r",      "health_r"),
  c("flu_shot",      "flu_shot"),
  c("presc_med",     "presc_med"),
  c("wellvisit_1yr", "wellvisit_1yr"),
  c("cholcheck_1yr", "cholcheck_1yr"),
  c("bpcheck_1yr",   "bpcheck_1yr")
)

ols_results <- do.call(rbind, lapply(outcomes, function(o) run_ols(o[1], o[2])))
print(ols_results)

ols_results$Outcome <- c(
  "Insured",
  "Employed",
  "Self-Reported Health",
  "Flu Shot",
  "Prescription Meds",
  "Well Visit",
  "Cholesterol Check",
  "BP Check"
)

ols_results <- ols_results %>%
  mutate(
    est = Estimate,
    se = SE
  )

est_lookup <- setNames(ols_results$Estimate, ols_results$Outcome)
se_lookup  <- setNames(ols_results$SE, ols_results$Outcome)

order_vec <- c(
  "Insured",
  "Employed",
  "Self-Reported Health",
  "Flu Shot",
  "Prescription Meds",
  "Well Visit",
  "Cholesterol Check",
  "BP Check"
)

est_row <- paste(est_lookup[order_vec], collapse = " & ")
se_row  <- paste(se_lookup[order_vec], collapse = " & ")
col_names <- paste(order_vec, collapse = " & ")
col_nums  <- paste(paste0("(", seq_along(order_vec), ")"), collapse = " & ")

latex_rows <- ""
for (outcome in order_vec) {
  sub <- ols_results %>% filter(Outcome == outcome)
  
  latex_rows <- paste0(
    latex_rows,
    outcome, " & ", sub$est, " \\\\\n",
    " & ", sub$se, " \\\\\n\n"
  )
}

latex_table <- paste0(
  "\\begin{table}[htbp]\n",
  "\\centering\n",
  "\\caption{Manual Regression RD Estimates}\n",
  "\\label{tab:manual_rd}\n",
  "\\begin{threeparttable}\n",
  "\\small\n",
  "\\resizebox{\\textwidth}{!}{%\n",
  "\\begin{tabular}{l", paste(rep("c", length(order_vec)), collapse = ""), "}\n",
  "\\hline\\hline\n",
  " & ", col_names, " \\\\\n",
  " & ", col_nums, " \\\\\n",
  "\\hline\n",
  "RD Estimate & ", est_row, " \\\\\n",
  " & ", se_row, " \\\\\n",
  "\\hline\n",
  "\\end{tabular}%\n",
  "}\n",
  "\n",
  "\\begin{tablenotes}[flushleft]\n",
  "\\footnotesize\n",
  "\\item \\textit{Notes:} Each column reports the estimated discontinuity at age 65 from a manual OLS regression with the CCT-optimal bandwidths. Robust standard errors are reported in parentheses *** $p<0.01$, ** $p<0.05$, * $p<0.1$.\n",
  "\\end{tablenotes}\n",
  "\\end{threeparttable}\n",
  "\\end{table}\n"
)

cat(latex_table)