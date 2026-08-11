source("Helpers.R") # import helper functions

ensure_package("haven")
ensure_package("rdrobust")
ensure_package("dplyr")
ensure_package("tidyr")

data_path <- find_data_file()
data <- read_dta(data_path)
data <- subset(data, astatflg == 1 & age >= 0 & age <= 100)

#### Helper RD Function ####

run_rd <- function(y_var, var_name, fuzzy_var) {
  rd <- rdrobust(
    y_var,
    data$age,
    c = 65,
    p = 1,
    bwselect = "mserd",
    masspoints = "adjust",
    weights = data$sampweight
  )
  
  est <- rd$coef[2]
  se  <- rd$se[3]
  z   <- est / se
  p   <- 2 * (1 - pnorm(abs(z)))
  
  stars <- ifelse(p < 0.01, "***",
                  ifelse(p < 0.05, "**",
                         ifelse(p < 0.1, "*", "")))
  
  return(data.frame(
    Outcome = var_name,
    Bandwidth = round(rd$bws[1, 1], 2),
    Estimate = round(est, 3),
    SE = round(se, 3),
    Stars = stars
  ))
}

#### Run First Stage + ITTs ####
outcomes <- list(
  insured = data$insured,
  employed = data$employed,
  health_r = data$health_r,
  flu_shot = data$flu_shot,
  presc_med = data$presc_med,
  wellvisit_1yr = data$wellvisit_1yr,
  cholcheck_1yr = data$cholcheck_1yr,
  bpcheck_1yr = data$bpcheck_1yr
)

results_list <- list()

for (name in names(outcomes)) {
  res <- run_rd(outcomes[[name]], name)
  results_list[[length(results_list) + 1]] <- res
}

results_df <- do.call(rbind, results_list)

#### Create table for ITTs excluding employment ####

# Rename variables
itt_labels <- c(
  insured = "Insured",
  health_r = "Self-Reported Health",
  flu_shot = "Flu Shot",
  presc_med = "Prescription Medication",
  wellvisit_1yr = "Well Visit",
  cholcheck_1yr = "Cholesterol Check",
  bpcheck_1yr = "Blood Pressure Check"
)

results_df_labeled <- results_df
results_df_labeled$Outcome <- ifelse(
  results_df_labeled$Outcome %in% names(itt_labels),
  itt_labels[results_df_labeled$Outcome],
  results_df_labeled$Outcome
)

itt_7 <- results_df_labeled %>%
  filter(Outcome != "employed") %>%
  mutate(
    est = paste0(Estimate, Stars),
    se  = paste0("(", SE, ")"),
    bw  = as.character(Bandwidth)
  )

order_vec <- c(
  "Insured",
  "Self-Reported Health",
  "Flu Shot",
  "Prescription Medication",
  "Well Visit",
  "Cholesterol Check",
  "Blood Pressure Check"
)

est_lookup <- setNames(itt_7$est, itt_7$Outcome)
se_lookup  <- setNames(itt_7$se, itt_7$Outcome)
bw_lookup  <- setNames(itt_7$bw, itt_7$Outcome)

est_row <- paste(est_lookup[order_vec], collapse = " & ")
se_row  <- paste(se_lookup[order_vec], collapse = " & ")
bw_row  <- paste(bw_lookup[order_vec], collapse = " & ")
col_names <- paste(order_vec, collapse = " & ")
col_nums  <- paste(paste0("(", seq_along(order_vec), ")"), collapse = " & ")

latex_table <- paste0(
  "\\begin{table}[htbp]\n",
  "\\centering\n",
  "\\caption{Intent-to-Treat Effects at Age 65 (CCT-Optimal Bandwidth)}\n",
  "\\resizebox{\\textwidth}{!}{%\n",
  "\\begin{tabular}{l", paste(rep("c", length(order_vec)), collapse = ""), "}\n",
  "\\hline\\hline\n",
  " & ", col_names, " \\\\\n",
  " & ", col_nums, " \\\\\n",
  "\\hline\n",
  "Medicare Eligibility & ", est_row, " \\\\\n",
  " & ", se_row, " \\\\\n",
  "Bandwidth (h) & ", bw_row, " \\\\\n",
  "\\hline\n",
  "\\end{tabular}%\n",
  "}\n",
  "\\end{table}\n"
)
cat(latex_table)

#### Create RD plots ####
# install.packages("showtext")

get_bw <- function(var_name) {
  results_df$Bandwidth[results_df$Outcome == var_name][1]
}

pdf("rd_plot_for_insurance.pdf", width = 10, height = 5)

rdplot(
  data$insured,
  data$age,
  c = 65,
  h = get_bw("insured"),
  p = 1,
  x.lim = c(45, 85),
  x.label = "Age",
  y.label = "Insurance Coverage",
  title = "Effect of Medicare Eligibility on Insurance at Age 65"
)

dev.off()

pdf("rd_plot_for_employment.pdf", width = 10, height = 5)
rdplot(
  data$employed,
  data$age,
  c = 65,
  h = get_bw("employed"),
  p = 1,
  x.lim = c(45, 85),
  x.label = "Age",
  y.label = "Employed",
  title = "Employment Discontinuity Plot"
)
dev.off()

pdf("rd_plot_for_health_status.pdf", width = 10, height = 5)

rdplot(
  data$health_r,
  data$age,
  c = 65,
  h = get_bw("health_r"),
  x.lim = c(45, 85),
  p = 1,
  x.label = "Age",
  y.label = "Self-Reported Health",
  title = "Effect of Medicare Eligibility on Self-Reported Health at Age 65"
)

dev.off()

pdf("rd_plot_for_flu_shot.pdf", width = 10, height = 5)

rdplot(
  data$flu_shot,
  data$age,
  c = 65,
  h = get_bw("flu_shot"),
  x.lim = c(45, 85),
  p = 1,
  x.label = "Age",
  y.label = "Flu Shot",
  title = "Effect of Medicare Eligibility on Flu Shot at Age 65"
)

dev.off()

pdf("rd_plot_for_meds.pdf", width = 10, height = 5)

rdplot(
  data$presc_med,
  data$age,
  c = 65,
  h = get_bw("presc_med"),
  x.lim = c(45, 85),
  p = 1,
  x.label = "Age",
  y.label = "Prescription Medication",
  title = "Effect of Medicare Eligibility on Prescription Medication at Age 65"
)

dev.off()

pdf("rd_plot_for_well_visit.pdf", width = 10, height = 5)

rdplot(
  data$wellvisit_1yr,
  data$age,
  c = 65,
  h = get_bw("wellvisit_1yr"),
  x.lim = c(45, 85),
  p = 1,
  x.label = "Age",
  y.label = "Wellness Visit",
  title = "Effect of Medicare Eligibility on Wellness Visit at Age 65"
)

dev.off()

pdf("rd_plot_for_chol_check.pdf", width = 10, height = 5)

rdplot(
  data$cholcheck_1yr,
  data$age,
  c = 65,
  h = get_bw("cholcheck_1yr"),
  x.lim = c(45, 85),
  p = 1,
  x.label = "Age",
  y.label = "Cholesterol Check",
  title = "Effect of Medicare Eligibility on Cholesterol Check at Age 65"
)

dev.off()

pdf("rd_plot_for_bp_check.pdf", width = 10, height = 5)

rdplot(
  data$bpcheck_1yr,
  data$age,
  c = 65,
  h = get_bw("bpcheck_1yr"),
  x.lim = c(45, 85),
  p = 1,
  x.label = "Age",
  y.label = "Blood Pressure Check",
  title = "Effect of Medicare Eligibility on Blood Pressure Check at Age 65"
)

dev.off()

#### Create table for employment (CCT-optimal bandwidth) ####
employment_df <- subset(results_df, Outcome == "employed")
print(employment_df)

est_row <- paste0(employment_df$Estimate, employment_df$Stars)
se_row  <- paste0("(", employment_df$SE, ")")
bw_row  <- as.character(employment_df$Bandwidth)

latex_table <- paste0(
  "\\begin{table}[ht]\n",
  "\\centering\n",
  "\\caption{Employment Effects at Age 65 (CCT-Optimal Bandwidth)}\n",
  "\\small\n",
  "\\begin{tabular}{lc}\n",
  "\\hline\\hline\n",
  "Outcome & Employed \\\\\n",
  "\\hline\n",
  "Medicare Eligibility & ", est_row, " \\\\\n",
  " & ", se_row, " \\\\\n",
  "Bandwidth (h) & ", bw_row, " \\\\\n",
  "\\hline\n",
  "\\end{tabular}\n",
  "\\begin{flushleft}\n",
  "\\footnotesize\n",
  "\\textit{Notes:} Reports the estimated regression discontinuity in employment at age 65 using the CCT MSE-optimal bandwidth. Standard errors are reported in parentheses. ***
$p<0.01$, ** $p<0.05$, * $p<0.1$.\n",
  "\\end{flushleft}\n",
  "\\end{table}"
)
cat(latex_table)

#### LATE ####

late_outcomes <- list(
  health_r      = data$health_r,
  flu_shot      = data$flu_shot,
  presc_med     = data$presc_med,
  wellvisit_1yr = data$wellvisit_1yr,
  cholcheck_1yr = data$cholcheck_1yr,
  bpcheck_1yr   = data$bpcheck_1yr
)

# Run LATE
late_results <- list()
for (name in names(late_outcomes)) {
  res <- run_fuzzy_rd(late_outcomes[[name]], name, data$insured)
  late_results[[length(late_results) + 1]] <- res
}

late_df <- do.call(rbind, late_results)

late_df$Estimate_SE <- paste0(
  late_df$Estimate,
  late_df$Stars,
  " (", late_df$SE, ")"
)

print(late_df)

#### LATE Table ####

late_labels <- c(
  health_r = "Self-Reported Health",
  flu_shot = "Flu Shot",
  presc_med = "Prescription Medication",
  wellvisit_1yr = "Well Visit",
  cholcheck_1yr = "Cholesterol Check",
  bpcheck_1yr = "Blood Pressure Check"
)

late_df$Outcome <- ifelse(
  late_df$Outcome %in% names(late_labels),
  late_labels[late_df$Outcome],
  late_df$Outcome
)

order_vec <- c(
  "Self-Reported Health",
  "Flu Shot",
  "Prescription Medication",
  "Well Visit",
  "Cholesterol Check",
  "Blood Pressure Check"
)

late_lookup <- setNames(late_df$Estimate_SE, late_df$Outcome)
late_row <- paste(late_lookup[order_vec], collapse = " & ")
col_names <- paste(order_vec, collapse = " & ")
col_nums <- paste(paste0("(", seq_along(order_vec), ")"), collapse = " & ")
bw_lookup_late <- setNames(as.character(late_df$Bandwidth), late_df$Outcome)
bw_row_late <- paste(bw_lookup_late[order_vec], collapse = " & ")

latex_table <- paste0(
  "\\begin{table}[htbp]\n",
  "\\centering\n",
  "\\caption{Local Average Treatment Effects at Age 65 (CCT-Optimal Bandwidth)}\n",
  "\\resizebox{\\textwidth}{!}{%\n",
  "\\begin{tabular}{l", paste(rep("c", length(order_vec)), collapse = ""), "}\n",
  "\\hline\\hline\n",
  " & ", col_names, " \\\\\n",
  " & ", col_nums, " \\\\\n",
  "\\hline\n",
  "Medicare Coverage & ", late_row, " \\\\\n",
  "Bandwidth (h) & ", bw_row_late, " \\\\\n",
  "\\hline\n",
  "\\end{tabular}%\n",
  "}\n",
  "\\begin{flushleft}\n",
  "\\footnotesize\n",
  "\\textit{Notes:} Robust standard errors are reported in parentheses. Bandwidths selected via the CCT MSE-optimal procedure. *** $p<0.01$, ** $p<0.05$, * $p<0.1$.\n",
  "\\end{flushleft}\n",
  "\\end{table}\n"
)
cat(latex_table)

late_table <- pivot_wider(
  late_df,
  id_cols = Outcome,
  names_from = Bandwidth,
  values_from = Estimate_SE
)

print(late_table)

#### SUMMARY STATISTICS ####

sumstat_vars <- c("insured", "health_r", "flu_shot", "presc_med",
                  "wellvisit_1yr", "cholcheck_1yr", "bpcheck_1yr")
sumstat_labels <- c("Insured", "Self-Reported Health", "Flu Shot",
                    "Prescription Medication", "Well Visit", "Cholesterol Check",
                    "Blood Pressure Check")

data_bw <- data %>% filter(abs(age_c) <= 7.28)

group_vals <- c("HSOrSomeCollege", "Bachelors", "NoHS")
group_labels <- c(HSOrSomeCollege = "HS/Some College", Bachelors = "Bachelors+", NoHS = "No HS")

#### COMPUTE MEANS BY GROUP x PERIOD ####
compute_stats <- function(df, period_label) {
  df %>%
    filter(edu3 %in% group_vals) %>%
    group_by(edu3) %>%
    summarise(across(all_of(sumstat_vars),
                     list(mean = ~round(mean(.x, na.rm = TRUE), 3),
                          sd   = ~round(sd(.x, na.rm = TRUE), 3)),
                     .names = "{.col}_{.fn}"),
              N = n(), .groups = "drop") %>%
    mutate(period = period_label)
}

stats_before  <- compute_stats(data_bw %>% filter(age65 == 0), "Before")
stats_after   <- compute_stats(data_bw %>% filter(age65 == 1), "After")
stats_overall <- compute_stats(data_bw, "Overall")

sumstat_wide <- bind_rows(stats_before, stats_after, stats_overall)

print(sumstat_wide)

#### BUILD VERTICAL LATEX TABLE ####
build_sumstat_latex_vertical_edu <- function(df, group_vals, period_order, caption_text, label_text) {
  # header: group names spanning 3 columns each
  group_header <- paste(sapply(group_vals, function(g) {
    n_val <- df %>% filter(edu3 == g, period == "Overall") %>% pull(N)
    paste0("\\multicolumn{3}{c}{", group_labels[[g]], " (N=", format(n_val, big.mark=","), ")}")
  }), collapse = " & ")
  
  period_header <- paste(rep(period_order, length(group_vals)), collapse = " & ")
  
  body <- ""
  for (i in seq_along(sumstat_vars)) {
    v <- sumstat_vars[i]
    mean_cells <- c()
    sd_cells <- c()
    for (g in group_vals) {
      for (p in period_order) {
        row <- df %>% filter(edu3 == g, period == p)
        mean_cells <- c(mean_cells, sprintf("%.3f", row[[paste0(v, "_mean")]]))
        sd_cells   <- c(sd_cells, sprintf("(%.3f)", row[[paste0(v, "_sd")]]))
      }
    }
    body <- paste0(body, sumstat_labels[i], " & ", paste(mean_cells, collapse = " & "), " \\\\\n")
    body <- paste0(body, " & ", paste(sd_cells, collapse = " & "), " \\\\\n")
  }
  
  n_cols <- length(group_vals) * length(period_order)
  cmidrules <- paste(sapply(seq_along(group_vals), function(i) {
    start <- (i - 1) * length(period_order) + 2
    end <- start + length(period_order) - 1
    paste0("\\cmidrule(lr){", start, "-", end, "}")
  }), collapse = " ")
  
  n_cells <- c()
  for (g in group_vals) {
    for (p in period_order) {
      row <- df %>% filter(edu3 == g, period == p)
      n_cells <- c(n_cells, format(row$N, big.mark = ","))
    }
  }
  body <- paste0(body, "N & ", paste(n_cells, collapse = " & "), " \\\\\n")
  
  latex_table <- paste0(
    "\\begin{table}[htbp]\n",
    "\\centering\n",
    "\\caption{", caption_text, "}\n",
    "\\label{", label_text, "}\n",
    "\\resizebox{\\textwidth}{!}{%\n",
    "\\small\n",
    "\\begin{tabular}{l", paste(rep("c", n_cols), collapse = ""), "}\n",
    "\\toprule\n",
    " & ", group_header, " \\\\\n",
    cmidrules, "\n",
    " & ", period_header, " \\\\\n",
    "\\midrule\n",
    body,
    "\\bottomrule\n",
    "\\end{tabular}%\n",
    "}\n",
    "\\begin{tablenotes}\n",
    "\\footnotesize\n",
    "\\item \\textit{Notes:} Sample restricted to observations within the widest CCT-optimal bandwidth (7.28 years) around age 65. Before/After columns split the sample at the Medicare eligibility threshold; Overall pools both periods.\n",
    "\\end{tablenotes}\n",
    "\\end{table}\n"
  )
  cat(latex_table)
}

build_sumstat_latex_vertical_edu(sumstat_wide, group_vals, c("Before", "After", "Overall"),
                                 "Summary Statistics by Educational Attainment", "tab:sumstat_edu")