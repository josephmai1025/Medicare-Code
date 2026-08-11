source("Helpers.R") # import helper functions

ensure_package("haven")
ensure_package("dplyr")
ensure_package("rdrobust")
ensure_package("knitr")
ensure_package("kableExtra")

data_path <- find_data_file()
data <- read_dta(data_path)
data <- subset(data, astatflg == 1 & age >= 0 & age <= 100)

#### LOOKUP BANDWIDTHS ####
bw_lookup <- c(
  insured        = 5.47,
  health_r       = 7.20,
  flu_shot       = 5.68,
  presc_med      = 4.72,
  wellvisit_1yr  = 5.54,
  cholcheck_1yr  = 7.28,
  bpcheck_1yr    = 5.54
)
get_bw <- function(var_name) bw_lookup[[var_name]]

#### LABEL AND SAVE OUTCOMES ####
outcomes <- c("health_r", "flu_shot", "presc_med", "wellvisit_1yr", "cholcheck_1yr", "bpcheck_1yr")
outcome_labels <- c("Self-Reported Health", "Flu Shot", "Prescription Medication",
                    "Well Visit", "Cholesterol Check", "Blood Pressure Check")

#### SHARED SANDWICH-VARIANCE HELPER ####
sandwich_pieces <- function(fit_y, fit_d, temp) {
  X <- model.matrix(fit_y)
  w <- temp$sampweight
  uy <- resid(fit_y)
  ud <- resid(fit_d)
  
  XtWX_inv <- solve(t(X) %*% (w * X))
  
  meat_yy <- matrix(0, ncol(X), ncol(X))
  meat_dd <- matrix(0, ncol(X), ncol(X))
  meat_yd <- matrix(0, ncol(X), ncol(X))
  
  for (i in seq_len(nrow(X))) {
    xi <- matrix(X[i, ], ncol = 1)
    meat_yy <- meat_yy + (w[i]^2) * (uy[i]^2) * (xi %*% t(xi))
    meat_dd <- meat_dd + (w[i]^2) * (ud[i]^2) * (xi %*% t(xi))
    meat_yd <- meat_yd + (w[i]^2) * (uy[i] * ud[i]) * (xi %*% t(xi))
  }
  
  list(
    Vyy = XtWX_inv %*% meat_yy %*% XtWX_inv,
    Vdd = XtWX_inv %*% meat_dd %*% XtWX_inv,
    Vyd = XtWX_inv %*% meat_yd %*% XtWX_inv
  )
}

get_stars <- function(b, se) {
  if (is.na(b) || is.na(se) || se == 0 || is.infinite(se)) return("")
  p <- 2 * (1 - pnorm(abs(b / se)))
  if (is.na(p)) "" else if (p < 0.01) "***" else if (p < 0.05) "**" else if (p < 0.1) "*" else ""
}
fmt <- function(est, se) paste0(sprintf("%.3f", est), get_stars(est, se), " (", sprintf("%.3f", se), ")")

#### RACE ANALYSIS (White / Black / Other) ####

fit_rd_multi_edu <- function(df, outcome_var, enroll_var = "insured", h = NULL) {
  if (is.null(h)) h <- get_bw(outcome_var)
  temp <- df %>%
    filter(edu3 %in% c("Bachelors", "HSOrSomeCollege", "NoHS")) %>%
    mutate(
      g_bach = as.numeric(edu3 == "Bachelors"),
      g_nohs = as.numeric(edu3 == "NoHS")
    ) %>%
    filter(
      !is.na(.data[[outcome_var]]), !is.na(.data[[enroll_var]]),
      !is.na(age), !is.na(sampweight),
      !is.na(age65), !is.na(age_c), !is.na(age_x_treat)
    ) %>%
    filter(abs(age_c) <= h)
  
  form_y <- as.formula(paste0(outcome_var,
                              " ~ age65 + age_c + age_x_treat + g_bach + g_nohs + age65:g_bach + age65:g_nohs"))
  form_d <- as.formula(paste0(enroll_var,
                              " ~ age65 + age_c + age_x_treat + g_bach + g_nohs + age65:g_bach + age65:g_nohs"))
  
  fit_y <- lm(form_y, data = temp, weights = sampweight)
  fit_d <- lm(form_d, data = temp, weights = sampweight)
  
  pieces <- sandwich_pieces(fit_y, fit_d, temp)
  
  idx <- c("age65", "age65:g_bach", "age65:g_nohs")
  theta <- c(
    coef(fit_y)["age65"], coef(fit_y)["age65:g_bach"], coef(fit_y)["age65:g_nohs"],
    coef(fit_d)["age65"], coef(fit_d)["age65:g_bach"], coef(fit_d)["age65:g_nohs"]
  )
  Vtheta <- rbind(
    cbind(pieces$Vyy[idx, idx], pieces$Vyd[idx, idx]),
    cbind(t(pieces$Vyd[idx, idx]), pieces$Vdd[idx, idx])
  )
  list(theta = theta, Vtheta = Vtheta)
}

delta_wald_edu <- function(df, outcome_var, enroll_var = "insured", outcome_label = outcome_var) {
  fit_obj <- fit_rd_multi_edu(df, outcome_var, enroll_var = enroll_var)
  th <- fit_obj$theta; V <- fit_obj$Vtheta
  a0 <- th[1]; aBa <- th[2]; aN <- th[3]
  b0 <- th[4]; bBa <- th[5]; bN <- th[6]
  
  w_hscollege <- a0 / b0
  w_bach <- (a0 + aBa) / (b0 + bBa)
  w_nohs <- (a0 + aN) / (b0 + bN)
  
  make_grad <- function(num_idx, den_idx, extra_num = NULL, extra_den = NULL) {
    g <- rep(0, 6)
    num <- th[num_idx] + if (!is.null(extra_num)) th[extra_num] else 0
    den <- th[den_idx] + if (!is.null(extra_den)) th[extra_den] else 0
    g[num_idx] <- 1/den
    if (!is.null(extra_num)) g[extra_num] <- 1/den
    g[den_idx] <- -num/(den^2)
    if (!is.null(extra_den)) g[extra_den] <- -num/(den^2)
    g
  }
  g_hscollege <- make_grad(1, 4)
  g_bach <- make_grad(1, 4, 2, 5)
  g_nohs <- make_grad(1, 4, 3, 6)
  
  se <- function(g) sqrt(as.numeric(t(g) %*% V %*% g))
  se_hscollege <- se(g_hscollege); se_bach <- se(g_bach); se_nohs <- se(g_nohs)
  
  diff_bach <- w_bach - w_hscollege; se_diff_bach <- se(g_bach - g_hscollege)
  diff_nohs <- w_nohs - w_hscollege; se_diff_nohs <- se(g_nohs - g_hscollege)
  
  data.frame(
    Outcome = outcome_label,
    HSOrSomeCollege = fmt(w_hscollege, se_hscollege),
    Bachelors = fmt(w_bach, se_bach),
    NoHS = fmt(w_nohs, se_nohs),
    Diff_Bachelors = fmt(diff_bach, se_diff_bach),
    Diff_NoHS = fmt(diff_nohs, se_diff_nohs),
    stringsAsFactors = FALSE
  )
}

#### FIRST STAGE TABLE ### 

first_stage_row_edu <- function(df, enroll_var = "insured") {
  fit_obj <- fit_rd_multi_edu(df, enroll_var, enroll_var = enroll_var)
  th <- fit_obj$theta; V <- fit_obj$Vtheta
  b0 <- th[4]; bBa <- th[5]; bN <- th[6]
  
  fs_hscollege <- b0
  fs_bach <- b0 + bBa
  fs_nohs <- b0 + bN
  
  g_hscollege <- c(0, 0, 0, 1, 0, 0)
  g_bach <- c(0, 0, 0, 1, 1, 0)
  g_nohs <- c(0, 0, 0, 1, 0, 1)
  
  se <- function(g) sqrt(as.numeric(t(g) %*% V %*% g))
  se_hscollege <- se(g_hscollege); se_bach <- se(g_bach); se_nohs <- se(g_nohs)
  
  diff_bach <- fs_bach - fs_hscollege; se_diff_bach <- se(g_bach - g_hscollege)
  diff_nohs <- fs_nohs - fs_hscollege; se_diff_nohs <- se(g_nohs - g_hscollege)
  
  data.frame(
    Outcome = "Insured (First Stage)",
    HSOrSomeCollege = fmt(fs_hscollege, se_hscollege),
    Bachelors = fmt(fs_bach, se_bach),
    NoHS = fmt(fs_nohs, se_nohs),
    Diff_Bachelors = fmt(diff_bach, se_diff_bach),
    Diff_NoHS = fmt(diff_nohs, se_diff_nohs),
    stringsAsFactors = FALSE
  )
}

make_first_stage_table_edu <- function(df, caption_text) {
  fs <- first_stage_row_edu(df)
  fs_row <- paste0(fs$Outcome, " & ", fs$HSOrSomeCollege, " & ", fs$Bachelors, " & ", fs$NoHS,
                   " & ", fs$Diff_Bachelors, " & ", fs$Diff_NoHS, " \\\\\n")
  
  latex_table <- paste0(
    "\\begin{table}[htbp]\n",
    "\\centering\n",
    "\\caption{", caption_text, "}\n",
    "\\label{tab:first_stage_edu}\n",
    "\\small\n",
    "\\begin{tabular}{lccccc}\n",
    "\\hline\\hline\n",
    " & HS/Some College & Bachelors+ & No HS & Diff (Bachelors+) & Diff (No HS) \\\\\n",
    "\\hline\n",
    fs_row,
    "\\hline\n",
    "\\end{tabular}\n",
    "\\begin{tablenotes}\n",
    "\\footnotesize\n",
    "\\item \\textit{Notes:} First-stage discontinuity in insurance coverage at age 65, by ",
    "educational attainment. Standard errors computed via the delta method in parentheses. ",
    "*** $p<0.01$, ** $p<0.05$, * $p<0.1$.\n",
    "\\end{tablenotes}\n",
    "\\end{table}\n"
  )
  cat(latex_table)
}

make_first_stage_table_edu(data, "First-Stage Discontinuity in Insurance Coverage by Educational Attainment")

#### RACE LATE ESTIMATES ####

make_edu_table_combined <- function(df, caption_text) {
  rows <- lapply(seq_along(outcomes), function(i) {
    delta_wald_edu(df, outcomes[i], outcome_label = outcome_labels[i])
  })
  table_out <- do.call(rbind, rows)
  
  kable(table_out, format = "html", caption = caption_text,
        align = c("l", rep("c", 5)), booktabs = TRUE, escape = FALSE,
        col.names = c("Outcome", "Bachelors+", "No HS", "HS/Some College",
                      "Diff (No HS)", "Diff (HS/Some College)")) %>%
    kable_styling(full_width = FALSE, position = "center", font_size = 13,
                  bootstrap_options = c("striped", "hover")) %>%
    row_spec(0, bold = TRUE)
}

make_edu_table_latex <- function(df, caption_text) {
  rows <- lapply(seq_along(outcomes), function(i) {
    delta_wald_edu(df, outcomes[i], outcome_label = outcome_labels[i])
  })
  table_df <- do.call(rbind, rows)
  
  body <- paste(sapply(1:nrow(table_df), function(i) {
    row <- table_df[i, ]
    est_row <- paste(row$HSOrSomeCollege, row$Bachelors, row$NoHS,
                     row$Diff_Bachelors, row$Diff_NoHS, sep = " & ")
    paste0(row$Outcome, " & ", est_row, " \\\\\n")
  }), collapse = "")
  
  latex_table <- paste0(
    "\\begin{table}[htbp]\n",
    "\\centering\n",
    "\\caption{", caption_text, "}\n",
    "\\label{tab:edu_wald}\n",
    "\\begin{threeparttable}\n",
    "\\small\n",
    "\\resizebox{\\textwidth}{!}{%\n",
    "\\begin{tabular}{lccccc}\n",
    "\\hline\\hline\n",
    "Outcome & HS/Some College & Bachelors+ & No HS & Diff (Bachelors+) & Diff (No HS) \\\\\n",
    "\\hline\n",
    body,
    "\\hline\n",
    "\\end{tabular}%\n",
    "}\n",
    "\\begin{tablenotes}\n",
    "\\footnotesize\n",
    "\\item \\textit{Notes:} LATE estimates by educational attainment at age 65, with ",
    "differences from HS/Some College reported using the delta method. ",
    "Robust standard errors in parentheses. *** $p<0.01$, ** $p<0.05$, * $p<0.1$.\n",
    "\\end{tablenotes}\n",
    "\\end{threeparttable}\n",
    "\\end{table}\n"
  )
  cat(latex_table)
}
make_edu_table_latex(data, "Wald Estimates at Age 65 by Educational Attainment")