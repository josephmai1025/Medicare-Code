source("Helpers.R")

ensure_package("dplyr")
ensure_package("knitr")
ensure_package("kableExtra")
ensure_package("rmarkdown")
ensure_package("tinytex")

# -----------------------------------------------------------------------------
# Load and prepare the analysis sample
# -----------------------------------------------------------------------------
data_path <- find_data_file()
data <- load_data_file(data_path)
data <- subset(data, astatflg == 1 & age >= 0 & age <= 100)
data$above_65 <- as.numeric(data$age >= 65)
data$age_c <- data$age - 65

# -----------------------------------------------------------------------------
# Bandwidths used in the project
# -----------------------------------------------------------------------------
bw_lookup <- c(
  insured        = 5.47,
  employed       = 5.47,
  health_r       = 7.20,
  flu_shot       = 5.68,
  presc_med      = 4.72,
  wellvisit_1yr  = 5.54,
  cholcheck_1yr  = 7.28,
  bpcheck_1yr    = 5.54,
  worked12m      = 5.47,
  worked12m_1yr  = 5.47
)
get_bw <- function(var_name) {
  if (is.null(bw_lookup[[var_name]])) {
    stop(paste("No bandwidth defined for variable:", var_name))
  }
  bw_lookup[[var_name]]
}

# -----------------------------------------------------------------------------
# Standardize education groups across project versions
# -----------------------------------------------------------------------------
standardize_edu3 <- function(df) {
  if ("edu3" %in% names(df)) {
    x <- df$edu3
    if (all(x %in% c("Bachelors", "HSOrSomeCollege", "NoHS"))) {
      return(df)
    }
    if (all(x %in% c("Bachelors+", "HS/Some College", "No HS"))) {
      df$edu3 <- dplyr::recode(
        df$edu3,
        "Bachelors+" = "Bachelors",
        "HS/Some College" = "HSOrSomeCollege",
        "No HS" = "NoHS"
      )
      return(df)
    }
    if (all(x %in% c("Bachelors+", "HS/Some College", "No HS", NA))) {
      df$edu3 <- dplyr::recode(
        df$edu3,
        "Bachelors+" = "Bachelors",
        "HS/Some College" = "HSOrSomeCollege",
        "No HS" = "NoHS"
      )
      return(df)
    }
  }

  if ("educ" %in% names(df)) {
    df <- df %>%
      mutate(
        edu3 = case_when(
          educ >= 100 & educ < 200 ~ "NoHS",
          educ >= 201 & educ < 400 ~ "HSOrSomeCollege",
          educ >= 400 & educ < 996 ~ "Bachelors",
          TRUE ~ NA_character_
        )
      )
    return(df)
  }

  df
}

data <- standardize_edu3(data)

# -----------------------------------------------------------------------------
# Helper for retrieving a full OLS coefficient table and fit metrics
# -----------------------------------------------------------------------------
get_ols_table <- function(fit, model_name, outcome_var, subgroup = NA_character_) {
  s <- summary(fit)
  coef_df <- as.data.frame(s$coefficients)
  coef_df$term <- rownames(coef_df)

  ci_df <- as.data.frame(confint(fit))
  ci_df$term <- rownames(ci_df)

  coef_df <- merge(coef_df, ci_df, by = "term", all.x = TRUE)
  coef_df$model_name <- model_name
  coef_df$outcome_var <- outcome_var
  coef_df$subgroup <- subgroup

  names(coef_df)[names(coef_df) == "Estimate"] <- "estimate"
  names(coef_df)[names(coef_df) == "Std. Error"] <- "std.error"
  names(coef_df)[names(coef_df) == "t value"] <- "t.value"
  names(coef_df)[names(coef_df) == "Pr(>|t|)"] <- "p.value"
  names(coef_df)[names(coef_df) == "2.5 %"] <- "conf.low"
  names(coef_df)[names(coef_df) == "97.5 %"] <- "conf.high"

  coef_df <- coef_df[, c(
    "model_name", "subgroup", "outcome_var", "term",
    "estimate", "std.error", "t.value", "p.value", "conf.low", "conf.high"
  )]

  y <- model.response(model.frame(fit))
  residuals_vec <- residuals(fit)
  n <- nobs(fit)
  p <- length(coef(fit))
  ss_res <- sum(residuals_vec^2, na.rm = TRUE)
  ss_tot <- sum((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE)

  r_squared <- if (!is.null(s$r.squared) && !is.na(s$r.squared)) {
    s$r.squared
  } else if (ss_tot > 0) {
    1 - ss_res / ss_tot
  } else {
    NA_real_
  }

  df_res <- if (!is.null(fit$df.residual)) fit$df.residual else max(n - p, 1)
  adj_r_squared <- if (!is.null(s$adj.r.squared) && !is.na(s$adj.r.squared)) {
    s$adj.r.squared
  } else if (!is.na(r_squared) && n > p + 1) {
    1 - (1 - r_squared) * (n - 1) / (n - p - 1)
  } else {
    NA_real_
  }

  sigma <- if (!is.null(s$sigma) && !is.na(s$sigma)) {
    s$sigma
  } else if (df_res > 0) {
    sqrt(ss_res / df_res)
  } else {
    NA_real_
  }

  f_stat <- if (!is.null(s$fstatistic) && length(s$fstatistic) >= 3) {
    as.numeric(s$fstatistic[1])
  } else {
    NA_real_
  }

  f_p <- if (!is.null(s$fstatistic) && length(s$fstatistic) >= 3) {
    pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3], lower.tail = FALSE)
  } else {
    NA_real_
  }

  fit_stats <- data.frame(
    model_name = model_name,
    subgroup = subgroup,
    outcome_var = outcome_var,
    nobs = n,
    df_model = if (!is.null(s$df[1])) s$df[1] else p - 1,
    df_residual = df_res,
    r_squared = r_squared,
    adj_r_squared = adj_r_squared,
    sigma = sigma,
    f_statistic = f_stat,
    f_p_value = f_p,
    aic = AIC(fit),
    bic = BIC(fit),
    stringsAsFactors = FALSE
  )

  list(coefficients = coef_df, fit_stats = fit_stats)
}

# -----------------------------------------------------------------------------
# Education regressions from Education.R
# -----------------------------------------------------------------------------
education_outcomes <- c("health_r", "flu_shot", "presc_med", "wellvisit_1yr", "cholcheck_1yr", "bpcheck_1yr")

run_education_ols <- function(df, outcome_var, enroll_var = "insured") {
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
    filter(abs(age_c) <= get_bw(outcome_var))

  if (nrow(temp) == 0L) {
    stop(paste("No observations left for education model:", outcome_var))
  }

  form_y <- as.formula(paste0(
    outcome_var,
    " ~ age65 + age_c + age_x_treat + g_bach + g_nohs + age65:g_bach + age65:g_nohs"
  ))
  form_d <- as.formula(paste0(
    enroll_var,
    " ~ age65 + age_c + age_x_treat + g_bach + g_nohs + age65:g_bach + age65:g_nohs"
  ))

  fit_y <- lm(form_y, data = temp, weights = sampweight)
  fit_d <- lm(form_d, data = temp, weights = sampweight)

  list(
    outcome = get_ols_table(fit_y, model_name = "education_outcome", outcome_var = outcome_var, subgroup = "all"),
    first_stage = get_ols_table(fit_d, model_name = "education_first_stage", outcome_var = enroll_var, subgroup = "all")
  )
}

education_tables <- lapply(education_outcomes, function(v) run_education_ols(data, v, enroll_var = "insured"))
education_coeffs <- do.call(rbind, lapply(education_tables, function(x) x$outcome$coefficients))
education_fit <- do.call(rbind, lapply(education_tables, function(x) x$outcome$fit_stats))
education_fs_coeffs <- do.call(rbind, lapply(education_tables, function(x) x$first_stage$coefficients))
education_fs_fit <- do.call(rbind, lapply(education_tables, function(x) x$first_stage$fit_stats))

# -----------------------------------------------------------------------------
# Manual RD OLS regressions from Robustness and Diagnostics.R
# -----------------------------------------------------------------------------
manual_outcomes <- list(
  c("insured", "insured"),
  c("employed", "employed"),
  c("health_r", "health_r"),
  c("flu_shot", "flu_shot"),
  c("presc_med", "presc_med"),
  c("wellvisit_1yr", "wellvisit_1yr"),
  c("cholcheck_1yr", "cholcheck_1yr"),
  c("bpcheck_1yr", "bpcheck_1yr")
)

run_manual_rd_ols <- function(y_var, var_name) {
  h <- get_bw(var_name)
  data_bw <- data[abs(data$age_c) <= h, ]
  if (nrow(data_bw) == 0L) {
    stop(paste("No observations for manual RD model:", var_name))
  }
  fml <- as.formula(paste(y_var, "~ above_65 + age_c + age_c:above_65"))
  fit <- lm(fml, data = data_bw, weights = sampweight)
  get_ols_table(fit, model_name = "manual_rd_ols", outcome_var = y_var, subgroup = var_name)
}

manual_rd_tables <- lapply(manual_outcomes, function(o) run_manual_rd_ols(o[1], o[2]))
manual_rd_coeffs <- do.call(rbind, lapply(manual_rd_tables, function(x) x$coefficients))
manual_rd_fit <- do.call(rbind, lapply(manual_rd_tables, function(x) x$fit_stats))

# -----------------------------------------------------------------------------
# Put the model results in a consistent order and print tables
# -----------------------------------------------------------------------------
pretty_ols_value <- function(x) {
  x <- as.character(x)
  replacements <- c(
    "education_outcome" = "Education outcome",
    "education_first_stage" = "Education first stage",
    "manual_rd_ols" = "Manual RD OLS",
    "Bachelors+" = "Bachelor's+",
    "Bachelors" = "Bachelor's",
    "HSOrSomeCollege" = "HS / Some College",
    "HS/Some College" = "HS / Some College",
    "NoHS" = "No HS",
    "No HS" = "No HS",
    "insured" = "Insured",
    "employed" = "Employed",
    "health_r" = "Self-reported health",
    "flu_shot" = "Influenza vaccination",
    "presc_med" = "Prescription medication use",
    "wellvisit_1yr" = "Annual wellness visit",
    "cholcheck_1yr" = "Cholesterol screening",
    "bpcheck_1yr" = "Blood pressure screening",
    "age65:g_bach" = "Age 65+ × Bachelor's",
    "age65:g_nohs" = "Age 65+ × No HS",
    "age65:g_mohs" = "Age 65+ × No HS",
    "g_bach" = "Bachelor's group",
    "g_nohs" = "No HS group",
    "g_mohs" = "No HS group",
    "age_x_treat" = "Age centered × 65+",
    "age_c:above_65" = "Age centered × 65+",
    "above_65:age_c" = "65+ × Age centered",
    "above_65" = "65+ indicator",
    "age65" = "Age 65+",
    "age_c" = "Age centered",
    "Intercept" = "Intercept"
  )

  key_order <- names(replacements)[order(nchar(names(replacements)), decreasing = TRUE)]
  for (nm in key_order) {
    x <- gsub(nm, replacements[[nm]], x, fixed = TRUE)
  }

  x <- gsub(":", " × ", x, fixed = TRUE)
  x <- gsub("_", " ", x)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

pretty_ols_table <- function(df) {
  out <- df
  col_map <- c(
    "model_name" = "Model Name",
    "subgroup" = "Subgroup",
    "outcome_var" = "Outcome Variable",
    "term" = "Term",
    "estimate" = "Estimate",
    "std.error" = "Std. Error",
    "t.value" = "t Statistic",
    "p.value" = "p-value",
    "conf.low" = "CI Lower",
    "conf.high" = "CI Upper",
    "nobs" = "N",
    "df_model" = "Model df",
    "df_residual" = "Residual df",
    "r_squared" = "R$^{2}$",
    "adj_r_squared" = "Adjusted $R^{2}$",
    "sigma" = "Residual SD",
    "f_statistic" = "F Statistic",
    "f_p_value" = "F $p$-value",
    "aic" = "AIC",
    "bic" = "BIC"
  )

  names(out) <- ifelse(names(out) %in% names(col_map), col_map[names(out)], names(out))

  for (nm in names(out)) {
    if (is.character(out[[nm]])) {
      out[[nm]] <- pretty_ols_value(out[[nm]])
    }
  }

  rownames(out) <- NULL
  out
}

print_ols_table <- function(df, digits = 4) {
  print(knitr::kable(pretty_ols_table(df), digits = digits, format = "markdown", row.names = FALSE, escape = FALSE))
}

cat("\n=== Education OLS coefficient table ===\n")
print_ols_table(education_coeffs)

cat("\n=== Education OLS fit statistics ===\n")
print_ols_table(education_fit)

cat("\n=== Education first-stage coefficient table ===\n")
print_ols_table(education_fs_coeffs)

cat("\n=== Education first-stage fit statistics ===\n")
print_ols_table(education_fs_fit)

cat("\n=== Manual RD OLS coefficient table ===\n")
print_ols_table(manual_rd_coeffs)

cat("\n=== Manual RD OLS fit statistics ===\n")
print_ols_table(manual_rd_fit)

# -----------------------------------------------------------------------------
# Publication-ready export: build a single PDF document with all tables
# -----------------------------------------------------------------------------
output_dir <- "output"
dir.create(output_dir, showWarnings = FALSE)

format_table_for_pdf <- function(df, caption_text, digits = 4) {
  df_out <- pretty_ols_table(df)
  numeric_cols <- intersect(names(df_out), c("Estimate", "Std. Error", "t Statistic", "p-value", "CI Lower", "CI Upper", "R$^{2}$", "Adjusted $R^{2}$", "Residual SD", "F Statistic", "F $p$-value", "AIC", "BIC", "N", "Model df", "Residual df"))
  if (length(numeric_cols) > 0) {
    df_out[numeric_cols] <- lapply(df_out[numeric_cols], function(x) if (is.numeric(x)) round(x, digits) else x)
  }

  knitr::kable(
    df_out,
    format = "latex",
    digits = digits,
    caption = caption_text,
    booktabs = TRUE,
    escape = FALSE,
    row.names = FALSE,
    longtable = FALSE
  ) %>%
    kableExtra::kable_styling(
      latex_options = c("hold_position"),
      font_size = 9
    )
}

report_path <- file.path(output_dir, "ols_regression_tables_report.Rmd")

report_lines <- c(
  "---",
  "title: \"Medicare OLS Regression Results\"",
  "output: pdf_document",
  "header-includes:",
  "  - \"\\usepackage{booktabs}\"",
  "  - \"\\usepackage{longtable}\"",
  "  - \"\\usepackage{float}\"",
  "---",
  "",
  "```{r setup, include=FALSE}",
  "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)",
  "library(knitr)",
  "library(kableExtra)",
  "```",
  "",
  "# OLS Regression Tables",
  "",
  "## Education outcome regressions",
  "",
  "```{r edu_coeffs, results='asis'}",
  "print(format_table_for_pdf(education_coeffs, \"Education outcome OLS coefficients\"))",
  "```",
  "",
  "```{r edu_fit, results='asis'}",
  "print(format_table_for_pdf(education_fit, \"Education outcome OLS fit statistics\"))",
  "```",
  "",
  "## Education first-stage regressions",
  "",
  "```{r edu_fs_coeffs, results='asis'}",
  "print(format_table_for_pdf(education_fs_coeffs, \"Education first-stage OLS coefficients\"))",
  "```",
  "",
  "```{r edu_fs_fit, results='asis'}",
  "print(format_table_for_pdf(education_fs_fit, \"Education first-stage OLS fit statistics\"))",
  "```",
  "",
  "## Manual RD robustness regressions",
  "",
  "```{r manual_coeffs, results='asis'}",
  "print(format_table_for_pdf(manual_rd_coeffs, \"Manual RD OLS coefficients\"))",
  "```",
  "",
  "```{r manual_fit, results='asis'}",
  "print(format_table_for_pdf(manual_rd_fit, \"Manual RD OLS fit statistics\"))",
  "```",
  ""
)

writeLines(report_lines, con = report_path)

export_success <- FALSE
tryCatch({
  rmarkdown::render(
    input = report_path,
    output_file = "ols_regression_tables_report.pdf",
    output_dir = output_dir,
    envir = new.env(parent = globalenv()),
    quiet = TRUE
  )
  export_success <- TRUE
}, error = function(e) {
  warning("PDF export failed via rmarkdown; writing a standalone LaTeX file instead.")
  export_success <- FALSE
})

if (!export_success) {
  tex_path <- file.path(output_dir, "ols_regression_tables_report.tex")
  tex_lines <- c(
    "\\documentclass{article}",
    "\\usepackage[margin=1in]{geometry}",
    "\\usepackage{booktabs}",
    "\\usepackage{longtable}",
    "\\usepackage{float}",
    "\\begin{document}",
    "\\section*{OLS Regression Tables}",
    "\\subsection*{Education outcome regressions}",
    knitr::kable(pretty_ols_table(education_coeffs), format = "latex", digits = 4, caption = "Education outcome OLS coefficients", booktabs = TRUE, escape = FALSE),
    knitr::kable(pretty_ols_table(education_fit), format = "latex", digits = 4, caption = "Education outcome OLS fit statistics", booktabs = TRUE, escape = FALSE),
    "\\subsection*{Education first-stage regressions}",
    knitr::kable(pretty_ols_table(education_fs_coeffs), format = "latex", digits = 4, caption = "Education first-stage OLS coefficients", booktabs = TRUE, escape = FALSE),
    knitr::kable(pretty_ols_table(education_fs_fit), format = "latex", digits = 4, caption = "Education first-stage OLS fit statistics", booktabs = TRUE, escape = FALSE),
    "\\subsection*{Manual RD robustness regressions}",
    knitr::kable(pretty_ols_table(manual_rd_coeffs), format = "latex", digits = 4, caption = "Manual RD OLS coefficients", booktabs = TRUE, escape = FALSE),
    knitr::kable(pretty_ols_table(manual_rd_fit), format = "latex", digits = 4, caption = "Manual RD OLS fit statistics", booktabs = TRUE, escape = FALSE),
    "\\end{document}",
    ""
  )
  writeLines(tex_lines, con = tex_path)
  cat("\nLaTeX export saved to:", tex_path, "\n")
}

# -----------------------------------------------------------------------------
# Legacy CSV export kept as commented-out reference for ad hoc data work
# -----------------------------------------------------------------------------
# output_dir <- "output"
# dir.create(output_dir, showWarnings = FALSE)
# write.csv(education_coeffs, file = file.path(output_dir, "education_ols_coefficients.csv"), row.names = FALSE)
# write.csv(education_fit, file = file.path(output_dir, "education_ols_fit_statistics.csv"), row.names = FALSE)
# write.csv(education_fs_coeffs, file = file.path(output_dir, "education_first_stage_ols_coefficients.csv"), row.names = FALSE)
# write.csv(education_fs_fit, file = file.path(output_dir, "education_first_stage_ols_fit_statistics.csv"), row.names = FALSE)
# write.csv(manual_rd_coeffs, file = file.path(output_dir, "manual_rd_ols_coefficients.csv"), row.names = FALSE)
# write.csv(manual_rd_fit, file = file.path(output_dir, "manual_rd_ols_fit_statistics.csv"), row.names = FALSE)

cat("\nPublication-ready regression tables have been exported to the output/ directory.\n")
