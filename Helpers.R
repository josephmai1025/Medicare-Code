### A set of functions used across the various analysis scripts.
### These are mostly for managing packages and files.
### Inserted here for better code modularity.


# A function to install a package if it hasn't been already. 
# Then load the library.
ensure_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

# A function to find the given data file anywhere under 
# the project directory tree.
find_data_file <- function(filename = "finalfile1.dta", start_dir = getwd()) {
  current_dir <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)
  repeat {
    file_pattern <- paste0("^", basename(filename), "$")
    matches <- list.files(
      path = current_dir,
      pattern = file_pattern,
      recursive = TRUE,
      full.names = TRUE
    )
    if (length(matches) > 0) {
      return(matches[1])
    }
    parent_dir <- dirname(current_dir)
    if (identical(parent_dir, current_dir)) {
      break
    }
    current_dir <- parent_dir
  }
  stop(paste0("Could not find ", filename, ". Expected it somewhere under the project folder."))
}

# A function for running a fuzzy RD regression
run_fuzzy_rd <- function(y_var, var_name, fuzzy_var) {
  rd <- rdrobust(
    y = y_var,
    x = data$age,
    c = 65,
    p = 1,
    bwselect = "mserd",
    fuzzy = fuzzy_var,
    masspoints = "adjust",
    weights = data$sampweight
  )
  
  est <- rd$coef[1]
  se  <- rd$se[3]
  z <- est / se
  p_val <- 2 * (1 - pnorm(abs(z)))
  stars <- ifelse(p_val < 0.01, "***",
                  ifelse(p_val < 0.05, "**",
                         ifelse(p_val < 0.1, "*", "")))
  
  return(data.frame(
    Outcome = var_name,
    Bandwidth = round(rd$bws[1, 1], 2),
    Estimate = round(est, 3),
    SE = round(se, 3),
    Stars = stars
  ))
}