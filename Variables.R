# ==============================================================================
# construct_variables.R
#
# Purpose: Reconstructs the analysis variables used in "How Does Medicare
# Interact with Educational Attainment? Evidence from a Regression
# Discontinuity Design" (Mai, Datta, Lee) from a raw IPUMS Health Surveys
# (NHIS) extract.
#
# This script does NOT include or distribute IPUMS microdata. IPUMS terms of
# use prohibit redistributing extracts. To use this script:
#
#   1. Create a free IPUMS Health Surveys account: https://nhis.ipums.org
#   2. Build an extract for NHIS survey years 2015-2024 with the variables
#      listed in the "REQUIRED IPUMS VARIABLES" block below.
#   3. Download the extract as .csv (or read a .dta via haven::read_dta) and
#      point RAW_DATA_PATH at it.
#   4. Run this script to produce the analysis-ready data frame.
#
# NOTE ON COLUMN NAME CASE: this script lowercases all column names
# immediately on load (see step 0), and every reference below is lowercase.
# That means it works the same way whether the input is a freshly pulled
# IPUMS extract (which uses uppercase names like AGE, HINOTCOVE, EDUC) or
# an already-processed file that went through this script before (which
# will already be lowercase). You never need to hand-edit case again.
#
# IMPORTANT - VARIABLE AVAILABILITY BY YEAR:
# Not every variable below is fielded in every NHIS year from 2015-2024.
# Per the IPUMS NHIS codebook (checked directly, Aug 2026):
#   - hinotcove, health, vacflu12m, empstat, educ: available every year
#     2015-2024.
#   - premedyr, dvintwell: available 2019-2024 only (not fielded 2015-2018).
#   - cholchekyrg, hypcheckgmo: available ONLY in 2019, 2021, and 2023 within
#     the 2015-2024 window (also fielded in some pre-2015 years, but those
#     are outside this paper's sample period).
# This means the cholesterol-check and blood-pressure-check outcomes in
# particular have a much smaller effective sample (three survey years, not
# ten) than insurance/health/flu-shot. Confirm this against your own
# extract's year coverage before comparing Ns to the paper's Table A1.
# ==============================================================================

library(tidyverse)

# ------------------------------------------------------------------------
# 0. SETUP
# ------------------------------------------------------------------------
# Load your extract ONE of these ways -- use whichever matches what you
# actually have, and leave the other two commented out. Do NOT run more
# than one of these blocks: each one reassigns `raw`, so running a second
# block after the first will silently overwrite/corrupt the correctly
# loaded data with whatever the second one produces (this is exactly what
# happened when a leftover read_csv() line ran a .dat fixed-width file
# through a comma-delimited parser and clobbered a working `raw` object).

# --- Option A: load directly via ipumsr (DDI + data file) -- RECOMMENDED,
# this is what actually works for the Sample Medicare extract ---
library(ipumsr)
file_path <- "nhis_00007.xml" # edit to your configuration
ddi <- read_ipums_ddi("nhis_00007.xml")
raw <- read_ipums_micro(ddi)

# --- Option B: load from a plain CSV (only if you exported one) ---
# RAW_DATA_PATH <- "nhis_2015_2024_extract.csv"
# raw <- read_csv(RAW_DATA_PATH)

# --- Option C: load from a Stata .dta file ---
# RAW_DATA_PATH <- "finalfile1.dta"
# raw <- haven::read_dta(RAW_DATA_PATH)

# Standardize all column names to lowercase immediately, before anything
# else touches this data. This is what makes the rest of the script
# case-agnostic -- see the note at the top of the file.
#
# Safety check first: if lowercasing would collide two distinct existing
# columns into the same name, do NOT silently overwrite one -- stop and
# surface it instead.
dupe_check <- tolower(names(raw))
if (any(duplicated(dupe_check))) {
  clashing <- names(raw)[dupe_check %in% dupe_check[duplicated(dupe_check)]]
  stop(
    "Lowercasing column names would collide the following existing columns ",
    "into duplicates -- resolve manually before proceeding: ",
    paste(clashing, collapse = ", ")
  )
}
names(raw) <- tolower(names(raw))

# ------------------------------------------------------------------------
# REQUIRED IPUMS VARIABLES (request these when building your extract)
# ------------------------------------------------------------------------
# year              survey year
# serial, pernum    household/person identifiers (for merging/weights only)
# sampweight        NHIS annual sample weight (person-level)
# strata, psu       design variables (only needed if replicating with
#                   survey-weighted SEs rather than the RD package's own)
# age               age in years (running variable)
# educ              detailed educational attainment recode
# hinotcove         insurance coverage status
# health            self-reported general health status (5-point scale)
# vacflu12m         had any flu vaccine, past 12 months
# cholchekyrg       duration since last cholesterol check (grouped years)
# hypcheckgmo       duration since last blood pressure check (grouped months)
# dvintwell         interval since last wellness visit
# premedyr          used prescription medication, past 12 months
# empstat           employment status in past 1-2 weeks
# astatflg          Sample Adult flag -- needed to restrict to the Sample
#                   Adult subsample, since cholchekyrg/hypcheckgmo/
#                   dvintwell/premedyr are only asked of Sample Adults

#### 1. RUNNING VARIABLE AND TREATMENT INDICATOR ####

df <- raw %>%
  mutate(
    age         = age,
    age_c       = age - 65,             # age centered at the RD cutoff
    age65       = as.integer(age >= 65),# treatment/eligibility indicator
    age_x_treat = age_c * age65         # RD interaction term
  )


#### 2. INSURANCE COVERAGE ####
# hinotcove codes: 0 NIU, 1 "No, has coverage" (i.e. covered),
# 2 "Yes, has no coverage" (i.e. NOT covered), 7/8/9 unknown.

df <- df %>%
  mutate(
    insured = case_when(
      hinotcove == 1 ~ 1L,   # covered
      hinotcove == 2 ~ 0L,   # not covered
      TRUE ~ NA_integer_     # 0 NIU, 7/8/9 unknown
    )
  )

#### 3. SELF-REPORTED HEALTH (monotone 1-5 scale, higher = better) ####
# health codes: 1 Excellent, 2 Very Good, 3 Good, 4 Fair, 5 Poor, 7/8/9
# unknown. Raw coding already has 1 = best; the paper describes Medicare as
# having a *positive* effect on this "monotone scale," so we reverse it so
# higher values = better self-reported health (5 = excellent, 1 = poor),
# consistent with Table A1's reported means (e.g. 3.79 for Bachelor's+).

df <- df %>%
  mutate(
    health_r = case_when(
      health %in% 1:5 ~ 6L - health,
      TRUE ~ NA_integer_    # 7/8/9 unknown
    )
  )

##### 4. PREVENTIVE CARE PROXIES ####

#### Flu shot ####
# vacflu12m codes: 0 NIU, 1 No, 2 Yes, 3 "both shot and spray" (only
# available pre-2015), 7/8/9 unknown.
df <- df %>%
  mutate(
    flu_shot = case_when(
      vacflu12m == 2 ~ 1L,
      vacflu12m == 1 ~ 0L,
      TRUE ~ NA_integer_   # 0 NIU, 7/8/9 unknown
    )
  )

#### Cholesterol check ####
# cholchekyrg codes (duration since last check, grouped years):
#   000 NIU, 100 Never,
#   200 "1 year or less", 210 "Less than 1 year ago"      -> within past yr
#   300 "1 to 2 yrs", 310 "1 yr to <2 yrs", 320 ">1 yr to 2 yrs"  -> 1-2 yrs
#   400+ "2+ years" and further breakdowns                -> 2+ yrs
#   997/998/999 unknown
# Binary "cholesterol check within the past 12 months":
df <- df %>%
  mutate(
    cholcheck_1yr = case_when(
      cholchekyrg %in% c(200, 210) ~ 1L,
      cholchekyrg %in% c(100, 300, 310, 320, 400, 410, 420, 430, 440, 450,
                         460, 470, 480, 481, 482) ~ 0L,
      TRUE ~ NA_integer_   # 000 NIU, 997/998/999 unknown
    )
  )

#### Blood pressure check ####
# hypcheckgmo codes (duration since last check, grouped months):
#   100-140 = various "within past year" groupings (1yr or less, 3mo or
#             less, <1mo, 1-3mo, 4-6mo, 7-12mo, within past yr <12mo)
#   200-220 = 1-2 years ago
#   300-330 = 2-5 years ago
#   400-420 = 5+ years ago
#   500     = Never
#   996 NIU, 997/998/999 unknown
# Binary "blood pressure check within the past 12 months":
df <- df %>%
  mutate(
    bpcheck_1yr = case_when(
      hypcheckgmo >= 100 & hypcheckgmo < 200 ~ 1L,
      hypcheckgmo >= 200 & hypcheckgmo < 996 ~ 0L,
      hypcheckgmo == 500 ~ 0L,
      TRUE ~ NA_integer_   # 996 NIU, 997/998/999 unknown
    )
  )

#### Wellness visit ####
# dvintwell codes: 00 Never, 01 "Within the past year", 02 "1-2 yrs ago",
# 03 "2-3 yrs ago", 04 "3-5 yrs ago", 05 "5-10 yrs ago", 06 "10+ yrs ago",
# 96 NIU, 97/98/99 unknown.
df <- df %>%
  mutate(
    wellvisit_1yr = case_when(
      dvintwell == 1 ~ 1L,
      dvintwell %in% c(0, 2, 3, 4, 5, 6) ~ 0L,
      TRUE ~ NA_integer_   # 96 NIU, 97/98/99 unknown
    )
  )

#### Prescription medication ####
# premedyr codes: 0 NIU, 1 No, 2 Yes, 7/8/9 unknown.
df <- df %>%
  mutate(
    presc_med = case_when(
      premedyr == 2 ~ 1L,
      premedyr == 1 ~ 0L,
      TRUE ~ NA_integer_   # 0 NIU, 7/8/9 unknown
    )
  )

#### 5. EMPLOYMENT (for the retirement/continuity check, Table 1) ####
# empstat codes: 000 NIU, 100-122 = employed (various work-status detail),
# 200-220 = not employed (unemployed / not in labor force),
# 900/997/998/999 unknown.
df <- df %>%
  mutate(
    employed = case_when(
      empstat >= 100 & empstat < 200 ~ 1L,
      empstat >= 200 & empstat < 900 ~ 0L,
      TRUE ~ NA_integer_   # 000 NIU, 900/997/998/999 unknown
    )
  )

#### 6. EDUCATIONAL ATTAINMENT GROUPS (heterogeneity dimension) ####
# educ codes (confirmed from IPUMS NHIS codebook):
#   100-116  = less than HS (Grade 12 or less, no diploma; grades 1-11;
#              never attended; 12th grade no diploma)
#   201-202  = HS diploma or GED
#   301-303  = Some college, no 4yr degree / AA degree
#   400      = Bachelor's degree
#   510-530  = Master's / Professional / Doctoral / Other degree
#   000      = NIU
#   996-999  = unknown / no degree-years unknown
# The paper's three buckets: No HS (<200), HS/Some College (200-399,
# reference group), Bachelor's+ (400+).

df <- df %>%
  mutate(
    edu3 = case_when(
      educ >= 100 & educ < 200 ~ "No HS",
      educ >= 201 & educ < 400 ~ "HS/Some College",
      educ >= 400 & educ < 996 ~ "Bachelors+",
      TRUE ~ NA_character_   # 000 NIU, 996-999 unknown
    )
  )

# ------------------------------------------------------------------------
# 7. FINAL ANALYSIS SAMPLE
# ------------------------------------------------------------------------
# Paper restricts to adults with non-missing age; RD estimation itself
# further restricts to each outcome's CCT-optimal bandwidth window inside
# the rdrobust/manual RD calls (see rd_estimation.R), not here.
#
# NOTE: because cholchekyrg/hypcheckgmo are only fielded in 2019/2021/2023
# and premedyr/dvintwell only from 2019 onward, rows from 2015-2018 will
# have NA on those specific outcomes by construction (not missing data in
# the usual sense) -- this is expected, not a data error.
#
# This keeps ALL columns -- every raw IPUMS variable from your extract
# plus every newly constructed variable (age_c, age65, age_x_treat,
# insured, health_r, flu_shot, cholcheck_1yr, bpcheck_1yr, wellvisit_1yr,
# presc_med, employed, edu3). If you only want the constructed variables
# plus identifiers, swap this for an explicit select() whitelist.

analysis_df <- df %>%
  filter(!is.na(age))

write_csv(analysis_df, "analysis_sample_full.csv")
# To write a Stata file instead:
# haven::write_dta(analysis_df, "analysis_sample_full.dta")

# ------------------------------------------------------------------------
# Sanity checks against Table A1 / summary statistics reported in the paper
# ------------------------------------------------------------------------
# nrow(analysis_df)                                   # expect ~301,290 pre-bandwidth
# mean(analysis_df$insured, na.rm = TRUE)              # sanity check vs Table A1
# table(analysis_df$edu3)
# analysis_df %>% filter(!is.na(cholcheck_1yr)) %>% count(year)  # should be 2019/2021/2023 only