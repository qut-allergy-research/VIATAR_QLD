# =============================================================================
# VIATAR_QLD
# Supplementary Table 2. Age-specific prevalence (%) of aeroallergen
#                        sensitisation
#
# PURPOSE
# Reproduce the age-specific sensitisation table reported in the supplementary
# material. Sensitisation prevalence is calculated separately for five age
# groups for each selected aeroallergen.
#
# AGE GROUPS
#   12–17 years
#   18–29 years
#   30–44 years
#   45–64 years
#   >=65 years
#
# SENSITISATION DEFINITION
#   sIgE >= 0.35 kU/L
#
# ANALYTICAL UNIT
# This script uses individual sIgE TEST RECORDS as observations, consistent
# with the current manuscript Methods and with the original working code.
# Individuals may therefore contribute multiple test records.
#
# STATISTICAL COMPARISON
# For each allergen, a 5 x 2 chi-square test compares the proportion of
# sensitised versus non-sensitised test records across the five age groups.
#
# MULTIPLE TESTING
# Benjamini-Hochberg (BH/FDR) correction is applied across the 21 allergen-level
# omnibus age-group tests displayed in this supplementary table.
#
# DATA PRIVACY
# The real pathology dataset is NOT included in the public repository.
#
# IMPORTANT
# Earlier working code calculated age-group summaries stratified by sex.
# This cleaned script removes the sex stratification because the current
# supplementary table reports overall age-specific prevalence.
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Packages
# -----------------------------------------------------------------------------

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(purrr)


# -----------------------------------------------------------------------------
# 2. Local paths
# -----------------------------------------------------------------------------

data_path <- file.path(
  "data",
  "VIATAR_filtered_QLD_01122025.xlsx"
)

output_dir <- file.path(
  "outputs",
  "supplementary_table2"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(data_path)) {
  stop(
    "Data file not found: ", data_path,
    "\nUpdate 'data_path' locally. Do not upload the real dataset to GitHub."
  )
}


# -----------------------------------------------------------------------------
# 3. Read data
# -----------------------------------------------------------------------------

df_raw <- read_excel(
  data_path,
  sheet = 1
)


# -----------------------------------------------------------------------------
# 4. Check required columns
# -----------------------------------------------------------------------------

required_columns <- c(
  "Age.in.Years",
  "Test.Name..Allergen.",
  "test.Value..SPT.or.spIgE.values.of.each.allergen.or.mmD.for.SPT"
)

missing_columns <- setdiff(
  required_columns,
  names(df_raw)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing required column(s): ",
    paste(missing_columns, collapse = ", ")
  )
}


# -----------------------------------------------------------------------------
# 5. Supplementary Table 2 allergen lookup
# -----------------------------------------------------------------------------
#
# The table contains 21 individual aeroallergens.
# Several raw spelling variants are accepted for a small number of allergens.
# -----------------------------------------------------------------------------

allergen_lookup <- tibble::tribble(
  ~order, ~allergen_group,              ~allergen_label,                         ~raw_key,

   1L,    "House dust mite (HDM)",      "Dermatophagoides farinae",              "dermatophagoides farinae",
   2L,    "House dust mite (HDM)",      "Dermatophagoides pteronyssinus",        "dermatophagoides pteronyssinus",
   3L,    "House dust mite (HDM)",      "Blomia tropicalis",                     "blomia tropicalis",

   4L,    "Moulds",                     "Alternaria alternata",                  "alternaria alternata",
   5L,    "Moulds",                     "Aspergillus fumigatus",                 "aspergillus fumigatus",
   6L,    "Moulds",                     "Cladosporium herbarum",                 "cladosporium herbarum",
   7L,    "Moulds",                     "Mould Mix 1",                           "mould mix 1",

   8L,    "Grass pollen",               "Bahia grass (Paspalum)",                "bahia grass (paspalum)",
   9L,    "Grass pollen",               "Bermuda (couch) grass",                 "bermuda (couch) grass",
  10L,    "Grass pollen",               "Johnson grass (Sorghum)",               "johnson grass (sorghum)",
  11L,    "Grass pollen",               "Perennial rye grass",                   "perennial rye grass",
  12L,    "Grass pollen",               "Grass Mix 2",                           "grass mix 2",

  13L,    "Weed pollen",                "Common ragweed",                        "common ragweed",
  14L,    "Weed pollen",                "Weed Mix 1",                            "weed mix 1",

  15L,    "Tree pollen",                "Eucalyptus (gum tree)",                 "eucalyptus, gum tree",
  15L,    "Tree pollen",                "Eucalyptus (gum tree)",                 "eucalyptus (gum tree)",
  15L,    "Tree pollen",                "Eucalyptus (gum tree)",                 "eucalyptus",
  16L,    "Tree pollen",                "Acacia",                                "acacia",
  17L,    "Tree pollen",                "Tree Mix 7",                            "tree mix 7",

  18L,    "Animal dander",              "Cat epithelium & dander",               "cat epithelium & dander",
  18L,    "Animal dander",              "Cat epithelium & dander",               "cat epithelium and dander",
  19L,    "Animal dander",              "Dog dander",                            "dog dander",
  20L,    "Animal dander",              "Animal Mix 1",                          "animal mix 1",

  21L,    "Cockroach",                  "Blattella germanica",                   "cockroach (blatella germanica)",
  21L,    "Cockroach",                  "Blattella germanica",                   "cockroach (blattella germanica)",
  21L,    "Cockroach",                  "Blattella germanica",                   "blattella germanica"
)


# -----------------------------------------------------------------------------
# 6. Prepare age, allergen, and sIgE variables
# -----------------------------------------------------------------------------
#
# IMPORTANT VALUE HANDLING
# parse_number("<0.35") returns 0.35. We therefore retain a below-limit flag so
# "<0.35" is correctly classified as negative.
# -----------------------------------------------------------------------------

analysis_data <- df_raw %>%
  transmute(
    age = suppressWarnings(
      as.numeric(`Age.in.Years`)
    ),

    allergen_raw = as.character(
      `Test.Name..Allergen.`
    ),

    sige_raw = as.character(
      `test.Value..SPT.or.spIgE.values.of.each.allergen.or.mmD.for.SPT`
    )
  ) %>%
  mutate(
    allergen_key = allergen_raw %>%
      str_to_lower() %>%
      str_squish(),

    sige_text = str_trim(
      sige_raw
    ),

    below_limit = str_detect(
      sige_text,
      "^<"
    ),

    sige = parse_number(
      sige_text
    ),

    age_group = case_when(
      age >= 12 & age <= 17 ~ "12–17 yrs",
      age >= 18 & age <= 29 ~ "18–29 yrs",
      age >= 30 & age <= 44 ~ "30–44 yrs",
      age >= 45 & age <= 64 ~ "45–64 yrs",
      age >= 65             ~ "≥65 yrs",
      TRUE                  ~ NA_character_
    )
  ) %>%
  filter(
    !is.na(age_group),
    !is.na(allergen_key),
    !is.na(sige)
  ) %>%
  left_join(
    allergen_lookup,
    by = c(
      "allergen_key" = "raw_key"
    )
  ) %>%
  filter(
    !is.na(allergen_label)
  ) %>%
  mutate(
    sensitised = case_when(
      below_limit ~ 0L,
      sige >= 0.35 ~ 1L,
      TRUE ~ 0L
    ),

    age_group = factor(
      age_group,
      levels = c(
        "12–17 yrs",
        "18–29 yrs",
        "30–44 yrs",
        "45–64 yrs",
        "≥65 yrs"
      )
    )
  )


# -----------------------------------------------------------------------------
# 7. Audit exact raw allergen labels included in the analysis
# -----------------------------------------------------------------------------

mapping_audit <- analysis_data %>%
  count(
    order,
    allergen_group,
    allergen_label,
    allergen_raw,
    name = "n_test_records"
  ) %>%
  arrange(
    order,
    desc(n_test_records)
  )

write_csv(
  mapping_audit,
  file.path(
    output_dir,
    "SupplementaryTable2_allergen_mapping_audit.csv"
  )
)


# -----------------------------------------------------------------------------
# 8. Confirm all 21 table allergens are present
# -----------------------------------------------------------------------------

expected_targets <- allergen_lookup %>%
  distinct(
    order,
    allergen_group,
    allergen_label
  )

observed_targets <- analysis_data %>%
  distinct(
    order,
    allergen_group,
    allergen_label
  )

missing_targets <- expected_targets %>%
  anti_join(
    observed_targets,
    by = c(
      "order",
      "allergen_group",
      "allergen_label"
    )
  )

if (nrow(missing_targets) > 0) {
  print(missing_targets)

  stop(
    "One or more Supplementary Table 2 allergens were not found. ",
    "Check the allergen mapping and raw labels."
  )
}


# -----------------------------------------------------------------------------
# 9. Age-specific prevalence for each allergen
# -----------------------------------------------------------------------------

age_summary_long <- analysis_data %>%
  group_by(
    order,
    allergen_group,
    allergen_label,
    age_group
  ) %>%
  summarise(
    n_tested = n(),

    n_sensitised = sum(
      sensitised == 1L,
      na.rm = TRUE
    ),

    prevalence = n_sensitised / n_tested,

    prevalence_pct = 100 * prevalence,

    .groups = "drop"
  ) %>%
  arrange(
    order,
    age_group
  )


# -----------------------------------------------------------------------------
# 10. Omnibus chi-square test across the five age groups
# -----------------------------------------------------------------------------
#
# For each allergen:
# rows    = 5 age groups
# columns = sensitised / not sensitised
#
# The current supplementary-table footnote specifies chi-square tests.
# We therefore use Pearson's chi-square test without continuity correction
# (continuity correction is not applicable to a 5 x 2 table).
# -----------------------------------------------------------------------------

age_chisq_one_allergen <- function(dat) {

  age_counts <- dat %>%
    group_by(age_group) %>%
    summarise(
      positive = sum(
        sensitised == 1L,
        na.rm = TRUE
      ),

      negative = sum(
        sensitised == 0L,
        na.rm = TRUE
      ),

      .groups = "drop"
    ) %>%
    complete(
      age_group = factor(
        c(
          "12–17 yrs",
          "18–29 yrs",
          "30–44 yrs",
          "45–64 yrs",
          "≥65 yrs"
        ),
        levels = c(
          "12–17 yrs",
          "18–29 yrs",
          "30–44 yrs",
          "45–64 yrs",
          "≥65 yrs"
        )
      ),
      fill = list(
        positive = 0,
        negative = 0
      )
    ) %>%
    arrange(age_group)

  tab <- as.matrix(
    age_counts %>%
      select(
        positive,
        negative
      )
  )

  rownames(tab) <- as.character(
    age_counts$age_group
  )

  # All five age groups should have data for the selected allergens.
  if (any(rowSums(tab) == 0)) {
    return(
      tibble(
        chi_square = NA_real_,
        df = NA_real_,
        p_value = NA_real_,
        min_expected = NA_real_,
        note = "At least one age group had no usable records"
      )
    )
  }

  test <- suppressWarnings(
    chisq.test(
      tab,
      correct = FALSE
    )
  )

  tibble(
    chi_square = unname(
      test$statistic
    ),

    df = unname(
      test$parameter
    ),

    p_value = test$p.value,

    min_expected = min(
      test$expected
    ),

    note = if_else(
      min(test$expected) < 5,
      "Warning: expected count <5",
      ""
    )
  )
}


age_tests <- analysis_data %>%
  group_by(
    order,
    allergen_group,
    allergen_label
  ) %>%
  group_modify(
    ~ age_chisq_one_allergen(.x)
  ) %>%
  ungroup() %>%
  arrange(order) %>%
  mutate(
    p_adj_BH = p.adjust(
      p_value,
      method = "BH"
    )
  )


# -----------------------------------------------------------------------------
# 11. Save numerical audit tables
# -----------------------------------------------------------------------------

write_csv(
  age_summary_long,
  file.path(
    output_dir,
    "SupplementaryTable2_age_prevalence_numeric.csv"
  )
)

write_csv(
  age_tests,
  file.path(
    output_dir,
    "SupplementaryTable2_age_chisquare_tests.csv"
  )
)


# -----------------------------------------------------------------------------
# 12. Reshape prevalences into manuscript format
# -----------------------------------------------------------------------------

age_prevalence_wide <- age_summary_long %>%
  select(
    order,
    allergen_group,
    allergen_label,
    age_group,
    prevalence_pct
  ) %>%
  pivot_wider(
    names_from = age_group,
    values_from = prevalence_pct
  )


# -----------------------------------------------------------------------------
# 13. Build manuscript-ready Supplementary Table 2
# -----------------------------------------------------------------------------

supp_table2 <- age_prevalence_wide %>%
  left_join(
    age_tests %>%
      select(
        order,
        allergen_group,
        allergen_label,
        p_value,
        p_adj_BH
      ),
    by = c(
      "order",
      "allergen_group",
      "allergen_label"
    )
  ) %>%
  arrange(order) %>%
  mutate(
    across(
      c(
        `12–17 yrs`,
        `18–29 yrs`,
        `30–44 yrs`,
        `45–64 yrs`,
        `≥65 yrs`
      ),
      ~ round(.x, 2)
    ),

    `p-value` = case_when(
      is.na(p_value) ~ NA_character_,
      p_value < 0.001 ~ "<0.001",
      TRUE ~ sprintf(
        "%.3f",
        p_value
      )
    ),

    `Adjusted P-value (FDR)` = case_when(
      is.na(p_adj_BH) ~ NA_character_,
      p_adj_BH < 0.001 ~ "<0.001",
      TRUE ~ sprintf(
        "%.3f",
        p_adj_BH
      )
    )
  ) %>%
  transmute(
    `Allergen group` = allergen_group,
    `Allergen` = allergen_label,
    `12–17 yrs` = .data[["12–17 yrs"]],
    `18–29 yrs` = .data[["18–29 yrs"]],
    `30–44 yrs` = .data[["30–44 yrs"]],
    `45–64 yrs` = .data[["45–64 yrs"]],
    `≥65 yrs` = .data[["≥65 yrs"]],
    `p-value` = .data[["p-value"]],
    `Adjusted P-value (FDR)` =
      .data[["Adjusted P-value (FDR)"]]
  )


# -----------------------------------------------------------------------------
# 14. Print and save final table
# -----------------------------------------------------------------------------

print(
  supp_table2,
  n = Inf
)

write_csv(
  supp_table2,
  file.path(
    output_dir,
    "SupplementaryTable2_age_specific_sensitisation.csv"
  )
)


# -----------------------------------------------------------------------------
# 15. Optional validation against the percentages currently reported
# -----------------------------------------------------------------------------
#
# These values are used ONLY as a reproducibility check.
# They are NOT used to calculate the analysis.
# -----------------------------------------------------------------------------

reported_values <- tibble::tribble(
  ~allergen_label,                        ~`12–17 yrs`, ~`18–29 yrs`, ~`30–44 yrs`, ~`45–64 yrs`, ~`≥65 yrs`,

  "Dermatophagoides farinae",                  57.46,        50.96,        46.06,        32.86,      25.25,
  "Dermatophagoides pteronyssinus",            50.38,        52.01,        49.01,        38.19,      26.84,
  "Blomia tropicalis",                         52.34,        43.20,        32.87,        23.34,      13.56,

  "Alternaria alternata",                      20.86,        12.35,        12.77,        10.29,       6.60,
  "Aspergillus fumigatus",                     12.35,         9.04,         6.92,        10.34,      17.02,
  "Cladosporium herbarum",                     11.24,         4.02,         2.48,         2.23,       1.50,
  "Mould Mix 1",                               14.42,        12.16,         9.86,         9.91,       9.39,

  "Bahia grass (Paspalum)",                    23.84,        23.19,        22.80,        22.40,      13.81,
  "Bermuda (couch) grass",                     23.52,        24.59,        23.59,        21.14,      13.66,
  "Johnson grass (Sorghum)",                   21.53,        21.07,        20.70,        19.51,      13.02,
  "Perennial rye grass",                       22.36,        21.42,        23.51,        20.64,      13.29,
  "Grass Mix 2",                               22.52,        23.01,        22.04,        20.94,      14.73,

  "Common ragweed",                            17.75,        12.35,        10.13,         9.08,       6.30,
  "Weed Mix 1",                                17.82,        14.50,        11.02,         9.38,       8.71,

  "Eucalyptus (gum tree)",                      6.92,         4.12,         5.10,         4.35,       3.78,
  "Acacia",                                     8.52,         4.26,         5.34,         4.43,       4.05,
  "Tree Mix 7",                                16.70,        12.66,         9.32,         7.69,       6.38,

  "Cat epithelium & dander",                   29.67,        31.17,        21.71,        14.73,       8.51,
  "Dog dander",                                28.46,        28.40,        19.16,        14.16,       7.38,
  "Animal Mix 1",                              26.40,        25.42,        17.85,        12.94,       7.98,

  "Blattella germanica",                       24.92,        21.18,        14.86,        14.18,       9.19
)


validation <- age_prevalence_wide %>%
  select(
    allergen_label,
    `12–17 yrs`,
    `18–29 yrs`,
    `30–44 yrs`,
    `45–64 yrs`,
    `≥65 yrs`
  ) %>%
  left_join(
    reported_values,
    by = "allergen_label",
    suffix = c("_calculated", "_reported")
  ) %>%
  mutate(
    diff_12_17 =
      `12–17 yrs_calculated` -
      `12–17 yrs_reported`,

    diff_18_29 =
      `18–29 yrs_calculated` -
      `18–29 yrs_reported`,

    diff_30_44 =
      `30–44 yrs_calculated` -
      `30–44 yrs_reported`,

    diff_45_64 =
      `45–64 yrs_calculated` -
      `45–64 yrs_reported`,

    diff_65plus =
      `≥65 yrs_calculated` -
      `≥65 yrs_reported`,

    max_abs_difference_pp = pmax(
      abs(diff_12_17),
      abs(diff_18_29),
      abs(diff_30_44),
      abs(diff_45_64),
      abs(diff_65plus),
      na.rm = TRUE
    )
  )

write_csv(
  validation,
  file.path(
    output_dir,
    "SupplementaryTable2_validation_against_reported_values.csv"
  )
)

if (
  any(
    validation$max_abs_difference_pp > 0.05,
    na.rm = TRUE
  )
) {
  warning(
    "At least one calculated age-specific prevalence differs from the ",
    "currently reported table by >0.05 percentage points. Inspect the ",
    "validation file before finalising the manuscript."
  )
}


# -----------------------------------------------------------------------------
# 16. Session information
# -----------------------------------------------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    output_dir,
    "SupplementaryTable2_sessionInfo.txt"
  )
)
