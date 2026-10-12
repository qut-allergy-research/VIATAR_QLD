# =============================================================================
# VIATAR_QLD
# Supplementary Table 1. Frequency of aeroallergen sensitisation by sex
#
# PURPOSE
# Reproduce the manuscript supplementary table comparing aeroallergen
# sensitisation frequency between female and male sIgE test records.
#
# ANALYTICAL UNIT
# This is a TEST-RECORD-LEVEL analysis.
# Each eligible sIgE test record is retained as one observation; individuals
# may contribute more than one test record.
#
# ELIGIBILITY
#   - Queensland analysis dataset
#   - age >= 12 years
#   - sex recorded as Female or Male
#   - selected aeroallergens shown in Supplementary Table 1
#
# SENSITISATION DEFINITION
#   sIgE >= 0.35 kU/L
#
# STATISTICAL COMPARISON
#   Female vs Male sensitisation frequency is compared separately for each
#   allergen using a 2 x 2 chi-square test with Yates continuity correction.
#   Fisher's exact test is used if any expected cell count is < 5.
#
# MULTIPLE TESTING
#   Benjamini-Hochberg (BH/FDR) correction is applied across the allergen
#   comparisons displayed in this supplementary table.
#
# DATA PRIVACY
# The patient-level pathology dataset is NOT included in the public repository.
#
# IMPORTANT
# The original working file also contains AGE-STRATIFIED sex analyses.
# Those analyses are not used here because this supplementary table compares
# Female vs Male overall, not Female vs Male within age strata.
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
  "supplementary_table1"
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
  "Sex.at.birth",
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
# 5. Supplementary Table 1 allergen lookup
# -----------------------------------------------------------------------------
#
# raw_key:
#   normalised value expected in Test.Name..Allergen.
#
# allergen_label:
#   publication label used in the supplementary table
#
# Multiple raw spelling variants may point to the same publication label.
# -----------------------------------------------------------------------------

allergen_lookup <- tibble::tribble(
  ~order, ~allergen_group,      ~allergen_label,                         ~raw_key,

   1L,    "House dust mite",    "Dermatophagoides pteronyssinus",        "dermatophagoides pteronyssinus",
   2L,    "House dust mite",    "Dermatophagoides farinae",              "dermatophagoides farinae",
   3L,    "House dust mite",    "Blomia tropicalis",                     "blomia tropicalis",

   4L,    "Moulds",             "Alternaria alternata",                  "alternaria alternata",
   5L,    "Moulds",             "Aspergillus fumigatus",                 "aspergillus fumigatus",
   6L,    "Moulds",             "Cladosporium herbarum",                 "cladosporium herbarum",
   7L,    "Moulds",             "Mould Mix 1",                           "mould mix 1",

   8L,    "Grass pollen",       "Bahia grass (Paspalum)",                "bahia grass (paspalum)",
   9L,    "Grass pollen",       "Bermuda (couch) grass",                 "bermuda (couch) grass",
  10L,    "Grass pollen",       "Johnson grass (Sorghum)",               "johnson grass (sorghum)",
  11L,    "Grass pollen",       "Perennial rye grass",                   "perennial rye grass",
  12L,    "Grass pollen",       "Grass Mix 2",                           "grass mix 2",

  13L,    "Weed pollen",        "Common ragweed",                        "common ragweed",
  14L,    "Weed pollen",        "Weed Mix 1",                            "weed mix 1",

  15L,    "Tree pollen",        "Eucalyptus (gum tree)",                 "eucalyptus, gum tree",
  15L,    "Tree pollen",        "Eucalyptus (gum tree)",                 "eucalyptus (gum tree)",
  15L,    "Tree pollen",        "Eucalyptus (gum tree)",                 "eucalyptus",
  16L,    "Tree pollen",        "Acacia",                                "acacia",
  17L,    "Tree pollen",        "Tree Mix 7",                            "tree mix 7",

  18L,    "Animal dander",      "Cat epithelium & dander",               "cat epithelium & dander",
  18L,    "Animal dander",      "Cat epithelium & dander",               "cat epithelium and dander",
  19L,    "Animal dander",      "Dog dander",                            "dog dander",
  20L,    "Animal dander",      "Animal Mix 1",                          "animal mix 1",

  21L,    "Cockroach",          "Blattella germanica",                   "cockroach (blatella germanica)",
  21L,    "Cockroach",          "Blattella germanica",                   "cockroach (blattella germanica)",
  21L,    "Cockroach",          "Blattella germanica",                   "blattella germanica"
)


# -----------------------------------------------------------------------------
# 6. Prepare variables
# -----------------------------------------------------------------------------
#
# IMPORTANT:
# parse_number("<0.35") returns 0.35. Therefore the '<' symbol is retained
# separately so "<0.35" is correctly classified as NEGATIVE.
# -----------------------------------------------------------------------------

analysis_data <- df_raw %>%
  transmute(
    age = suppressWarnings(
      as.numeric(`Age.in.Years`)
    ),

    sex_raw = as.character(
      `Sex.at.birth`
    ),

    allergen_raw = as.character(
      `Test.Name..Allergen.`
    ),

    sige_raw = as.character(
      `test.Value..SPT.or.spIgE.values.of.each.allergen.or.mmD.for.SPT`
    )
  ) %>%
  mutate(
    sex = case_when(
      str_to_lower(str_trim(sex_raw)) %in% c(
        "female", "f", "woman", "girl"
      ) ~ "Female",

      str_to_lower(str_trim(sex_raw)) %in% c(
        "male", "m", "man", "boy"
      ) ~ "Male",

      TRUE ~ NA_character_
    ),

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
    )
  ) %>%
  filter(
    !is.na(age),
    age >= 12,
    sex %in% c("Female", "Male"),
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

    sex = factor(
      sex,
      levels = c(
        "Female",
        "Male"
      )
    )
  )


# -----------------------------------------------------------------------------
# 7. Audit exact raw allergen labels included in the table
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
    "SupplementaryTable1_allergen_mapping_audit.csv"
  )
)


# -----------------------------------------------------------------------------
# 8. Confirm all 21 publication allergens were found
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
    "One or more Supplementary Table 1 allergens were not found. ",
    "Check SupplementaryTable1_allergen_mapping_audit.csv and update the lookup."
  )
}


# -----------------------------------------------------------------------------
# 9. Count tested and sensitised test records by allergen and sex
# -----------------------------------------------------------------------------

summary_long <- analysis_data %>%
  group_by(
    order,
    allergen_group,
    allergen_label,
    sex
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
  )


# -----------------------------------------------------------------------------
# 10. Reshape Female and Male counts side-by-side
# -----------------------------------------------------------------------------

summary_wide <- summary_long %>%
  select(
    order,
    allergen_group,
    allergen_label,
    sex,
    n_tested,
    n_sensitised,
    prevalence_pct
  ) %>%
  pivot_wider(
    names_from = sex,
    values_from = c(
      n_tested,
      n_sensitised,
      prevalence_pct
    ),
    names_sep = "_"
  ) %>%
  arrange(order)


# -----------------------------------------------------------------------------
# 11. Female vs Male statistical test
# -----------------------------------------------------------------------------
#
# This follows the approach used in the working analysis code:
#   - chi-square test for ordinary 2 x 2 tables
#   - Fisher's exact test if any expected count is < 5
#
# R's standard 2 x 2 chisq.test applies Yates continuity correction by default.
# -----------------------------------------------------------------------------

compare_sexes <- function(
  female_tested,
  female_positive,
  male_tested,
  male_positive
) {

  tab <- matrix(
    c(
      female_positive,
      female_tested - female_positive,
      male_positive,
      male_tested - male_positive
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      Sex = c(
        "Female",
        "Male"
      ),
      Sensitisation = c(
        "Positive",
        "Negative"
      )
    )
  )

  chi <- suppressWarnings(
    chisq.test(
      tab,
      correct = TRUE
    )
  )

  if (any(chi$expected < 5)) {

    fisher <- fisher.test(
      tab
    )

    list(
      test_used = "Fisher exact",
      p_value = fisher$p.value
    )

  } else {

    list(
      test_used = "Chi-square (Yates)",
      p_value = chi$p.value
    )
  }
}


test_results <- pmap(
  list(
    summary_wide$n_tested_Female,
    summary_wide$n_sensitised_Female,
    summary_wide$n_tested_Male,
    summary_wide$n_sensitised_Male
  ),
  compare_sexes
)

summary_wide <- summary_wide %>%
  mutate(
    test_used = map_chr(
      test_results,
      "test_used"
    ),

    p_value_raw = map_dbl(
      test_results,
      "p_value"
    ),

    # BH/FDR correction across the 21 allergen comparisons displayed here.
    p_value_BH = p.adjust(
      p_value_raw,
      method = "BH"
    )
  )


# -----------------------------------------------------------------------------
# 12. Save unformatted statistical audit table
# -----------------------------------------------------------------------------

write_csv(
  summary_wide,
  file.path(
    output_dir,
    "SupplementaryTable1_sex_specific_numeric_audit.csv"
  )
)


# -----------------------------------------------------------------------------
# 13. Build manuscript-ready Supplementary Table 1
# -----------------------------------------------------------------------------

supp_table1 <- summary_wide %>%
  arrange(order) %>%
  mutate(
    `Female: Sensitised (n)` = sprintf(
      "%s (%.2f%%)",
      format(
        n_sensitised_Female,
        big.mark = ",",
        scientific = FALSE
      ),
      prevalence_pct_Female
    ),

    `Male: Sensitised (n)` = sprintf(
      "%s (%.2f%%)",
      format(
        n_sensitised_Male,
        big.mark = ",",
        scientific = FALSE
      ),
      prevalence_pct_Male
    ),

    `p-value (adjusted)` = case_when(
      is.na(p_value_BH) ~ NA_character_,

      p_value_BH < 0.001 ~ "<0.001",

      TRUE ~ sprintf(
        "%.3f",
        p_value_BH
      )
    )
  ) %>%
  transmute(
    `Allergen group` = allergen_group,
    `Allergen` = allergen_label,

    `Female: Tested (n)` = n_tested_Female,

    `Female: Sensitised (n)` =
      .data[["Female: Sensitised (n)"]],

    `Male: Tested (n)` = n_tested_Male,

    `Male: Sensitised (n)` =
      .data[["Male: Sensitised (n)"]],

    `p-value (adjusted)` =
      .data[["p-value (adjusted)"]]
  )


# -----------------------------------------------------------------------------
# 14. Print and save final table
# -----------------------------------------------------------------------------

print(
  supp_table1,
  n = Inf
)

write_csv(
  supp_table1,
  file.path(
    output_dir,
    "SupplementaryTable1_sex_specific_sensitisation.csv"
  )
)


# -----------------------------------------------------------------------------
# 15. Reproducibility information
# -----------------------------------------------------------------------------

cat(
  "\nSupplementary Table 1 completed.\n",
  "Analytical unit: sIgE test record\n",
  "Age: >=12 years\n",
  "Sex comparison: Female vs Male\n",
  "Sensitisation threshold: sIgE >= 0.35 kU/L\n",
  "Multiplicity correction: Benjamini-Hochberg across table allergens\n",
  sep = ""
)

capture.output(
  sessionInfo(),
  file = file.path(
    output_dir,
    "SupplementaryTable1_sessionInfo.txt"
  )
)
