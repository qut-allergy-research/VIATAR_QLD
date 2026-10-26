# =============================================================================
# VIATAR_QLD
# Supplementary Table 3. Frequency of aeroallergen sensitisation in urban
#                        and regional areas
#
# PURPOSE
# Reproduce the urban-versus-regional sensitisation analysis reported in
# Supplementary Table 3.
#
# ANALYTICAL UNIT
# TEST RECORD.
# Each eligible sIgE test record is retained as one observation; individuals
# may contribute multiple records.
#
# ELIGIBILITY
#   - Queensland residential postcode: 4000–4999
#   - age >= 12 years
#   - postcode successfully classified as Urban or Regional
#   - selected aeroallergens shown in Supplementary Table 3
#   - usable sIgE result
#
# SENSITISATION DEFINITION
#   sIgE >= 0.35 kU/L
#
# EFFECT DIRECTION
#   Urban minus Regional for risk difference (RD)
#   Urban / Regional for odds ratio (OR)
#
# STATISTICAL METHODS
#   - 2 x 2 Pearson chi-square test without Yates correction
#   - Fisher's exact test only if an expected cell count is < 5
#   - RD 95% CI: Newcombe method using Wilson score intervals
#   - OR 95% CI: Wald interval on log(OR), with 0.5 correction if a zero cell
#   - Benjamini-Hochberg FDR correction across the 22 allergens in this table
#
# SIGNIFICANCE LABELS
#   BH-adjusted P < 0.001 : ***
#   BH-adjusted P < 0.01  : **
#   BH-adjusted P < 0.05  : *
#   otherwise             : ns
#
# DATA PRIVACY
# The real pathology dataset and postcode lookup are NOT included in the
# public GitHub repository.
#
# IMPORTANT
# Earlier working scripts included patient-level collapsing and Bonferroni
# correction. Those versions do NOT match the current Supplementary Table 3
# footnote. This cleaned script follows the final table specification.
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

postcode_lookup_path <- file.path(
  "data",
  "postcode_urbanRegion.xlsx"
)

output_dir <- file.path(
  "outputs",
  "supplementary_table3"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(data_path)) {
  stop(
    "Main sIgE data file not found: ", data_path,
    "\nUpdate 'data_path' locally. Do not upload the real dataset to GitHub."
  )
}

if (!file.exists(postcode_lookup_path)) {
  stop(
    "Urban/regional postcode lookup not found: ", postcode_lookup_path,
    "\nUpdate 'postcode_lookup_path' locally."
  )
}


# -----------------------------------------------------------------------------
# 3. Analysis parameters
# -----------------------------------------------------------------------------

positive_threshold <- 0.35
qld_postcode_min <- 4000L
qld_postcode_max <- 4999L


# -----------------------------------------------------------------------------
# 4. Supplementary Table 3 target allergens
# -----------------------------------------------------------------------------

allergen_lookup <- tibble::tribble(
  ~order, ~allergen_group,      ~allergen_label,              ~raw_key,

   1L,    "House dust mites",   "D. pteronyssinus",           "dermatophagoides pteronyssinus",
   2L,    "House dust mites",   "D. farinae",                 "dermatophagoides farinae",
   3L,    "House dust mites",   "Blomia tropicalis",          "blomia tropicalis",

   4L,    "Grass pollen",       "Grass Mix 2",                "grass mix 2",
   5L,    "Grass pollen",       "Perennial ryegrass",         "perennial rye grass",
   5L,    "Grass pollen",       "Perennial ryegrass",         "perennial ryegrass",
   6L,    "Grass pollen",       "Bahia grass",                "bahia grass (paspalum)",
   6L,    "Grass pollen",       "Bahia grass",                "bahia grass",
   7L,    "Grass pollen",       "Bermuda grass",              "bermuda (couch) grass",
   7L,    "Grass pollen",       "Bermuda grass",              "bermuda grass",
   8L,    "Grass pollen",       "Johnson grass",              "johnson grass (sorghum)",
   8L,    "Grass pollen",       "Johnson grass",              "johnson grass",

   9L,    "Tree pollen",        "Acacia",                     "acacia",
  10L,    "Tree pollen",        "Eucalyptus / Gum tree",      "eucalyptus, gum tree",
  10L,    "Tree pollen",        "Eucalyptus / Gum tree",      "eucalyptus (gum tree)",
  10L,    "Tree pollen",        "Eucalyptus / Gum tree",      "eucalyptus",
  10L,    "Tree pollen",        "Eucalyptus / Gum tree",      "gum tree",
  11L,    "Tree pollen",        "Tree Mix 7",                 "tree mix 7",

  12L,    "Weed pollen",        "Common ragweed",             "common ragweed",
  13L,    "Weed pollen",        "Weed Mix 1",                 "weed mix 1",

  14L,    "Animal dander",      "Cat dander",                 "cat epithelium & dander",
  14L,    "Animal dander",      "Cat dander",                 "cat epithelium and dander",
  14L,    "Animal dander",      "Cat dander",                 "cat dander",
  15L,    "Animal dander",      "Dog dander",                 "dog dander",
  16L,    "Animal dander",      "Animal Mix 1",               "animal mix 1",
  17L,    "Animal dander",      "Feather Mix",                "feather mix",

  18L,    "Moulds",             "Alternaria alternata",       "alternaria alternata",
  19L,    "Moulds",             "Aspergillus fumigatus",      "aspergillus fumigatus",
  20L,    "Moulds",             "Cladosporium herbarum",      "cladosporium herbarum",
  21L,    "Moulds",             "Mould Mix 1",                "mould mix 1",

  22L,    "Cockroach",          "Cockroach",                  "cockroach (blatella germanica)",
  22L,    "Cockroach",          "Cockroach",                  "cockroach (blattella germanica)",
  22L,    "Cockroach",          "Cockroach",                  "blattella germanica"
)


# -----------------------------------------------------------------------------
# 5. Helper functions
# -----------------------------------------------------------------------------

clean_postcode <- function(x) {
  x <- as.character(x)
  p <- str_extract(
    x,
    "\\d{4}"
  )
  suppressWarnings(
    as.integer(p)
  )
}


normalise_lookup_names <- function(x) {
  x %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9]+", "_") %>%
    str_replace_all("^_|_$", "")
}


wilson_ci <- function(x, n, z = 1.96) {

  if (
    is.na(x) ||
    is.na(n) ||
    n <= 0
  ) {
    return(
      c(
        lower = NA_real_,
        upper = NA_real_
      )
    )
  }

  p <- x / n

  denominator <- 1 + z^2 / n

  centre <- (
    p + z^2 / (2 * n)
  ) / denominator

  half_width <- (
    z *
      sqrt(
        (
          p * (1 - p) +
            z^2 / (4 * n)
        ) /
          n
      )
  ) / denominator

  c(
    lower = centre - half_width,
    upper = centre + half_width
  )
}


newcombe_rd_ci <- function(
  urban_pos,
  urban_n,
  regional_pos,
  regional_n
) {

  p_urban <- urban_pos / urban_n
  p_regional <- regional_pos / regional_n

  urban_ci <- wilson_ci(
    urban_pos,
    urban_n
  )

  regional_ci <- wilson_ci(
    regional_pos,
    regional_n
  )

  c(
    rd = p_urban - p_regional,
    lower = urban_ci["lower"] -
      regional_ci["upper"],
    upper = urban_ci["upper"] -
      regional_ci["lower"]
  )
}


odds_ratio_ci <- function(
  urban_pos,
  urban_neg,
  regional_pos,
  regional_neg
) {

  cells <- c(
    urban_pos,
    urban_neg,
    regional_pos,
    regional_neg
  )

  if (any(cells == 0)) {
    urban_pos <- urban_pos + 0.5
    urban_neg <- urban_neg + 0.5
    regional_pos <- regional_pos + 0.5
    regional_neg <- regional_neg + 0.5
  }

  OR <- (
    urban_pos / urban_neg
  ) / (
    regional_pos / regional_neg
  )

  se_log_or <- sqrt(
    1 / urban_pos +
      1 / urban_neg +
      1 / regional_pos +
      1 / regional_neg
  )

  c(
    OR = OR,
    lower = exp(
      log(OR) - 1.96 * se_log_or
    ),
    upper = exp(
      log(OR) + 1.96 * se_log_or
    )
  )
}


test_2x2 <- function(
  urban_pos,
  urban_neg,
  regional_pos,
  regional_neg
) {

  tab <- matrix(
    c(
      urban_pos,
      urban_neg,
      regional_pos,
      regional_neg
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      Area = c(
        "Urban",
        "Regional"
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
      correct = FALSE
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
      test_used = "Chi-square",
      p_value = chi$p.value
    )
  }
}


# -----------------------------------------------------------------------------
# 6. Read data
# -----------------------------------------------------------------------------

raw_data <- read_excel(
  data_path,
  guess_max = 20000
)

postcode_lookup_raw <- read_excel(
  postcode_lookup_path
)


# -----------------------------------------------------------------------------
# 7. Prepare postcode -> Urban/Regional lookup
# -----------------------------------------------------------------------------

lookup_names_clean <- normalise_lookup_names(
  names(postcode_lookup_raw)
)

names(postcode_lookup_raw) <- lookup_names_clean

postcode_candidates <- names(
  postcode_lookup_raw
)[
  str_detect(
    names(postcode_lookup_raw),
    "postcode|post_code|poa"
  )
]

urban_candidates <- names(
  postcode_lookup_raw
)[
  str_detect(
    names(postcode_lookup_raw),
    "urban.*reg|urban_reg|urbanreg|urban.*regional"
  )
]

if (length(postcode_candidates) == 0) {
  stop(
    "Could not identify the postcode column in postcode_urbanRegion.xlsx."
  )
}

if (length(urban_candidates) == 0) {
  stop(
    "Could not identify the Urban/Regional classification column in postcode_urbanRegion.xlsx."
  )
}

postcode_lookup <- postcode_lookup_raw %>%
  transmute(
    postcode = clean_postcode(
      .data[[postcode_candidates[1]]]
    ),

    urban_regional_raw = as.character(
      .data[[urban_candidates[1]]]
    )
  ) %>%
  mutate(
    urban_regional = case_when(
      str_to_lower(
        str_trim(
          urban_regional_raw
        )
      ) == "urban" ~ "Urban",

      str_to_lower(
        str_trim(
          urban_regional_raw
        )
      ) == "regional" ~ "Regional",

      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    !is.na(postcode),
    !is.na(urban_regional)
  ) %>%
  distinct(
    postcode,
    .keep_all = TRUE
  ) %>%
  select(
    postcode,
    urban_regional
  )


# -----------------------------------------------------------------------------
# 8. Check main-data columns
# -----------------------------------------------------------------------------

required_columns <- c(
  "Post.Code",
  "Age.in.Years",
  "Test.Name..Allergen.",
  "test.Value..SPT.or.spIgE.values.of.each.allergen.or.mmD.for.SPT"
)

missing_columns <- setdiff(
  required_columns,
  names(raw_data)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing required main-data column(s): ",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}


# -----------------------------------------------------------------------------
# 9. Prepare TEST-RECORD-LEVEL analysis data
# -----------------------------------------------------------------------------
#
# IMPORTANT
# No participant-level collapsing is performed.
#
# Values such as "<0.35" are retained as valid NEGATIVE tests.
# -----------------------------------------------------------------------------

analysis_data <- raw_data %>%
  transmute(
    postcode = clean_postcode(
      `Post.Code`
    ),

    age = suppressWarnings(
      as.numeric(
        `Age.in.Years`
      )
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
    )
  ) %>%
  filter(
    !is.na(postcode),
    postcode >= qld_postcode_min,
    postcode <= qld_postcode_max,
    !is.na(age),
    age >= 12,
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
  left_join(
    postcode_lookup,
    by = "postcode"
  ) %>%
  filter(
    urban_regional %in% c(
      "Urban",
      "Regional"
    )
  ) %>%
  mutate(
    sensitised = case_when(
      below_limit ~ 0L,
      sige >= positive_threshold ~ 1L,
      TRUE ~ 0L
    ),

    urban_regional = factor(
      urban_regional,
      levels = c(
        "Urban",
        "Regional"
      )
    )
  )


# -----------------------------------------------------------------------------
# 10. Audit raw allergen mapping
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
    "SupplementaryTable3_allergen_mapping_audit.csv"
  )
)


# -----------------------------------------------------------------------------
# 11. Confirm all 22 publication allergens were found
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

  print(
    missing_targets
  )

  stop(
    "One or more Supplementary Table 3 allergens were not found. ",
    "Inspect the allergen mapping audit before continuing."
  )
}


# -----------------------------------------------------------------------------
# 12. Counts and percentages by allergen x Urban/Regional
# -----------------------------------------------------------------------------

counts_long <- analysis_data %>%
  group_by(
    order,
    allergen_group,
    allergen_label,
    urban_regional
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
    urban_regional
  )


counts_wide <- counts_long %>%
  select(
    order,
    allergen_group,
    allergen_label,
    urban_regional,
    n_tested,
    n_sensitised,
    prevalence_pct
  ) %>%
  pivot_wider(
    names_from = urban_regional,
    values_from = c(
      n_tested,
      n_sensitised,
      prevalence_pct
    ),
    names_sep = "_"
  ) %>%
  arrange(order)


# -----------------------------------------------------------------------------
# 13. Effect sizes and p-values
# -----------------------------------------------------------------------------

calculate_one_row <- function(
  urban_tested,
  urban_pos,
  regional_tested,
  regional_pos
) {

  urban_neg <- urban_tested - urban_pos
  regional_neg <- regional_tested - regional_pos

  rd <- newcombe_rd_ci(
    urban_pos = urban_pos,
    urban_n = urban_tested,
    regional_pos = regional_pos,
    regional_n = regional_tested
  )

  or_result <- odds_ratio_ci(
    urban_pos = urban_pos,
    urban_neg = urban_neg,
    regional_pos = regional_pos,
    regional_neg = regional_neg
  )

  test_result <- test_2x2(
    urban_pos = urban_pos,
    urban_neg = urban_neg,
    regional_pos = regional_pos,
    regional_neg = regional_neg
  )

  tibble(
    rd_pp = 100 * rd["rd"],
    rd_low_pp = 100 * rd["lower"],
    rd_high_pp = 100 * rd["upper"],

    OR = or_result["OR"],
    OR_low = or_result["lower"],
    OR_high = or_result["upper"],

    test_used = test_result$test_used,
    p_value = test_result$p_value
  )
}


statistics <- pmap_dfr(
  list(
    counts_wide$n_tested_Urban,
    counts_wide$n_sensitised_Urban,
    counts_wide$n_tested_Regional,
    counts_wide$n_sensitised_Regional
  ),
  calculate_one_row
)


results_numeric <- bind_cols(
  counts_wide,
  statistics
) %>%
  mutate(
    p_adj_BH = p.adjust(
      p_value,
      method = "BH"
    ),

    significance = case_when(
      is.na(p_adj_BH) ~ NA_character_,
      p_adj_BH < 0.001 ~ "***",
      p_adj_BH < 0.01 ~ "**",
      p_adj_BH < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )


# -----------------------------------------------------------------------------
# 14. Save full numerical audit output
# -----------------------------------------------------------------------------

write_csv(
  results_numeric,
  file.path(
    output_dir,
    "SupplementaryTable3_urban_regional_numeric_audit.csv"
  )
)


# -----------------------------------------------------------------------------
# 15. Build manuscript-ready Supplementary Table 3
# -----------------------------------------------------------------------------

supp_table3 <- results_numeric %>%
  arrange(order) %>%
  mutate(
    `RD (pp, 95% CI)` = sprintf(
      "%+.2f (%+.2f to %+.2f)",
      rd_pp,
      rd_low_pp,
      rd_high_pp
    ),

    `OR (95% CI)` = sprintf(
      "%.2f (%.2f–%.2f)",
      OR,
      OR_low,
      OR_high
    ),

    `p-value` = case_when(
      is.na(p_value) ~ NA_character_,
      p_value < 0.001 ~ "<0.001",
      TRUE ~ sprintf(
        "%.3f",
        p_value
      )
    ),

    `BH/FDR` = case_when(
      is.na(p_adj_BH) ~ NA_character_,
      p_adj_BH < 0.001 ~ "<0.001",
      TRUE ~ sprintf(
        "%.3f",
        p_adj_BH
      )
    )
  ) %>%
  transmute(
    `Allergen Group` = allergen_group,
    `Allergen` = allergen_label,

    `Urban tested` = n_tested_Urban,
    `Urban sensitised` = n_sensitised_Urban,
    `Urban %` = round(
      prevalence_pct_Urban,
      2
    ),

    `Regional tested` = n_tested_Regional,
    `Regional sensitised` = n_sensitised_Regional,
    `Regional %` = round(
      prevalence_pct_Regional,
      2
    ),

    `RD (pp, 95% CI)` =
      .data[["RD (pp, 95% CI)"]],

    `OR (95% CI)` =
      .data[["OR (95% CI)"]],

    `p-value` =
      .data[["p-value"]],

    `BH/FDR` =
      .data[["BH/FDR"]],

    `Sig` = significance
  )


# -----------------------------------------------------------------------------
# 16. Print and save final table
# -----------------------------------------------------------------------------

print(
  supp_table3,
  n = Inf
)

write_csv(
  supp_table3,
  file.path(
    output_dir,
    "SupplementaryTable3_urban_regional_sensitisation.csv"
  )
)


# -----------------------------------------------------------------------------
# 17. Validation against the currently reported manuscript counts
# -----------------------------------------------------------------------------
#
# These values are used ONLY to verify reproducibility.
# They do not drive any calculation.
# -----------------------------------------------------------------------------

reported_counts <- tibble::tribble(
  ~allergen_label,             ~urban_tested_reported, ~urban_pos_reported, ~regional_tested_reported, ~regional_pos_reported,

  "D. pteronyssinus",                        70234L, 33879L, 55792L, 24182L,
  "D. farinae",                               1484L,   732L,  2198L,   898L,
  "Blomia tropicalis",                        2041L,   676L,  2975L,   974L,

  "Grass Mix 2",                             61253L, 13566L, 48972L,  9579L,
  "Perennial ryegrass",                       5427L,  1221L,  4596L,   916L,
  "Bahia grass",                              4864L,  1122L,  4525L,   990L,
  "Bermuda grass",                            4959L,  1103L,  4583L,  1024L,
  "Johnson grass",                            4720L,   930L,  4497L,   892L,

  "Acacia",                                   1991L,    87L,  1500L,    98L,
  "Eucalyptus / Gum tree",                    1967L,    81L,  1418L,    83L,
  "Tree Mix 7",                              12065L,  1210L, 10745L,  1000L,

  "Common ragweed",                           1923L,   199L,  1378L,   167L,
  "Weed Mix 1",                               7886L,   902L,  7434L,   851L,

  "Cat dander",                              10243L,  2544L,  8672L,  1626L,
  "Dog dander",                               9900L,  2228L,  9387L,  1653L,
  "Animal Mix 1",                            47234L,  8416L, 39192L,  6167L,
  "Feather Mix",                               580L,     7L,   734L,     5L,

  "Alternaria alternata",                     2745L,   342L,  1698L,   209L,
  "Aspergillus fumigatus",                   17140L,  2316L, 11183L,  1284L,
  "Cladosporium herbarum",                    1928L,    75L,  1377L,    58L,
  "Mould Mix 1",                              9786L,  1103L,  9808L,   981L,

  "Cockroach",                                1121L,   181L,  1480L,   242L
)


validation <- results_numeric %>%
  select(
    allergen_label,
    n_tested_Urban,
    n_sensitised_Urban,
    n_tested_Regional,
    n_sensitised_Regional,
    rd_pp,
    rd_low_pp,
    rd_high_pp,
    OR,
    OR_low,
    OR_high,
    p_value,
    p_adj_BH
  ) %>%
  left_join(
    reported_counts,
    by = "allergen_label"
  ) %>%
  mutate(
    urban_tested_match =
      n_tested_Urban ==
      urban_tested_reported,

    urban_positive_match =
      n_sensitised_Urban ==
      urban_pos_reported,

    regional_tested_match =
      n_tested_Regional ==
      regional_tested_reported,

    regional_positive_match =
      n_sensitised_Regional ==
      regional_pos_reported,

    all_reported_counts_match =
      urban_tested_match &
      urban_positive_match &
      regional_tested_match &
      regional_positive_match
  )


write_csv(
  validation,
  file.path(
    output_dir,
    "SupplementaryTable3_validation_against_reported_counts.csv"
  )
)

if (
  any(
    !validation$all_reported_counts_match,
    na.rm = TRUE
  )
) {
  warning(
    "At least one calculated Urban/Regional count differs from the ",
    "currently reported Supplementary Table 3. Inspect the validation CSV."
  )
}


# -----------------------------------------------------------------------------
# 18. Session information
# -----------------------------------------------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    output_dir,
    "SupplementaryTable3_sessionInfo.txt"
  )
)
