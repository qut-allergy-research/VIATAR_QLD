# =============================================================================
# VIATAR_QLD
# Table 2. Aeroallergen sensitisation according to coastal proximity
#          of residence
#
# PURPOSE
# Reproduce the coastal-proximity analysis reported in Table 2 of the manuscript.
#
# ANALYTICAL UNIT
# IMPORTANT: this is a TEST-RECORD-LEVEL analysis.
# Each eligible sIgE test record is retained as a separate observation.
# Individuals may therefore contribute multiple test records.
#
# COASTAL GROUPS
#   <15 km
#   15–30 km
#   >30 km
#
# POSITIVE sIgE
#   sIgE >= 0.35 kU/L
#
# TABLE OUTPUT
# For each selected aeroallergen:
#   - n sensitised / N tested (%) in each coastal group
#   - relative percentage change: >30 km versus <15 km
#   - risk difference (RD, percentage points) with 95% CI
#   - odds ratio (OR) with 95% CI
#   - overall 3-group Chi-square/Fisher p-value
#   - Benjamini-Hochberg adjusted p-value
#   - pairwise z-test grouping letters with Bonferroni correction
#
# REFERENCE GROUP FOR EFFECT SIZES
#   <15 km
#
# DATA PRIVACY
# Patient-level pathology data are NOT included in the public repository.
#
# NOTE
# Earlier working code collapsed repeated records to one patient x allergen
# observation. That does NOT match the current manuscript Table 2 footnote.
# This cleaned script does not perform patient-level collapsing.
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
library(multcompView)


# -----------------------------------------------------------------------------
# 2. Local paths
# -----------------------------------------------------------------------------
#
# Keep both input files outside the public GitHub repository, or inside a local
# data/ directory excluded by .gitignore.
# -----------------------------------------------------------------------------

sige_path <- file.path(
  "data",
  "VIATAR_filtered_QLD_01122025.xlsx"
)

postcode_coast_path <- file.path(
  "data",
  "postcode_coastline.xlsx"
)

output_dir <- file.path(
  "outputs",
  "table2"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(sige_path)) {
  stop(
    "sIgE file not found: ", sige_path,
    "\nUpdate 'sige_path' locally."
  )
}

if (!file.exists(postcode_coast_path)) {
  stop(
    "Postcode coastal-proximity file not found: ",
    postcode_coast_path,
    "\nUpdate 'postcode_coast_path' locally."
  )
}


# -----------------------------------------------------------------------------
# 3. Analysis parameters
# -----------------------------------------------------------------------------

positive_threshold <- 0.35

coast_levels <- c(
  "<15 km",
  "15–30 km",
  ">30 km"
)

reference_group <- "<15 km"
comparison_group <- ">30 km"


# -----------------------------------------------------------------------------
# 4. Table 2 aeroallergens
# -----------------------------------------------------------------------------
#
# regex is matched against lower-case, squished Test.Name..Allergen. values.
# The order is deliberate: e.g. Bahia is matched before Paspalum so that
# "Bahia grass (Paspalum)" is not accidentally classified as Paspalum grass.
#
# The mapping audit written later allows the exact raw labels used for every
# publication label to be checked before manuscript submission.
# -----------------------------------------------------------------------------

table2_targets <- tibble::tribble(
  ~order, ~allergen_group,  ~allergen_label,                   ~regex,

   1L,    "Animal dander",  "Cat epithelium & dander",         "^cat epithelium.*dander$",
   2L,    "Animal dander",  "Dog dander",                      "^dog dander$",

   3L,    "Cockroach",      "Cockroach",                       "cockroach|blattella|blatella",

   4L,    "Grass pollen",   "Bahia grass",                     "^bahia grass",
   5L,    "Grass pollen",   "Bermuda grass",                   "^bermuda|^couch grass",
   6L,    "Grass pollen",   "Johnson grass",                   "^johnson grass",
   7L,    "Grass pollen",   "Paspalum grass",                  "^paspalum",
   8L,    "Grass pollen",   "Rye grass",                       "perennial rye|^rye grass",
   9L,    "Grass pollen",   "Timothy grass",                   "timothy",

  10L,    "House dust mite","Dermatophagoides farinae",        "dermatophagoides farinae|\\bd\\.\\s*farinae\\b",
  11L,    "House dust mite","Dermatophagoides pteronyssinus",  "dermatophagoides pteronyssinus|\\bd\\.\\s*pteronyssinus\\b",

  12L,    "Mould",          "Alternaria alternata",             "alternaria alternata",
  13L,    "Mould",          "Aspergillus fumigatus",            "aspergillus fumigatus",
  14L,    "Mould",          "Cladosporium herbarum",            "cladosporium herbarum",

  15L,    "Tree pollen",    "Acacia",                           "^acacia$",
  16L,    "Tree pollen",    "Birch",                            "birch|betula",
  17L,    "Tree pollen",    "Casuarina",                        "casuarina",
  18L,    "Tree pollen",    "Olive pollen",                     "olive",

  19L,    "Weed pollen",    "Parietaria judaica",               "parietaria judaica",
  20L,    "Weed pollen",    "Plantain",                         "plantain",
  21L,    "Weed pollen",    "Ragweed",                          "ragweed"
)


# -----------------------------------------------------------------------------
# 5. Helper: match a raw allergen name to one Table 2 allergen
# -----------------------------------------------------------------------------

match_table2_allergen <- function(x) {

  if (is.na(x) || x == "") {
    return(NA_integer_)
  }

  hit <- which(
    str_detect(
      x,
      table2_targets$regex
    )
  )

  if (length(hit) == 0) {
    return(NA_integer_)
  }

  # First match wins; target order is intentionally defined above.
  hit[1]
}


# -----------------------------------------------------------------------------
# 6. Read input files
# -----------------------------------------------------------------------------

sige_raw <- read_excel(
  sige_path
)

postcode_raw <- read_excel(
  postcode_coast_path
)


# -----------------------------------------------------------------------------
# 7. Prepare postcode -> coastal-proximity lookup
# -----------------------------------------------------------------------------
#
# Expected columns from the original analysis:
#   postcode
#   CoastGroup
#
# Normalisation below accepts common hyphen/en-dash variants.
# -----------------------------------------------------------------------------

required_postcode_columns <- c(
  "postcode",
  "CoastGroup"
)

missing_postcode_columns <- setdiff(
  required_postcode_columns,
  names(postcode_raw)
)

if (length(missing_postcode_columns) > 0) {
  stop(
    "Missing postcode lookup column(s): ",
    paste(missing_postcode_columns, collapse = ", ")
  )
}

postcode_lookup <- postcode_raw %>%
  transmute(
    postcode = str_trim(
      as.character(postcode)
    ),

    coast_raw = str_squish(
      as.character(CoastGroup)
    )
  ) %>%
  mutate(
    postcode = str_pad(
      postcode,
      width = 4,
      side = "left",
      pad = "0"
    ),

    CoastGroup = case_when(
      str_detect(coast_raw, "^<\\s*15") ~ "<15 km",

      str_detect(
        coast_raw,
        "15\\s*[–—-]\\s*30"
      ) ~ "15–30 km",

      str_detect(coast_raw, "^>\\s*30") ~ ">30 km",

      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    !is.na(postcode)
  ) %>%
  distinct(
    postcode,
    .keep_all = TRUE
  ) %>%
  select(
    postcode,
    CoastGroup
  )


# -----------------------------------------------------------------------------
# 8. Prepare test-record-level sIgE data
# -----------------------------------------------------------------------------
#
# IMPORTANT VALUE HANDLING
# Values such as "<0.35" must remain NEGATIVE.
# parse_number("<0.35") returns 0.35, so the '<' flag is retained explicitly.
# -----------------------------------------------------------------------------

required_sige_columns <- c(
  "Post.Code",
  "Test.Name..Allergen.",
  "test.Value..SPT.or.spIgE.values.of.each.allergen.or.mmD.for.SPT"
)

missing_sige_columns <- setdiff(
  required_sige_columns,
  names(sige_raw)
)

if (length(missing_sige_columns) > 0) {
  stop(
    "Missing sIgE dataset column(s): ",
    paste(missing_sige_columns, collapse = ", ")
  )
}

test_records <- sige_raw %>%
  transmute(
    postcode = str_pad(
      str_trim(
        as.character(`Post.Code`)
      ),
      width = 4,
      side = "left",
      pad = "0"
    ),

    allergen_raw = as.character(
      `Test.Name..Allergen.`
    ),

    sige_raw = as.character(
      `test.Value..SPT.or.spIgE.values.of.each.allergen.or.mmD.for.SPT`
    )
  ) %>%
  mutate(
    allergen_clean = allergen_raw %>%
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

    target_index = vapply(
      allergen_clean,
      match_table2_allergen,
      integer(1)
    )
  ) %>%
  filter(
    !is.na(target_index),
    !is.na(sige),
    !is.na(postcode)
  ) %>%
  mutate(
    allergen_group = table2_targets$allergen_group[target_index],
    allergen_label = table2_targets$allergen_label[target_index],
    display_order = table2_targets$order[target_index],

    sensitised = case_when(
      below_limit ~ 0L,
      sige >= positive_threshold ~ 1L,
      TRUE ~ 0L
    )
  ) %>%
  left_join(
    postcode_lookup,
    by = "postcode"
  ) %>%
  filter(
    !is.na(CoastGroup)
  ) %>%
  mutate(
    CoastGroup = factor(
      CoastGroup,
      levels = coast_levels
    )
  )


# -----------------------------------------------------------------------------
# 9. Mapping audit
# -----------------------------------------------------------------------------
#
# This should be inspected after the first run. It shows which exact raw
# laboratory allergen labels were assigned to each publication label.
# -----------------------------------------------------------------------------

mapping_audit <- test_records %>%
  count(
    allergen_group,
    allergen_label,
    allergen_raw,
    name = "n_test_records"
  ) %>%
  arrange(
    display_order = match(
      allergen_label,
      table2_targets$allergen_label
    ),
    desc(n_test_records)
  ) %>%
  select(
    -display_order
  )

write_csv(
  mapping_audit,
  file.path(
    output_dir,
    "Table2_allergen_mapping_audit.csv"
  )
)


# -----------------------------------------------------------------------------
# 10. Descriptive counts by allergen x coastal group
# -----------------------------------------------------------------------------

prevalence_long <- test_records %>%
  group_by(
    display_order,
    allergen_group,
    allergen_label,
    CoastGroup
  ) %>%
  summarise(
    tested = n(),
    sensitised = sum(
      sensitised,
      na.rm = TRUE
    ),
    prevalence = sensitised / tested,
    prevalence_pct = 100 * prevalence,
    .groups = "drop"
  ) %>%
  arrange(
    display_order,
    CoastGroup
  )


# -----------------------------------------------------------------------------
# 11. Overall 3-group test for each allergen
# -----------------------------------------------------------------------------

overall_test_one_allergen <- function(dat) {

  tab <- table(
    dat$CoastGroup,
    dat$sensitised
  )

  # Need all 3 coastal groups and both outcome levels to test.
  if (
    nrow(tab) < 2 ||
    ncol(tab) < 2
  ) {
    return(
      tibble(
        test_used = "Not testable",
        p_value = NA_real_
      )
    )
  }

  chi <- suppressWarnings(
    chisq.test(
      tab,
      correct = FALSE
    )
  )

  if (any(chi$expected < 5)) {

    p <- tryCatch(
      fisher.test(tab)$p.value,
      error = function(e) NA_real_
    )

    tibble(
      test_used = "Fisher",
      p_value = p
    )

  } else {

    tibble(
      test_used = "Chi-square",
      p_value = chi$p.value
    )
  }
}

overall_tests <- test_records %>%
  group_by(
    display_order,
    allergen_group,
    allergen_label
  ) %>%
  group_modify(
    ~ overall_test_one_allergen(.x)
  ) %>%
  ungroup() %>%
  mutate(
    p_adj_BH = p.adjust(
      p_value,
      method = "BH"
    ),

    significant_BH =
      !is.na(p_adj_BH) &
      p_adj_BH < 0.05
  )


# -----------------------------------------------------------------------------
# 12. Pairwise proportion z-tests + Bonferroni grouping letters
# -----------------------------------------------------------------------------

pairwise_letters_one_allergen <- function(dat) {

  groups_present <- coast_levels[
    coast_levels %in% as.character(
      unique(dat$CoastGroup)
    )
  ]

  # Return blanks if fewer than 2 coastal groups exist.
  if (length(groups_present) < 2) {
    return(
      tibble(
        CoastGroup = groups_present,
        group_letter = NA_character_
      )
    )
  }

  pairs <- combn(
    groups_present,
    2,
    simplify = FALSE
  )

  pairwise <- map_dfr(
    pairs,
    function(pair) {

      d1 <- dat %>%
        filter(
          CoastGroup == pair[1]
        )

      d2 <- dat %>%
        filter(
          CoastGroup == pair[2]
        )

      x <- c(
        sum(d1$sensitised == 1),
        sum(d2$sensitised == 1)
      )

      n <- c(
        nrow(d1),
        nrow(d2)
      )

      p <- if (all(n > 0)) {
        suppressWarnings(
          prop.test(
            x = x,
            n = n,
            correct = FALSE
          )$p.value
        )
      } else {
        NA_real_
      }

      tibble(
        group1 = pair[1],
        group2 = pair[2],
        p_value = p
      )
    }
  ) %>%
    mutate(
      p_adj_bonferroni = p.adjust(
        p_value,
        method = "bonferroni"
      )
    )

  # multcompLetters accepts a named vector such as:
  # "<15 km-15–30 km" = adjusted p-value
  valid_pairs <- pairwise %>%
    filter(
      !is.na(p_adj_bonferroni)
    )

  if (nrow(valid_pairs) == 0) {

    letters <- setNames(
      rep("", length(groups_present)),
      groups_present
    )

  } else {

    p_vector <- valid_pairs$p_adj_bonferroni

    names(p_vector) <- paste(
      valid_pairs$group1,
      valid_pairs$group2,
      sep = "-"
    )

    letters <- multcompView::multcompLetters(
      p_vector,
      threshold = 0.05
    )$Letters
  }

  tibble(
    CoastGroup = names(letters),
    group_letter = unname(letters)
  )
}

letters_table <- test_records %>%
  group_by(
    display_order,
    allergen_group,
    allergen_label
  ) %>%
  group_modify(
    ~ pairwise_letters_one_allergen(.x)
  ) %>%
  ungroup()

write_csv(
  letters_table,
  file.path(
    output_dir,
    "Table2_pairwise_group_letters.csv"
  )
)


# -----------------------------------------------------------------------------
# 13. Effect size helpers: >30 km versus <15 km
# -----------------------------------------------------------------------------

risk_difference_ci <- function(
  x_ref,
  n_ref,
  x_cmp,
  n_cmp
) {

  p_ref <- x_ref / n_ref
  p_cmp <- x_cmp / n_cmp

  rd <- p_cmp - p_ref

  se <- sqrt(
    p_ref * (1 - p_ref) / n_ref +
      p_cmp * (1 - p_cmp) / n_cmp
  )

  tibble(
    rd_pp = 100 * rd,
    rd_low_pp = 100 * (rd - 1.96 * se),
    rd_high_pp = 100 * (rd + 1.96 * se)
  )
}


odds_ratio_ci <- function(
  x_ref,
  n_ref,
  x_cmp,
  n_cmp
) {

  # 2x2 cells:
  # comparison positive / negative vs reference positive / negative
  a <- x_cmp
  b <- n_cmp - x_cmp
  c <- x_ref
  d <- n_ref - x_ref

  # Haldane-Anscombe continuity correction only if a zero cell occurs.
  if (any(c(a, b, c, d) == 0)) {
    a <- a + 0.5
    b <- b + 0.5
    c <- c + 0.5
    d <- d + 0.5
  }

  OR <- (a / b) / (c / d)

  se_log_or <- sqrt(
    1 / a +
      1 / b +
      1 / c +
      1 / d
  )

  tibble(
    OR = OR,
    OR_low = exp(
      log(OR) - 1.96 * se_log_or
    ),
    OR_high = exp(
      log(OR) + 1.96 * se_log_or
    )
  )
}


# -----------------------------------------------------------------------------
# 14. Compute >30 km vs <15 km effects
# -----------------------------------------------------------------------------

prevalence_wide_numeric <- prevalence_long %>%
  select(
    display_order,
    allergen_group,
    allergen_label,
    CoastGroup,
    tested,
    sensitised,
    prevalence,
    prevalence_pct
  ) %>%
  pivot_wider(
    names_from = CoastGroup,
    values_from = c(
      tested,
      sensitised,
      prevalence,
      prevalence_pct
    ),
    names_sep = "__"
  )

effect_sizes <- prevalence_wide_numeric %>%
  rowwise() %>%
  mutate(
    relative_change_pct = 100 * (
      .data[["prevalence__>30 km"]] -
        .data[["prevalence__<15 km"]]
    ) /
      .data[["prevalence__<15 km"]],

    rd = list(
      risk_difference_ci(
        x_ref = .data[["sensitised__<15 km"]],
        n_ref = .data[["tested__<15 km"]],
        x_cmp = .data[["sensitised__>30 km"]],
        n_cmp = .data[["tested__>30 km"]]
      )
    ),

    or_result = list(
      odds_ratio_ci(
        x_ref = .data[["sensitised__<15 km"]],
        n_ref = .data[["tested__<15 km"]],
        x_cmp = .data[["sensitised__>30 km"]],
        n_cmp = .data[["tested__>30 km"]]
      )
    )
  ) %>%
  unnest_wider(rd) %>%
  unnest_wider(or_result) %>%
  ungroup()

write_csv(
  effect_sizes,
  file.path(
    output_dir,
    "Table2_effect_sizes_gt30_vs_lt15.csv"
  )
)


# -----------------------------------------------------------------------------
# 15. Add pairwise letters to prevalence summaries
# -----------------------------------------------------------------------------

prevalence_with_letters <- prevalence_long %>%
  left_join(
    letters_table,
    by = c(
      "display_order",
      "allergen_group",
      "allergen_label",
      "CoastGroup"
    )
  ) %>%
  mutate(
    n_N_pct = sprintf(
      "%s/%s (%.1f%%)%s",
      format(
        sensitised,
        big.mark = ",",
        scientific = FALSE
      ),
      format(
        tested,
        big.mark = ",",
        scientific = FALSE
      ),
      prevalence_pct,
      if_else(
        is.na(group_letter) |
          group_letter == "",
        "",
        paste0(" ", group_letter)
      )
    )
  )


# -----------------------------------------------------------------------------
# 16. Create one manuscript-ready Table 2 row per allergen
# -----------------------------------------------------------------------------

table2_prevalence_wide <- prevalence_with_letters %>%
  select(
    display_order,
    allergen_group,
    allergen_label,
    CoastGroup,
    n_N_pct
  ) %>%
  pivot_wider(
    names_from = CoastGroup,
    values_from = n_N_pct
  )

table2 <- table2_prevalence_wide %>%
  left_join(
    effect_sizes %>%
      select(
        display_order,
        allergen_group,
        allergen_label,
        relative_change_pct,
        rd_pp,
        rd_low_pp,
        rd_high_pp,
        OR,
        OR_low,
        OR_high
      ),
    by = c(
      "display_order",
      "allergen_group",
      "allergen_label"
    )
  ) %>%
  left_join(
    overall_tests %>%
      select(
        display_order,
        allergen_group,
        allergen_label,
        test_used,
        p_value,
        p_adj_BH,
        significant_BH
      ),
    by = c(
      "display_order",
      "allergen_group",
      "allergen_label"
    )
  ) %>%
  arrange(
    display_order
  ) %>%
  mutate(
    `l %` = if_else(
      relative_change_pct >= 0,
      sprintf(
        "+%.1f%%",
        relative_change_pct
      ),
      sprintf(
        "%.1f%%",
        relative_change_pct
      )
    ),

    `RD (pp, 95% CI)` = sprintf(
      "%+.2f (%+.2f, %+.2f)",
      rd_pp,
      rd_low_pp,
      rd_high_pp
    ),

    `OR (95% CI)` = sprintf(
      "%.2f (%.2f, %.2f)",
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
    ),

    `Sig.` = if_else(
      significant_BH,
      "*",
      ""
    )
  ) %>%
  transmute(
    `Allergen Group` = allergen_group,
    `Allergen` = allergen_label,
    `<15 km n/N (%)` = .data[["<15 km"]],
    `15–30 km n/N (%)` = .data[["15–30 km"]],
    `>30 km n/N (%)` = .data[[">30 km"]],
    `l %` = .data[["l %"]],
    `RD (pp, 95% CI)` = .data[["RD (pp, 95% CI)"]],
    `OR (95% CI)` = .data[["OR (95% CI)"]],
    `p-value` = .data[["p-value"]],
    `BH/FDR` = .data[["BH/FDR"]],
    `Sig.` = .data[["Sig."]]
  )


# -----------------------------------------------------------------------------
# 17. Print and save final table
# -----------------------------------------------------------------------------

print(
  table2,
  n = Inf
)

write_csv(
  table2,
  file.path(
    output_dir,
    "Table2_coastal_proximity_final.csv"
  )
)

write_csv(
  prevalence_long,
  file.path(
    output_dir,
    "Table2_prevalence_numeric.csv"
  )
)

write_csv(
  overall_tests,
  file.path(
    output_dir,
    "Table2_overall_tests.csv"
  )
)


# -----------------------------------------------------------------------------
# 18. Reproducibility checks
# -----------------------------------------------------------------------------

cat(
  "\nTEST-RECORD-LEVEL analysis completed.\n",
  "Positive threshold: sIgE >= ",
  positive_threshold,
  " kU/L\n",
  "Reference group for RD and OR: ",
  reference_group,
  "\nComparison group: ",
  comparison_group,
  "\n",
  sep = ""
)

capture.output(
  sessionInfo(),
  file = file.path(
    output_dir,
    "Table2_sessionInfo.txt"
  )
)
