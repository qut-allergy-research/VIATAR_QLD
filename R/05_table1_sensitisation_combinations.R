# =============================================================================
# VIATAR_QLD
# Table 1. Frequency of sensitisation combinations among participants with
#          >=1 positive aeroallergen group
#
# PURPOSE
# Create participant-level sensitisation combinations across the seven
# aeroallergen groups used in the manuscript.
#
# A participant is positive for an aeroallergen group if they have at least
# one positive sIgE test (>= 0.35 kU/L) to any allergen within that group.
#
# TABLE DENOMINATOR
# Percentages in Table 1 are calculated among participants with at least one
# positive aeroallergen group, NOT among all participants and NOT only among
# participants sensitised to >=2 groups.
#
# The script also produces a supporting distribution of participants sensitised
# to exactly 0, 1, 2, or >=3 aeroallergen groups, using all eligible
# participants as the denominator.
#
# DATA PRIVACY
# The patient-level pathology dataset is NOT included in this public repository.
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Packages
# -----------------------------------------------------------------------------

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)


# -----------------------------------------------------------------------------
# 2. Paths
# -----------------------------------------------------------------------------

data_path <- file.path(
  "data",
  "VIATAR_filtered_QLD_01122025.xlsx"
)

output_dir <- file.path("outputs", "table1")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(data_path)) {
  stop(
    "Data file not found: ", data_path,
    "\nUpdate 'data_path' locally. Do not upload the patient-level dataset to GitHub."
  )
}


# -----------------------------------------------------------------------------
# 3. Read data
# -----------------------------------------------------------------------------

df_raw <- read_excel(data_path)


# -----------------------------------------------------------------------------
# 4. Select participant identifier
# -----------------------------------------------------------------------------
#
# The original combination-analysis working code preferred the SLK identifier
# when available and otherwise used the VIATAR patient ID. The selected field is
# printed so the analysis is auditable.
#
# If the manuscript Table 1 was generated with a different identifier, set
# patient_id_col explicitly and rerun the script.

if ("Patients.Statistical.Linkage.Key..SLK...SLK581" %in% names(df_raw)) {

  patient_id_col <- "Patients.Statistical.Linkage.Key..SLK...SLK581"

} else if ("Viatar..Patient..ID" %in% names(df_raw)) {

  patient_id_col <- "Viatar..Patient..ID"

} else {

  stop(
    "No expected participant identifier was found in the dataset."
  )
}

message("Participant identifier used: ", patient_id_col)


# -----------------------------------------------------------------------------
# 5. Exact aeroallergen-group definitions
# -----------------------------------------------------------------------------
#
# Only the seven aeroallergen groups used in the manuscript are included.
# Penicillin/drug allergens and other miscellaneous tests are intentionally
# excluded from this analysis.

allergen_lookup <- tibble::tribble(
  ~allergen_key,                         ~allergen_group,

  # House dust mites
  "dermatophagoides farinae",            "Dust Mites",
  "dermatophagoides pteronyssinus",      "Dust Mites",
  "blomia tropicalis",                   "Dust Mites",

  # Moulds
  "alternaria alternata",                "Moulds",
  "aspergillus fumigatus",               "Moulds",
  "cladosporium herbarum",               "Moulds",
  "mould mix 1",                         "Moulds",

  # Grass pollens
  "bahia grass (paspalum)",              "Grass Pollens",
  "bermuda (couch) grass",               "Grass Pollens",
  "johnson grass (sorghum)",             "Grass Pollens",
  "perennial rye grass",                 "Grass Pollens",
  "grass mix 2",                         "Grass Pollens",

  # Weed pollens
  "common ragweed",                      "Weed Pollens",
  "weed mix 1",                          "Weed Pollens",

  # Tree pollens
  "eucalyptus",                          "Tree Pollens",
  "gum tree",                            "Tree Pollens",
  "eucalyptus, gum tree",                "Tree Pollens",
  "acacia",                              "Tree Pollens",
  "tree mix 7",                          "Tree Pollens",

  # Animal dander / feathers
  "cat epithelium & dander",             "Animal Dander",
  "cat epithelium and dander",           "Animal Dander",
  "dog dander",                          "Animal Dander",
  "chicken feathers",                    "Animal Dander",
  "animal mix 1",                        "Animal Dander",
  "feather mix",                         "Animal Dander",

  # Cockroach
  "cockroach (blatella germanica)",      "Cockroach",
  "cockroach (blattella germanica)",     "Cockroach",
  "blattella germanica",                 "Cockroach"
)

group_order <- c(
  "Animal Dander",
  "Cockroach",
  "Dust Mites",
  "Grass Pollens",
  "Moulds",
  "Tree Pollens",
  "Weed Pollens"
)


# -----------------------------------------------------------------------------
# 6. Prepare participant, allergen, and sIgE variables
# -----------------------------------------------------------------------------
#
# Values reported as "<0.35" are explicitly treated as negative.
# Values such as ">100" remain positive after numeric parsing.

df <- df_raw %>%
  transmute(
    participant_id = as.character(.data[[patient_id_col]]),
    age = suppressWarnings(as.numeric(`Age.in.Years`)),
    allergen_raw = as.character(`Test.Name..Allergen.`),
    sige_raw = as.character(
      `test.Value..SPT.or.spIgE.values.of.each.allergen.or.mmD.for.SPT`
    )
  ) %>%
  mutate(
    allergen_key = allergen_raw %>%
      str_to_lower() %>%
      str_squish(),

    sige_text = str_trim(sige_raw),
    below_limit = str_detect(sige_text, "^<"),
    sige = parse_number(sige_text)
  ) %>%
  filter(
    !is.na(participant_id),
    participant_id != "",
    !is.na(age),
    age >= 12
  ) %>%
  left_join(
    allergen_lookup,
    by = "allergen_key"
  ) %>%
  mutate(
    positive = case_when(
      is.na(sige) ~ NA_integer_,
      below_limit ~ 0L,
      sige >= 0.35 ~ 1L,
      TRUE ~ 0L
    )
  )


# -----------------------------------------------------------------------------
# 7. Audit unmatched allergen labels
# -----------------------------------------------------------------------------

unmatched_allergens <- df %>%
  filter(
    is.na(allergen_group),
    !is.na(allergen_raw),
    allergen_raw != ""
  ) %>%
  count(
    allergen_raw,
    sort = TRUE,
    name = "n_records"
  )

write_csv(
  unmatched_allergens,
  file.path(
    output_dir,
    "Table1_unmatched_allergen_labels.csv"
  )
)


# -----------------------------------------------------------------------------
# 8. Create one participant x aeroallergen-group positivity record
# -----------------------------------------------------------------------------
#
# If a participant has repeated tests within the same group, the group is
# considered positive if ANY mapped test in that group is positive.

participant_group <- df %>%
  filter(!is.na(allergen_group)) %>%
  group_by(
    participant_id,
    allergen_group
  ) %>%
  summarise(
    group_positive = as.integer(
      any(positive == 1L, na.rm = TRUE)
    ),
    .groups = "drop"
  )


# -----------------------------------------------------------------------------
# 9. Retain positive participant-group combinations only
# -----------------------------------------------------------------------------

positive_participant_groups <- participant_group %>%
  filter(group_positive == 1L) %>%
  distinct(
    participant_id,
    allergen_group
  )


# -----------------------------------------------------------------------------
# 10. Build one sensitisation combination per positive participant
# -----------------------------------------------------------------------------

participant_combinations <- positive_participant_groups %>%
  mutate(
    allergen_group = factor(
      allergen_group,
      levels = group_order
    )
  ) %>%
  arrange(
    participant_id,
    allergen_group
  ) %>%
  group_by(participant_id) %>%
  summarise(
    n_positive_groups = n_distinct(allergen_group),

    combination = {
      groups <- sort(as.character(allergen_group))

      if (length(groups) == 1) {
        paste0(groups, " only")
      } else {
        paste(groups, collapse = " + ")
      }
    },

    .groups = "drop"
  )


# -----------------------------------------------------------------------------
# 11. Table 1 denominator
# -----------------------------------------------------------------------------

n_participants_any_positive <- n_distinct(
  participant_combinations$participant_id
)

message(
  "Participants with >=1 positive aeroallergen group: ",
  format(n_participants_any_positive, big.mark = ",")
)


# -----------------------------------------------------------------------------
# 12. Full combination-frequency table
# -----------------------------------------------------------------------------
#
# IMPORTANT:
# Percentages use ALL participants with >=1 positive aeroallergen group as
# the denominator. They are not recalculated within the displayed top 10 rows.

combination_table_full <- participant_combinations %>%
  count(
    combination,
    name = "n"
  ) %>%
  mutate(
    percent = 100 * n / n_participants_any_positive
  ) %>%
  arrange(
    desc(n),
    combination
  )


# -----------------------------------------------------------------------------
# 13. Manuscript Table 1: 10 most frequent combinations
# -----------------------------------------------------------------------------

table1 <- combination_table_full %>%
  slice_head(n = 10) %>%
  transmute(
    `Combination of Allergen Groups` = combination,
    n = n,
    `%` = round(percent, 1)
  )

print(table1)


# -----------------------------------------------------------------------------
# 14. Save Table 1 and the complete combination table
# -----------------------------------------------------------------------------

write_csv(
  table1,
  file.path(
    output_dir,
    "Table1_top10_sensitisation_combinations.csv"
  )
)

write_csv(
  combination_table_full,
  file.path(
    output_dir,
    "Table1_all_sensitisation_combinations.csv"
  )
)


# -----------------------------------------------------------------------------
# 15. Supporting analysis: number of positive groups per participant
# -----------------------------------------------------------------------------
#
# This supports the Results statement describing the proportion of all included
# participants sensitised to exactly one, exactly two, or >=3 allergen groups.

all_eligible_participants <- df %>%
  distinct(participant_id)

participant_positive_group_count <- all_eligible_participants %>%
  left_join(
    participant_combinations %>%
      select(
        participant_id,
        n_positive_groups
      ),
    by = "participant_id"
  ) %>%
  mutate(
    n_positive_groups = replace_na(
      n_positive_groups,
      0L
    ),

    category = case_when(
      n_positive_groups == 0 ~ "0 groups",
      n_positive_groups == 1 ~ "1 group",
      n_positive_groups == 2 ~ "2 groups",
      n_positive_groups >= 3 ~ "3+ groups"
    )
  )

n_all_participants <- n_distinct(
  participant_positive_group_count$participant_id
)

group_count_distribution <- participant_positive_group_count %>%
  count(
    category,
    name = "n"
  ) %>%
  complete(
    category = c(
      "0 groups",
      "1 group",
      "2 groups",
      "3+ groups"
    ),
    fill = list(n = 0)
  ) %>%
  mutate(
    category = factor(
      category,
      levels = c(
        "0 groups",
        "1 group",
        "2 groups",
        "3+ groups"
      )
    ),
    percent_all_participants =
      100 * n / n_all_participants
  ) %>%
  arrange(category)

print(group_count_distribution)

write_csv(
  group_count_distribution,
  file.path(
    output_dir,
    "Table1_supporting_number_of_positive_groups.csv"
  )
)


# -----------------------------------------------------------------------------
# 16. Reproducibility summary
# -----------------------------------------------------------------------------

cat(
  "\nParticipants included in denominator for Table 1:",
  format(n_participants_any_positive, big.mark = ","),
  "\nTotal eligible participants:",
  format(n_all_participants, big.mark = ","),
  "\n"
)

capture.output(
  sessionInfo(),
  file = file.path(
    output_dir,
    "Table1_sessionInfo.txt"
  )
)
