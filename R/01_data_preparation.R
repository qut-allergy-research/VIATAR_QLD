# =============================================================================
# VIATAR_QLD
# 01. Data preparation and cohort description
#
# PURPOSE
# Prepare the Queensland VIATAR analysis dataset and generate the descriptive
# cohort summaries reported in the manuscript.
#
# DATA PREPARATION RULES
#   1. Retain individuals aged >=12 years.
#   2. Exclude records with missing/unusable residential postcode.
#   3. Retain Queensland postcodes only (4000–4999).
#   4. Exclude records with missing suburb/locality.
#   5. Exclude records with missing allergen name.
#   6. Exclude records with missing/unusable sIgE result.
#
# IMPORTANT
# Values such as "<0.35" are valid laboratory results and must NOT be treated
# as missing. They are retained in the prepared dataset.
#
# COHORT DESCRIPTION PRODUCED BY THIS SCRIPT
#   - total number of test records
#   - number of unique individuals
#   - median age, Q1, Q3 and IQR
#   - female/male counts and percentages
#   - age-group distribution
#
# Both test-record-level and participant-level demographic summaries are saved
# because the study contains repeated tests from some individuals.
#
# DATA PRIVACY
# Raw and filtered pathology datasets must NOT be committed to the public
# GitHub repository. Keep data/ excluded using .gitignore.
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Packages
# -----------------------------------------------------------------------------

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(openxlsx)


# -----------------------------------------------------------------------------
# 2. Local file paths
# -----------------------------------------------------------------------------
#
# Change only these filenames locally if your raw data file has another name.
# -----------------------------------------------------------------------------

input_path <- file.path(
  "data",
  "VIATAR_raw.xlsx"
)

output_path <- file.path(
  "data",
  "VIATAR_filtered_QLD_01122025.xlsx"
)

output_dir <- file.path(
  "outputs",
  "data_preparation"
)

dir.create(
  dirname(output_path),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(input_path)) {
  stop(
    "Raw data file not found: ", input_path,
    "\nUpdate 'input_path' locally. Do not upload the real dataset to GitHub."
  )
}


# -----------------------------------------------------------------------------
# 3. Read raw data
# -----------------------------------------------------------------------------

raw_data <- read_excel(
  input_path,
  sheet = 1,
  guess_max = 20000
) %>%
  as_tibble()

message(
  "Raw data loaded: ",
  format(nrow(raw_data), big.mark = ","),
  " test records"
)


# -----------------------------------------------------------------------------
# 4. Column names
# -----------------------------------------------------------------------------

patient_id_col <- "Viatar..Patient..ID"
postcode_col <- "Post.Code"
age_col <- "Age.in.Years"
sex_col <- "Sex.at.birth"
allergen_col <- "Test.Name..Allergen."
sige_col <- "test.Value..SPT.or.spIgE.values.of.each.allergen.or.mmD.for.SPT"

required_columns <- c(
  patient_id_col,
  postcode_col,
  age_col,
  sex_col,
  allergen_col,
  sige_col
)

missing_columns <- setdiff(
  required_columns,
  names(raw_data)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing required column(s): ",
    paste(missing_columns, collapse = ", ")
  )
}


# -----------------------------------------------------------------------------
# 5. Detect suburb/locality column
# -----------------------------------------------------------------------------
#
# The exact suburb/locality field name may vary between source exports.
# The script tries common names automatically.
#
# If automatic detection fails, replace NULL below with the exact raw-data
# column name, for example:
#
# suburb_col <- "Suburb"
# -----------------------------------------------------------------------------

suburb_col <- NULL

normalise_column_name <- function(x) {
  x %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9]+", "_") %>%
    str_replace_all("^_|_$", "")
}

normalised_names <- normalise_column_name(
  names(raw_data)
)

preferred_suburb_names <- c(
  "suburb",
  "suburb_name",
  "residential_suburb",
  "patient_suburb",
  "suburb_town_city",
  "town_suburb",
  "locality",
  "locality_name",
  "residential_locality"
)

if (is.null(suburb_col)) {

  candidate_index <- which(
    normalised_names %in% preferred_suburb_names
  )

  if (length(candidate_index) == 0) {
    candidate_index <- which(
      str_detect(
        normalised_names,
        "suburb|locality"
      )
    )
  }

  if (length(candidate_index) == 0) {
    stop(
      "Could not automatically identify the suburb/locality column.\n",
      "Set `suburb_col` in Section 5 to the exact column name."
    )
  }

  suburb_col <- names(raw_data)[candidate_index[1]]
}

message(
  "Suburb/locality column used: ",
  suburb_col
)


# -----------------------------------------------------------------------------
# 6. Helper functions
# -----------------------------------------------------------------------------

clean_postcode <- function(x) {

  x <- as.character(x)

  postcode <- str_extract(
    x,
    "\\d{4}"
  )

  suppressWarnings(
    as.integer(postcode)
  )
}


clean_age <- function(x) {

  x <- as.character(x) %>%
    str_trim() %>%
    str_replace_all(",", "") %>%
    str_replace_all("[<>≤≥]", "") %>%
    str_replace_all(
      regex(
        "years|year|yrs|yr",
        ignore_case = TRUE
      ),
      ""
    )

  x[x == ""] <- NA_character_

  suppressWarnings(
    as.numeric(x)
  )
}


clean_sex <- function(x) {

  x <- str_to_lower(
    str_trim(
      as.character(x)
    )
  )

  case_when(
    x %in% c(
      "female",
      "f",
      "woman",
      "girl"
    ) ~ "Female",

    x %in% c(
      "male",
      "m",
      "man",
      "boy"
    ) ~ "Male",

    TRUE ~ NA_character_
  )
}


valid_text <- function(x) {

  x <- str_squish(
    as.character(x)
  )

  !is.na(x) &
    x != "" &
    !str_to_lower(x) %in% c(
      "na",
      "n/a",
      "none",
      "unknown",
      "missing",
      "not recorded",
      "not available"
    )
}


parse_sige_for_missingness <- function(x) {

  x <- str_squish(
    as.character(x)
  )

  x[
    str_to_lower(x) %in% c(
      "",
      "na",
      "n/a",
      "nd",
      "not done",
      "not tested",
      "missing",
      "not available"
    )
  ] <- NA_character_

  # This retains valid censored results such as <0.35 and >100.
  parse_number(x)
}


# -----------------------------------------------------------------------------
# 7. Prepare temporary cleaning variables
# -----------------------------------------------------------------------------

prepared <- raw_data %>%
  mutate(
    .patient_id = str_trim(
      as.character(
        .data[[patient_id_col]]
      )
    ),

    .postcode = clean_postcode(
      .data[[postcode_col]]
    ),

    .suburb_valid = valid_text(
      .data[[suburb_col]]
    ),

    .age = clean_age(
      .data[[age_col]]
    ),

    .sex = clean_sex(
      .data[[sex_col]]
    ),

    .allergen_valid = valid_text(
      .data[[allergen_col]]
    ),

    .sige_numeric_check = parse_sige_for_missingness(
      .data[[sige_col]]
    )
  )


# -----------------------------------------------------------------------------
# 8. Dataset-flow audit
# -----------------------------------------------------------------------------
#
# Each filter is applied sequentially so the number excluded at each step can
# be reported transparently.
# -----------------------------------------------------------------------------

flow <- tibble(
  step = "Raw data",
  test_records_remaining = nrow(prepared),
  test_records_excluded_at_step = 0L,
  unique_individuals_remaining = n_distinct(
    prepared$.patient_id[
      !is.na(prepared$.patient_id) &
        prepared$.patient_id != ""
    ]
  )
)

add_flow_step <- function(
  flow_table,
  label,
  before_data,
  after_data
) {

  bind_rows(
    flow_table,
    tibble(
      step = label,

      test_records_remaining =
        nrow(after_data),

      test_records_excluded_at_step =
        nrow(before_data) -
        nrow(after_data),

      unique_individuals_remaining =
        n_distinct(
          after_data$.patient_id[
            !is.na(after_data$.patient_id) &
              after_data$.patient_id != ""
          ]
        )
    )
  )
}


# -----------------------------------------------------------------------------
# 9. Apply eligibility and data-quality filters
# -----------------------------------------------------------------------------

current <- prepared


# 9.1 Age >=12 years
before <- current

current <- current %>%
  filter(
    !is.na(.age),
    .age >= 12
  )

flow <- add_flow_step(
  flow,
  "Age >=12 years",
  before,
  current
)


# 9.2 Postcode must be present
before <- current

current <- current %>%
  filter(
    !is.na(.postcode)
  )

flow <- add_flow_step(
  flow,
  "Residential postcode present",
  before,
  current
)


# 9.3 Queensland postcode only
before <- current

current <- current %>%
  filter(
    .postcode >= 4000,
    .postcode <= 4999
  )

flow <- add_flow_step(
  flow,
  "Queensland postcode 4000–4999",
  before,
  current
)


# 9.4 Suburb/locality must be present
before <- current

current <- current %>%
  filter(
    .suburb_valid
  )

flow <- add_flow_step(
  flow,
  "Suburb/locality present",
  before,
  current
)


# 9.5 Allergen name must be present
before <- current

current <- current %>%
  filter(
    .allergen_valid
  )

flow <- add_flow_step(
  flow,
  "Allergen information present",
  before,
  current
)


# 9.6 Usable sIgE result must be present
before <- current

current <- current %>%
  filter(
    !is.na(.sige_numeric_check)
  )

flow <- add_flow_step(
  flow,
  "Usable sIgE result present",
  before,
  current
)


# -----------------------------------------------------------------------------
# 10. Final analysis dataset
# -----------------------------------------------------------------------------

analysis_data <- current %>%
  select(
    -.patient_id,
    -.postcode,
    -.suburb_valid,
    -.age,
    -.sex,
    -.allergen_valid,
    -.sige_numeric_check
  )


# -----------------------------------------------------------------------------
# 11. Overall dataset counts
# -----------------------------------------------------------------------------

final_record_count <- nrow(current)

final_unique_individuals <- n_distinct(
  current$.patient_id[
    !is.na(current$.patient_id) &
      current$.patient_id != ""
  ]
)

overall_counts <- tibble(
  metric = c(
    "Final eligible test records",
    "Final unique individuals"
  ),

  value = c(
    final_record_count,
    final_unique_individuals
  )
)


# -----------------------------------------------------------------------------
# 12. AGE SUMMARY — test-record level
# -----------------------------------------------------------------------------
#
# This describes the age attached to each retained test record.
# -----------------------------------------------------------------------------

age_summary_test_records <- tibble(
  analytical_level = "Test records",

  n_with_age = sum(
    !is.na(current$.age)
  ),

  median_age = median(
    current$.age,
    na.rm = TRUE
  ),

  q1_age = quantile(
    current$.age,
    0.25,
    na.rm = TRUE,
    names = FALSE
  ),

  q3_age = quantile(
    current$.age,
    0.75,
    na.rm = TRUE,
    names = FALSE
  ),

  iqr_age = IQR(
    current$.age,
    na.rm = TRUE
  )
)


# -----------------------------------------------------------------------------
# 13. PARTICIPANT-LEVEL demographic dataset
# -----------------------------------------------------------------------------
#
# Participants may have repeated tests. For participant-level descriptive
# summaries:
#
#   age = median of that participant's recorded ages
#   sex = first non-missing cleaned sex value
#
# An audit of participants with inconsistent recorded sex is also produced.
# -----------------------------------------------------------------------------

participant_demographics <- current %>%
  filter(
    !is.na(.patient_id),
    .patient_id != ""
  ) %>%
  group_by(
    .patient_id
  ) %>%
  summarise(
    participant_age = median(
      .age,
      na.rm = TRUE
    ),

    participant_sex = {
      valid_sex <- .sex[
        !is.na(.sex)
      ]

      if (length(valid_sex) == 0) {
        NA_character_
      } else {
        valid_sex[1]
      }
    },

    n_distinct_recorded_sexes =
      n_distinct(
        .sex[
          !is.na(.sex)
        ]
      ),

    n_test_records = n(),

    .groups = "drop"
  )


# -----------------------------------------------------------------------------
# 14. AGE SUMMARY — participant level
# -----------------------------------------------------------------------------

age_summary_participants <- tibble(
  analytical_level = "Unique individuals",

  n_with_age = sum(
    !is.na(
      participant_demographics$participant_age
    )
  ),

  median_age = median(
    participant_demographics$participant_age,
    na.rm = TRUE
  ),

  q1_age = quantile(
    participant_demographics$participant_age,
    0.25,
    na.rm = TRUE,
    names = FALSE
  ),

  q3_age = quantile(
    participant_demographics$participant_age,
    0.75,
    na.rm = TRUE,
    names = FALSE
  ),

  iqr_age = IQR(
    participant_demographics$participant_age,
    na.rm = TRUE
  )
)


age_summary <- bind_rows(
  age_summary_test_records,
  age_summary_participants
)


# -----------------------------------------------------------------------------
# 15. SEX DISTRIBUTION — test-record level
# -----------------------------------------------------------------------------

sex_distribution_test_records <- current %>%
  mutate(
    sex_summary = if_else(
      is.na(.sex),
      "Missing/other",
      .sex
    )
  ) %>%
  count(
    sex_summary,
    name = "n_test_records"
  ) %>%
  mutate(
    percent_test_records =
      100 *
      n_test_records /
      sum(n_test_records)
  ) %>%
  arrange(
    factor(
      sex_summary,
      levels = c(
        "Female",
        "Male",
        "Missing/other"
      )
    )
  )


# -----------------------------------------------------------------------------
# 16. SEX DISTRIBUTION — participant level
# -----------------------------------------------------------------------------

sex_distribution_participants <- participant_demographics %>%
  mutate(
    sex_summary = if_else(
      is.na(participant_sex),
      "Missing/other",
      participant_sex
    )
  ) %>%
  count(
    sex_summary,
    name = "n_individuals"
  ) %>%
  mutate(
    percent_individuals =
      100 *
      n_individuals /
      sum(n_individuals)
  ) %>%
  arrange(
    factor(
      sex_summary,
      levels = c(
        "Female",
        "Male",
        "Missing/other"
      )
    )
  )


# -----------------------------------------------------------------------------
# 17. AGE-GROUP DISTRIBUTION
# -----------------------------------------------------------------------------

age_group_levels <- c(
  "12–17",
  "18–29",
  "30–44",
  "45–64",
  "≥65"
)

current <- current %>%
  mutate(
    .age_group = case_when(
      .age >= 12 & .age <= 17 ~ "12–17",
      .age >= 18 & .age <= 29 ~ "18–29",
      .age >= 30 & .age <= 44 ~ "30–44",
      .age >= 45 & .age <= 64 ~ "45–64",
      .age >= 65 ~ "≥65",
      TRUE ~ NA_character_
    ),

    .age_group = factor(
      .age_group,
      levels = age_group_levels
    )
  )


age_group_distribution_test_records <- current %>%
  count(
    .age_group,
    name = "n_test_records",
    .drop = FALSE
  ) %>%
  mutate(
    percent_test_records =
      100 *
      n_test_records /
      sum(n_test_records)
  ) %>%
  rename(
    age_group = .age_group
  )


participant_age_groups <- participant_demographics %>%
  mutate(
    age_group = case_when(
      participant_age >= 12 &
        participant_age <= 17 ~ "12–17",

      participant_age >= 18 &
        participant_age <= 29 ~ "18–29",

      participant_age >= 30 &
        participant_age <= 44 ~ "30–44",

      participant_age >= 45 &
        participant_age <= 64 ~ "45–64",

      participant_age >= 65 ~ "≥65",

      TRUE ~ NA_character_
    ),

    age_group = factor(
      age_group,
      levels = age_group_levels
    )
  )


age_group_distribution_participants <- participant_age_groups %>%
  count(
    age_group,
    name = "n_individuals",
    .drop = FALSE
  ) %>%
  mutate(
    percent_individuals =
      100 *
      n_individuals /
      sum(n_individuals)
  )


# -----------------------------------------------------------------------------
# 18. Demographic consistency audit
# -----------------------------------------------------------------------------

sex_inconsistency_audit <- participant_demographics %>%
  filter(
    n_distinct_recorded_sexes > 1
  )

message(
  "Participants with inconsistent non-missing sex entries: ",
  nrow(sex_inconsistency_audit)
)


# -----------------------------------------------------------------------------
# 19. Save the filtered dataset LOCALLY
# -----------------------------------------------------------------------------
#
# WARNING:
# This is individual-level pathology data. Do not commit this file to GitHub.
# -----------------------------------------------------------------------------

openxlsx::write.xlsx(
  analysis_data,
  file = output_path,
  overwrite = TRUE
)

message(
  "Filtered Queensland dataset saved locally to: ",
  output_path
)


# -----------------------------------------------------------------------------
# 20. Save aggregate/audit outputs
# -----------------------------------------------------------------------------

write_csv(
  flow,
  file.path(
    output_dir,
    "01_dataset_flow.csv"
  )
)

write_csv(
  overall_counts,
  file.path(
    output_dir,
    "02_overall_counts.csv"
  )
)

write_csv(
  age_summary,
  file.path(
    output_dir,
    "03_age_summary.csv"
  )
)

write_csv(
  sex_distribution_test_records,
  file.path(
    output_dir,
    "04_sex_distribution_test_records.csv"
  )
)

write_csv(
  sex_distribution_participants,
  file.path(
    output_dir,
    "05_sex_distribution_unique_individuals.csv"
  )
)

write_csv(
  age_group_distribution_test_records,
  file.path(
    output_dir,
    "06_age_group_distribution_test_records.csv"
  )
)

write_csv(
  age_group_distribution_participants,
  file.path(
    output_dir,
    "07_age_group_distribution_unique_individuals.csv"
  )
)

write_csv(
  sex_inconsistency_audit,
  file.path(
    output_dir,
    "08_participant_sex_consistency_audit.csv"
  )
)


# -----------------------------------------------------------------------------
# 21. Write a readable cohort summary
# -----------------------------------------------------------------------------

female_test_records <- sex_distribution_test_records %>%
  filter(
    sex_summary == "Female"
  )

male_test_records <- sex_distribution_test_records %>%
  filter(
    sex_summary == "Male"
  )

participant_age_row <- age_summary %>%
  filter(
    analytical_level == "Unique individuals"
  )

test_record_age_row <- age_summary %>%
  filter(
    analytical_level == "Test records"
  )


summary_lines <- c(
  "VIATAR_QLD — DATA PREPARATION AND COHORT SUMMARY",
  "================================================",
  "",
  paste0(
    "Final eligible test records: ",
    format(
      final_record_count,
      big.mark = ","
    )
  ),
  paste0(
    "Final unique individuals: ",
    format(
      final_unique_individuals,
      big.mark = ","
    )
  ),
  "",
  "Age summary — unique individuals:",
  paste0(
    "  Median age: ",
    round(
      participant_age_row$median_age,
      2
    ),
    " years"
  ),
  paste0(
    "  Q1: ",
    round(
      participant_age_row$q1_age,
      2
    ),
    " years"
  ),
  paste0(
    "  Q3: ",
    round(
      participant_age_row$q3_age,
      2
    ),
    " years"
  ),
  paste0(
    "  IQR: ",
    round(
      participant_age_row$iqr_age,
      2
    ),
    " years"
  ),
  "",
  "Age summary — test records:",
  paste0(
    "  Median age: ",
    round(
      test_record_age_row$median_age,
      2
    ),
    " years"
  ),
  paste0(
    "  Q1: ",
    round(
      test_record_age_row$q1_age,
      2
    ),
    " years"
  ),
  paste0(
    "  Q3: ",
    round(
      test_record_age_row$q3_age,
      2
    ),
    " years"
  ),
  paste0(
    "  IQR: ",
    round(
      test_record_age_row$iqr_age,
      2
    ),
    " years"
  ),
  "",
  "Sex distribution — test records:",
  if (nrow(female_test_records) == 1) {
    paste0(
      "  Female: ",
      format(
        female_test_records$n_test_records,
        big.mark = ","
      ),
      " (",
      round(
        female_test_records$percent_test_records,
        2
      ),
      "%)"
    )
  } else {
    "  Female: not available"
  },
  if (nrow(male_test_records) == 1) {
    paste0(
      "  Male: ",
      format(
        male_test_records$n_test_records,
        big.mark = ","
      ),
      " (",
      round(
        male_test_records$percent_test_records,
        2
      ),
      "%)"
    )
  } else {
    "  Male: not available"
  },
  "",
  "Filtering rules:",
  "  - age >=12 years",
  "  - residential postcode present",
  "  - Queensland postcode 4000–4999",
  "  - suburb/locality present",
  "  - allergen information present",
  "  - usable sIgE result present"
)

writeLines(
  summary_lines,
  con = file.path(
    output_dir,
    "cohort_summary.txt"
  )
)

cat(
  paste(
    summary_lines,
    collapse = "\n"
  ),
  "\n"
)


# -----------------------------------------------------------------------------
# 22. Session information
# -----------------------------------------------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    output_dir,
    "data_preparation_sessionInfo.txt"
  )
)
