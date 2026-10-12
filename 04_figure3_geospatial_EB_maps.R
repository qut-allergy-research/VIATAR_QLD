# =============================================================================
# VIATAR_QLD
# Figure 3. Geographic variation in aeroallergen sensitisation frequency
#           and level across Queensland, Australia
#
# PURPOSE
# Reproduce the two postcode-level map panels shown in Figure 3:
#
#   Panel A
#   Empirical Bayes (EB)-smoothed proportion of sIgE test records positive
#   for each aeroallergen group.
#
#   Panel B
#   EB-smoothed geometric mean sIgE concentration among positive test records
#   (sIgE >= 0.35 kU/L). Postcode-level mean log(sIgE) values are smoothed
#   and then back-transformed to the original kU/L scale.
#
# MAIN EB PRIOR WEIGHT
#   k = 20
#
# MAP MASKING
#   Grey = postcode has test data for the aeroallergen group but fewer than
#          5 positive test records.
#   Blue = no available test records for that aeroallergen group/postcode.
#
# ANALYTICAL UNIT
# Each sIgE test record is retained as a separate observation. Spatial
# summaries are aggregated at postcode x aeroallergen-group level.
#
# DATA PRIVACY
# The patient-level pathology dataset is NOT included in this public
# repository. Keep all real data outside GitHub.
#
# SPATIAL DATA
# Uses Australian Statistical Geography Standard (ASGS) Edition 3,
# 2021 Postal Areas (POA) boundaries.
#
# IMPORTANT
# The original working files contained several exploratory mapping versions:
# individual-allergen maps, unsmoothed maps, ordinal "severity" maps, arithmetic
# mean sIgE maps, and multiple EB k values. Those exploratory versions are not
# used here. This script matches the current manuscript Figure 3 and legend:
#   - seven aeroallergen GROUPS, not individual allergens
#   - fixed k = 20 for the displayed maps
#   - n < 5 based on POSITIVE test records
#   - geometric mean sIgE after smoothing mean log(sIgE)
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Packages
# -----------------------------------------------------------------------------

library(readxl)
library(dplyr)
library(stringr)
library(readr)
library(tidyr)
library(sf)
library(ggplot2)
library(RColorBrewer)
library(patchwork)


# -----------------------------------------------------------------------------
# 2. Local paths
# -----------------------------------------------------------------------------
#
# Edit these paths locally if needed.
# Keep the real pathology dataset and spatial files outside the public repo,
# or place them in a local data/ directory that is excluded by .gitignore.
# -----------------------------------------------------------------------------

data_path <- file.path(
  "data",
  "VIATAR_filtered_QLD_01122025.xlsx"
)

shapefile_path <- file.path(
  "data",
  "spatial",
  "POA_2021_AUST_GDA2020.shp"
)

output_dir <- file.path("outputs", "figure3")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "figures"), recursive = TRUE, showWarnings = FALSE)

if (!file.exists(data_path)) {
  stop(
    "Pathology data file not found: ", data_path,
    "\nUpdate 'data_path' locally. Do not upload the real dataset to GitHub."
  )
}

if (!file.exists(shapefile_path)) {
  stop(
    "POA shapefile not found: ", shapefile_path,
    "\nUpdate 'shapefile_path' locally."
  )
}


# -----------------------------------------------------------------------------
# 3. Analysis parameters
# -----------------------------------------------------------------------------

positive_threshold <- 0.35
small_n_threshold  <- 5
k_main             <- 20

# Sensitivity values described in the manuscript.
# They are calculated and exported, but only k = 20 is plotted in Figure 3.
k_sensitivity <- c(10, 30, 50)

# Figure 3A frequency bins used in the final map.
frequency_breaks <- c(
  0, 7, 14, 22, 28, 36, 43, 50,
  56, 64, 72, 79, 86, 93, Inf
)

frequency_labels <- c(
  "0–<7%",
  "7–<14%",
  "14–<22%",
  "22–<28%",
  "28–<36%",
  "36–<43%",
  "43–<50%",
  "50–<56%",
  "56–<64%",
  "64–<72%",
  "72–<79%",
  "79–<86%",
  "86–<93%",
  "93–100%"
)

frequency_colours <- colorRampPalette(
  brewer.pal(9, "YlOrRd")
)(length(frequency_labels))

# Figure 3B categories for the EB-smoothed geometric mean sIgE.
sige_level_labels <- c(
  "0.35–<0.70",
  "0.70–<3.50",
  "3.50–<17.50",
  "17.50–<50.00",
  "50.00–100.00",
  ">100.00"
)

sige_level_colours <- colorRampPalette(
  brewer.pal(9, "YlOrRd")
)(length(sige_level_labels))

# Colours described in the Figure 3 legend.
small_n_colour <- "grey80"
no_data_colour <- "#cfe8f3"


# -----------------------------------------------------------------------------
# 4. Exact aeroallergen-group lookup
# -----------------------------------------------------------------------------
#
# The lookup follows the aeroallergen groups defined in the manuscript.
# A small number of spelling variants are included only to accommodate labels
# observed in the working analysis files.
# -----------------------------------------------------------------------------

allergen_lookup <- tibble::tribble(
  ~allergen_key,                         ~allergen_group,

  # House dust mites
  "dermatophagoides farinae",            "Dust mites",
  "dermatophagoides pteronyssinus",      "Dust mites",
  "blomia tropicalis",                   "Dust mites",

  # Moulds
  "alternaria alternata",                "Mould",
  "aspergillus fumigatus",               "Mould",
  "cladosporium herbarum",               "Mould",
  "mould mix 1",                         "Mould",

  # Grass pollens
  "bahia grass (paspalum)",              "Grass pollen",
  "bermuda (couch) grass",               "Grass pollen",
  "johnson grass (sorghum)",             "Grass pollen",
  "perennial rye grass",                 "Grass pollen",
  "grass mix 2",                         "Grass pollen",

  # Weed pollens
  "common ragweed",                      "Weed pollen",
  "weed mix 1",                          "Weed pollen",

  # Tree pollens
  "eucalyptus",                          "Tree pollen",
  "gum tree",                            "Tree pollen",
  "acacia",                              "Tree pollen",
  "tree mix 7",                          "Tree pollen",

  # Animal-derived allergens
  "cat epithelium & dander",             "Animal dander",
  "cat epithelium and dander",           "Animal dander",
  "dog dander",                          "Animal dander",
  "chicken feathers",                    "Animal dander",
  "animal mix 1",                        "Animal dander",
  "feather mix",                         "Animal dander",

  # Cockroach
  "cockroach (blatella germanica)",      "Cockroach",
  "cockroach (blattella germanica)",     "Cockroach",
  "blattella germanica",                 "Cockroach"
)

group_order <- c(
  "Dust mites",
  "Animal dander",
  "Cockroach",
  "Mould",
  "Tree pollen",
  "Weed pollen",
  "Grass pollen"
)


# -----------------------------------------------------------------------------
# 5. Read pathology data
# -----------------------------------------------------------------------------

raw_data <- read_xlsx(
  data_path,
  guess_max = 20000
) %>%
  as_tibble()


# -----------------------------------------------------------------------------
# 6. Prepare pathology data
# -----------------------------------------------------------------------------
#
# IMPORTANT:
# If a result is stored as "<0.35", simply using as.numeric() would return NA.
# That would remove a negative test from the Panel A denominator.
#
# We therefore:
#   - parse the numeric boundary,
#   - keep an explicit "<" flag,
#   - classify "<0.35" as negative.
#
# For concentration analyses, only positive tests are retained.
# -----------------------------------------------------------------------------

analysis_data <- raw_data %>%
  transmute(
    postcode = suppressWarnings(
      as.integer(as.character(`Post.Code`))
    ),
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
  left_join(allergen_lookup, by = "allergen_key") %>%
  filter(
    !is.na(postcode),
    postcode >= 4000,
    postcode <= 4999,
    !is.na(allergen_group),
    !is.na(sige)
  ) %>%
  mutate(
    positive = if_else(
      below_limit,
      0L,
      as.integer(sige >= positive_threshold)
    ),
    allergen_group = factor(
      allergen_group,
      levels = group_order
    )
  )


# -----------------------------------------------------------------------------
# 7. Audit allergen mapping
# -----------------------------------------------------------------------------

unmapped_allergens <- raw_data %>%
  transmute(
    allergen_raw = as.character(`Test.Name..Allergen.`),
    allergen_key = allergen_raw %>%
      str_to_lower() %>%
      str_squish()
  ) %>%
  distinct(allergen_raw, allergen_key) %>%
  anti_join(allergen_lookup, by = "allergen_key")

write_csv(
  unmapped_allergens,
  file.path(
    output_dir,
    "tables",
    "Figure3_unmapped_allergen_labels.csv"
  )
)

group_record_counts <- analysis_data %>%
  count(allergen_group, name = "n_test_records")

write_csv(
  group_record_counts,
  file.path(
    output_dir,
    "tables",
    "Figure3_group_test_record_counts.csv"
  )
)


# -----------------------------------------------------------------------------
# 8. Read and prepare Queensland postcode polygons
# -----------------------------------------------------------------------------

poa <- st_read(
  shapefile_path,
  quiet = TRUE
)

possible_postcode_columns <- c(
  "POA_CODE21",
  "POA_CODE_2021",
  "POA21_CODE",
  "POA_CODE",
  "POA_CODE_2016"
)

postcode_column <- intersect(
  possible_postcode_columns,
  names(poa)
)

if (length(postcode_column) == 0) {
  stop(
    "Could not identify the postcode code column in the POA shapefile."
  )
}

postcode_column <- postcode_column[1]

poa_qld <- poa %>%
  rename(
    postcode_raw = all_of(postcode_column)
  ) %>%
  mutate(
    postcode = suppressWarnings(
      as.integer(as.character(postcode_raw))
    )
  ) %>%
  filter(
    !is.na(postcode),
    postcode >= 4000,
    postcode <= 4999
  ) %>%
  select(postcode, geometry)


# =============================================================================
# PANEL A
# EB-SMOOTHED SENSITISATION FREQUENCY
# =============================================================================


# -----------------------------------------------------------------------------
# 9. Raw postcode x group frequency summaries
# -----------------------------------------------------------------------------

frequency_summary <- analysis_data %>%
  group_by(
    postcode,
    allergen_group
  ) %>%
  summarise(
    n_total = n(),
    n_positive = sum(positive),
    proportion_raw = n_positive / n_total,
    .groups = "drop"
  )


# -----------------------------------------------------------------------------
# 10. Global prior frequency for each aeroallergen group
# -----------------------------------------------------------------------------
#
# The prior mean is the overall positive proportion for that aeroallergen group.
# -----------------------------------------------------------------------------

frequency_priors <- frequency_summary %>%
  group_by(allergen_group) %>%
  summarise(
    prior_frequency = sum(n_positive) / sum(n_total),
    .groups = "drop"
  )


# -----------------------------------------------------------------------------
# 11. EB smoothing for frequency
# -----------------------------------------------------------------------------
#
# p_EB = (x + k * m) / (n + k)
#
# x = number of positive test records in postcode
# n = total test records in postcode
# m = global positive proportion for that aeroallergen group
# k = prior weight
# -----------------------------------------------------------------------------

frequency_eb <- frequency_summary %>%
  left_join(
    frequency_priors,
    by = "allergen_group"
  ) %>%
  mutate(
    proportion_eb_k20 =
      (n_positive + k_main * prior_frequency) /
      (n_total + k_main),

    percent_eb_k20 = 100 * proportion_eb_k20,

    insufficient_positive_n = n_positive < small_n_threshold
  )

# Optional sensitivity estimates: k = 10, 30, 50
for (k in k_sensitivity) {
  frequency_eb[[paste0("percent_eb_k", k)]] <-
    100 * (
      frequency_eb$n_positive +
        k * frequency_eb$prior_frequency
    ) / (
      frequency_eb$n_total + k
    )
}


# -----------------------------------------------------------------------------
# 12. Create complete postcode x group map data
# -----------------------------------------------------------------------------

all_groups <- tibble(
  allergen_group = factor(
    group_order,
    levels = group_order
  )
)

map_template <- crossing(
  postcode = poa_qld$postcode,
  allergen_group = all_groups$allergen_group
)

frequency_complete <- map_template %>%
  left_join(
    frequency_eb,
    by = c("postcode", "allergen_group")
  ) %>%
  mutate(
    no_data = is.na(n_total),

    # A postcode with tests but 0–4 positive records is grey.
    insufficient_positive_n = if_else(
      no_data,
      FALSE,
      replace_na(insufficient_positive_n, FALSE)
    ),

    frequency_bin = cut(
      percent_eb_k20,
      breaks = frequency_breaks,
      labels = frequency_labels,
      include.lowest = TRUE,
      right = FALSE
    )
  )

frequency_map <- poa_qld %>%
  left_join(
    frequency_complete,
    by = "postcode"
  )


# -----------------------------------------------------------------------------
# 13. Save Panel A postcode summaries
# -----------------------------------------------------------------------------

write_csv(
  frequency_complete,
  file.path(
    output_dir,
    "tables",
    "Figure3A_postcode_EB_frequency.csv"
  )
)


# -----------------------------------------------------------------------------
# 14. Plot Panel A
# -----------------------------------------------------------------------------

p_frequency <- ggplot() +

  # Valid mapped values: >= 5 positive test records
  geom_sf(
    data = frequency_map %>%
      filter(
        !no_data,
        !insufficient_positive_n
      ),
    aes(fill = frequency_bin),
    colour = NA
  ) +

  # Insufficient positive test records
  geom_sf(
    data = frequency_map %>%
      filter(
        !no_data,
        insufficient_positive_n
      ),
    fill = small_n_colour,
    colour = NA
  ) +

  # No test records
  geom_sf(
    data = frequency_map %>%
      filter(no_data),
    fill = no_data_colour,
    colour = NA
  ) +

  facet_wrap(
    ~ allergen_group,
    ncol = 4
  ) +

  scale_fill_manual(
    values = setNames(
      frequency_colours,
      frequency_labels
    ),
    drop = FALSE,
    name = "Percent sensitised\n(EB k = 20)"
  ) +

  coord_sf(
    expand = FALSE
  ) +

  theme_void(base_size = 11) +

  theme(
    strip.text = element_text(
      face = "plain",
      size = 10
    ),
    legend.position = "right",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    panel.spacing = unit(0.5, "lines")
  )


# =============================================================================
# PANEL B
# EB-SMOOTHED GEOMETRIC MEAN sIgE AMONG POSITIVE TESTS
# =============================================================================


# -----------------------------------------------------------------------------
# 15. Positive tests only
# -----------------------------------------------------------------------------

positive_data <- analysis_data %>%
  filter(positive == 1) %>%
  mutate(
    log_sige = log(sige)
  )


# -----------------------------------------------------------------------------
# 16. Postcode x group mean log(sIgE)
# -----------------------------------------------------------------------------

level_summary <- positive_data %>%
  group_by(
    postcode,
    allergen_group
  ) %>%
  summarise(
    n_positive = n(),
    sum_log_sige = sum(log_sige),
    mean_log_sige_raw = mean(log_sige),
    geometric_mean_raw = exp(mean_log_sige_raw),
    .groups = "drop"
  )


# -----------------------------------------------------------------------------
# 17. Global prior mean log(sIgE) for each aeroallergen group
# -----------------------------------------------------------------------------

level_priors <- positive_data %>%
  group_by(allergen_group) %>%
  summarise(
    prior_mean_log_sige = mean(log_sige),
    .groups = "drop"
  )


# -----------------------------------------------------------------------------
# 18. EB smoothing on the LOG scale, then back-transform
# -----------------------------------------------------------------------------
#
# mean_log_EB = (sum(log(sIgE)) + k * global_mean_log) / (n_positive + k)
#
# geometric_mean_EB = exp(mean_log_EB)
# -----------------------------------------------------------------------------

level_eb <- level_summary %>%
  left_join(
    level_priors,
    by = "allergen_group"
  ) %>%
  mutate(
    mean_log_sige_eb_k20 =
      (
        sum_log_sige +
          k_main * prior_mean_log_sige
      ) /
      (
        n_positive + k_main
      ),

    geometric_mean_eb_k20 =
      exp(mean_log_sige_eb_k20)
  )

# Optional sensitivity estimates: k = 10, 30, 50
for (k in k_sensitivity) {

  mean_log_col <- paste0("mean_log_sige_eb_k", k)
  geom_col <- paste0("geometric_mean_eb_k", k)

  level_eb[[mean_log_col]] <-
    (
      level_eb$sum_log_sige +
        k * level_eb$prior_mean_log_sige
    ) /
    (
      level_eb$n_positive + k
    )

  level_eb[[geom_col]] <-
    exp(level_eb[[mean_log_col]])
}


# -----------------------------------------------------------------------------
# 19. Preserve total-test information so blue and grey have distinct meanings
# -----------------------------------------------------------------------------

total_test_counts <- frequency_summary %>%
  select(
    postcode,
    allergen_group,
    n_total
  )

level_complete <- map_template %>%
  left_join(
    total_test_counts,
    by = c("postcode", "allergen_group")
  ) %>%
  left_join(
    level_eb,
    by = c("postcode", "allergen_group")
  ) %>%
  mutate(
    no_data = is.na(n_total),

    n_positive = replace_na(
      n_positive,
      0L
    ),

    # Grey if there are group test records but fewer than 5 positive tests.
    insufficient_positive_n =
      !no_data &
      n_positive < small_n_threshold,

    sige_level_bin = case_when(
      is.na(geometric_mean_eb_k20) ~ NA_character_,
      geometric_mean_eb_k20 < 0.70 ~ "0.35–<0.70",
      geometric_mean_eb_k20 < 3.50 ~ "0.70–<3.50",
      geometric_mean_eb_k20 < 17.50 ~ "3.50–<17.50",
      geometric_mean_eb_k20 < 50.00 ~ "17.50–<50.00",
      geometric_mean_eb_k20 <= 100.00 ~ "50.00–100.00",
      geometric_mean_eb_k20 > 100.00 ~ ">100.00"
    ),

    sige_level_bin = factor(
      sige_level_bin,
      levels = sige_level_labels
    )
  )

level_map <- poa_qld %>%
  left_join(
    level_complete,
    by = "postcode"
  )


# -----------------------------------------------------------------------------
# 20. Save Panel B postcode summaries
# -----------------------------------------------------------------------------

write_csv(
  level_complete,
  file.path(
    output_dir,
    "tables",
    "Figure3B_postcode_EB_geometric_mean_sIgE.csv"
  )
)


# -----------------------------------------------------------------------------
# 21. Plot Panel B
# -----------------------------------------------------------------------------

p_level <- ggplot() +

  # Valid mapped values: >= 5 positive test records
  geom_sf(
    data = level_map %>%
      filter(
        !no_data,
        !insufficient_positive_n
      ),
    aes(fill = sige_level_bin),
    colour = NA
  ) +

  # Insufficient positive test records
  geom_sf(
    data = level_map %>%
      filter(
        !no_data,
        insufficient_positive_n
      ),
    fill = small_n_colour,
    colour = NA
  ) +

  # No test records
  geom_sf(
    data = level_map %>%
      filter(no_data),
    fill = no_data_colour,
    colour = NA
  ) +

  facet_wrap(
    ~ allergen_group,
    ncol = 4
  ) +

  scale_fill_manual(
    values = setNames(
      sige_level_colours,
      sige_level_labels
    ),
    drop = FALSE,
    name = "EB-smoothed geometric\nmean sIgE (kU/L)"
  ) +

  coord_sf(
    expand = FALSE
  ) +

  theme_void(base_size = 11) +

  theme(
    strip.text = element_text(
      face = "plain",
      size = 10
    ),
    legend.position = "right",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    panel.spacing = unit(0.5, "lines")
  )


# =============================================================================
# 22. Combine Figure 3A and Figure 3B
# =============================================================================

figure3 <- (
  p_frequency /
    p_level
) +
  plot_annotation(
    tag_levels = "A"
  ) &
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 13
    )
  )

print(figure3)

ggsave(
  filename = file.path(
    output_dir,
    "figures",
    "Figure3_geographic_EB_frequency_and_sIgE.png"
  ),
  plot = figure3,
  width = 12,
  height = 14,
  dpi = 600
)


# -----------------------------------------------------------------------------
# 23. Session information
# -----------------------------------------------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    output_dir,
    "Figure3_sessionInfo.txt"
  )
)
