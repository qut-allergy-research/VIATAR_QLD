# =============================================================================
# VIATAR_QLD
# Figure 1. Sensitisation frequency and sIgE levels across aeroallergen groups
#
# Panel A: sensitisation frequency across aeroallergen groups
# Panel B: distribution of sIgE concentrations among sensitised tests
#
# IMPORTANT:
# - Raw/de-identified patient data are NOT included in the public repository.
# - Panel A below currently reproduces the published/reporting plot from the
#   final summary values in the original analysis script. The raw-data code
#   that generated those percentages and pairwise frequency comparisons should
#   replace that section when available.
# - Panel B is generated from test-level sIgE records and includes the
#   Kruskal-Wallis test plus BH-adjusted pairwise Dunn tests.
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Packages
# -----------------------------------------------------------------------------

library(readxl)
library(dplyr)
library(stringr)
library(readr)
library(janitor)
library(ggplot2)
library(ggsignif)
library(rstatix)
library(ggpubr)


# -----------------------------------------------------------------------------
# 2. Output directories
# -----------------------------------------------------------------------------

dir.create("outputs", showWarnings = FALSE)
dir.create(file.path("outputs", "figures"), showWarnings = FALSE)
dir.create(file.path("outputs", "tables"), showWarnings = FALSE)


# -----------------------------------------------------------------------------
# 3. Common figure settings
# -----------------------------------------------------------------------------

group_levels <- c(
  "HDM",
  "Mould",
  "Grass Pollen",
  "Weed Pollen",
  "Tree Pollen",
  "Animal Dander",
  "Cockroach"
)

group_colours <- c(
  "HDM"           = "#1F77B4",
  "Mould"         = "#2CA02C",
  "Grass Pollen"  = "#FF7F0E",
  "Weed Pollen"   = "#8B4513",
  "Tree Pollen"   = "#9467BD",
  "Animal Dander" = "#D62728",
  "Cockroach"     = "#C51B7D"
)

sige_breaks <- c(0.35, 0.69, 3.49, 17.49, 52.49, 100)


# =============================================================================
# PANEL A — SENSITISATION FREQUENCY
# =============================================================================
#
# NOTE:
# The original plotting script entered these summary values manually.
# This section therefore reproduces the final Figure 1A plot, but it does NOT
# independently calculate the percentages from raw data.
#
# When the original frequency-analysis code is located, replace this table with
# code that calculates n tested, n sensitised, percentages, and BH-adjusted
# pairwise comparisons directly from the analysis dataset.
# =============================================================================

figure1a_summary <- tibble::tibble(
  group = factor(group_levels, levels = group_levels),
  percent = c(42.45, 11.41, 20.71, 11.39, 9.11, 17.79, 16.22),
  n_reported = c(129826, 55734, 148506, 18637, 29714, 126017, 129826)
)

# Non-significant pairwise comparisons shown in the final figure.
# These were hard-coded in the original plotting script.
figure1a_ns <- list(
  c("Tree Pollen", "Weed Pollen"),
  c("Cockroach", "Tree Pollen"),
  c("Cockroach", "Weed Pollen")
)

p1a <- ggplot(
  figure1a_summary,
  aes(x = group, y = percent, fill = group)
) +
  geom_col(width = 0.72) +
  geom_text(
    aes(label = sprintf("%.1f%%", percent)),
    vjust = -0.35,
    size = 4.5,
    fontface = "bold"
  ) +
  geom_text(
    aes(y = 2.5, label = paste0("n=", scales::comma(n_reported))),
    colour = "white",
    size = 3.5
  ) +
  scale_fill_manual(values = group_colours) +
  scale_y_continuous(
    limits = c(0, 75),
    breaks = seq(0, 70, 10),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = NULL,
    y = "Proportion sensitised (%)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold")
  ) +
  geom_signif(
    comparisons = figure1a_ns,
    annotations = rep("ns", length(figure1a_ns)),
    y_position = c(25, 32, 38),
    tip_length = 0.01,
    textsize = 4
  )

print(p1a)

ggsave(
  filename = file.path("outputs", "figures", "Figure1A_sensitisation_frequency.png"),
  plot = p1a,
  width = 8.5,
  height = 5.5,
  dpi = 300
)


# =============================================================================
# PANEL B — sIgE LEVELS AMONG SENSITISED TESTS
# =============================================================================

# -----------------------------------------------------------------------------
# 4. Data path
# -----------------------------------------------------------------------------
#
# Keep the real dataset outside the public GitHub repository.
# This relative path is only a local placeholder.
# Add data/ to .gitignore before placing any real data there.
# -----------------------------------------------------------------------------

data_path <- file.path("data", "VIATAR_filtered_QLD.xlsx")

if (!file.exists(data_path)) {
  stop(
    "Data file not found at: ", data_path,
    "\nKeep the real dataset locally and do not upload it to GitHub."
  )
}


# -----------------------------------------------------------------------------
# 5. Read and prepare the variables required for Figure 1B
# -----------------------------------------------------------------------------

raw_data <- read_excel(data_path) %>%
  clean_names()

figure1b_data <- raw_data %>%
  transmute(
    allergen = test_name_allergen,
    sige_raw = test_value_spt_or_sp_ig_e_values_of_each_allergen_or_mm_d_for_spt
  ) %>%
  mutate(
    sige_text = str_trim(as.character(sige_raw)),
    below_detection = str_detect(sige_text, "^<"),
    sige_detected = parse_number(sige_text),

    # Preserves the approach used in the original analysis code:
    # values reported as "<x" are assigned x/2.
    # Values such as ">100", if present, are parsed to their numeric boundary.
    sige = if_else(
      below_detection,
      sige_detected / 2,
      sige_detected
    ),

    allergen_clean = allergen %>%
      as.character() %>%
      str_to_lower() %>%
      str_replace_all("_|\\.", " ") %>%
      str_squish(),

    allergen_group = case_when(
      str_detect(allergen_clean, "dermatophagoides|blomia") ~ "HDM",

      str_detect(
        allergen_clean,
        "alternaria|aspergillus|cladosporium|mould"
      ) ~ "Mould",

      str_detect(
        allergen_clean,
        "bahia|bermuda|couch|johnson|perennial rye|rye grass|grass mix"
      ) ~ "Grass Pollen",

      str_detect(
        allergen_clean,
        "ragweed|weed mix"
      ) ~ "Weed Pollen",

      str_detect(
        allergen_clean,
        "eucalyptus|gum tree|acacia|tree mix"
      ) ~ "Tree Pollen",

      str_detect(
        allergen_clean,
        "cat|dog|chicken|animal mix|feather"
      ) ~ "Animal Dander",

      str_detect(
        allergen_clean,
        "cockroach|blattella"
      ) ~ "Cockroach",

      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    !is.na(allergen_group),
    !is.na(sige)
  ) %>%
  mutate(
    allergen_group = factor(allergen_group, levels = group_levels)
  )


# -----------------------------------------------------------------------------
# 6. Restrict Figure 1B to sensitised tests
# -----------------------------------------------------------------------------

positive_tests <- figure1b_data %>%
  filter(sige >= 0.35)


# -----------------------------------------------------------------------------
# 7. Descriptive statistics: median and IQR
# -----------------------------------------------------------------------------

sige_summary <- positive_tests %>%
  group_by(allergen_group) %>%
  summarise(
    n_tests = n(),
    median_sige = median(sige, na.rm = TRUE),
    q1_sige = quantile(sige, 0.25, na.rm = TRUE),
    q3_sige = quantile(sige, 0.75, na.rm = TRUE),
    iqr_sige = IQR(sige, na.rm = TRUE),
    .groups = "drop"
  )

print(sige_summary)

write_csv(
  sige_summary,
  file.path("outputs", "tables", "Figure1B_sIgE_summary_statistics.csv")
)


# -----------------------------------------------------------------------------
# 8. Overall comparison: Kruskal-Wallis test
# -----------------------------------------------------------------------------

kruskal_result <- positive_tests %>%
  kruskal_test(sige ~ allergen_group)

print(kruskal_result)

write_csv(
  kruskal_result,
  file.path("outputs", "tables", "Figure1B_kruskal_wallis.csv")
)


# -----------------------------------------------------------------------------
# 9. Post-hoc pairwise Dunn tests with Benjamini-Hochberg correction
# -----------------------------------------------------------------------------

dunn_results <- positive_tests %>%
  dunn_test(
    sige ~ allergen_group,
    p.adjust.method = "BH"
  )

print(dunn_results)

write_csv(
  dunn_results,
  file.path("outputs", "tables", "Figure1B_dunn_BH_pairwise.csv")
)


# -----------------------------------------------------------------------------
# 10. Keep only non-significant pairwise comparisons for annotation
# -----------------------------------------------------------------------------

ns_results <- dunn_results %>%
  filter(p.adj >= 0.05) %>%
  arrange(group1, group2) %>%
  mutate(
    annotation = "ns"
  )

# Place brackets within the displayed log-scale range.
if (nrow(ns_results) > 0) {
  ns_results$y.position <- exp(
    seq(
      log(18),
      log(85),
      length.out = nrow(ns_results)
    )
  )
}


# -----------------------------------------------------------------------------
# 11. Figure 1B: violin plots with embedded boxplots
# -----------------------------------------------------------------------------

p1b <- ggplot(
  positive_tests,
  aes(x = allergen_group, y = sige, fill = allergen_group)
) +
  geom_violin(
    scale = "width",
    trim = TRUE,
    adjust = 0.35,
    width = 0.70,
    alpha = 0.35,
    colour = "black",
    linewidth = 0.35
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    fill = "white",
    colour = "black",
    linewidth = 0.40
  ) +
  scale_fill_manual(values = group_colours) +
  scale_y_log10(
    breaks = sige_breaks,
    labels = sige_breaks
  ) +
  coord_cartesian(
    ylim = c(0.35, 100),
    clip = "off"
  ) +
  labs(
    x = NULL,
    y = "sIgE (kU/L, log scale)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(
      angle = 25,
      hjust = 1,
      face = "bold"
    )
  )

if (nrow(ns_results) > 0) {
  p1b <- p1b +
    stat_pvalue_manual(
      ns_results,
      label = "annotation",
      xmin = "group1",
      xmax = "group2",
      y.position = "y.position",
      tip.length = 0.01,
      bracket.size = 0.4
    )
}

print(p1b)

ggsave(
  filename = file.path("outputs", "figures", "Figure1B_sIgE_violin_boxplot.png"),
  plot = p1b,
  width = 8.5,
  height = 5.5,
  dpi = 300
)


# -----------------------------------------------------------------------------
# 12. Session information for reproducibility
# -----------------------------------------------------------------------------

capture.output(
  sessionInfo(),
  file = file.path("outputs", "sessionInfo.txt")
)
