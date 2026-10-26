# =============================================================================
# VIATAR_QLD
# Figure 2. Age- and sex-specific patterns in aeroallergen sensitisation
#           frequency and sIgE levels
#
# PURPOSE
# Reproduce the two analyses shown in Figure 2 of the manuscript:
#
#   Panel A — Sensitisation frequency
#   Separate binomial GAMs are fitted for each aeroallergen group to estimate
#   age- and sex-specific probability of a positive sIgE test (>= 0.35 kU/L).
#
#   Panel B — sIgE level among sensitised tests
#   Separate Gaussian GAMs are fitted for each aeroallergen group using
#   log-transformed sIgE concentrations. Predictions and 95% confidence
#   intervals are back-transformed to the original kU/L scale.
#
# ANALYTICAL UNIT
# Each sIgE test record is retained as a separate observation, consistent with
# the current manuscript. No patient-level random effect or clustered standard
# error is used in this Figure 2 script.
#
# DATA PRIVACY
# The real pathology dataset must NOT be committed to the public repository.
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Packages
# -----------------------------------------------------------------------------

library(readxl)
library(dplyr)
library(stringr)
library(readr)
library(tidyr)
library(purrr)
library(mgcv)
library(ggplot2)


# -----------------------------------------------------------------------------
# 2. Output directories
# -----------------------------------------------------------------------------

dir.create("outputs", showWarnings = FALSE)
dir.create(file.path("outputs", "figures"), showWarnings = FALSE)
dir.create(file.path("outputs", "tables"), showWarnings = FALSE)
dir.create(file.path("outputs", "model_diagnostics"), showWarnings = FALSE)


# -----------------------------------------------------------------------------
# 3. Local data path
# -----------------------------------------------------------------------------
# Keep the real dataset outside GitHub.
# If you keep a local copy in data/, make sure "data/" is in .gitignore.

data_path <- file.path("data", "VIATAR_filtered_QLD_01122025.xlsx")

if (!file.exists(data_path)) {
  stop(
    "Data file not found: ", data_path,
    "\nUpdate 'data_path' locally. Do not upload the patient-level dataset to GitHub."
  )
}


# -----------------------------------------------------------------------------
# 4. Read data
# -----------------------------------------------------------------------------

df_raw <- read_excel(data_path)


# -----------------------------------------------------------------------------
# 5. Exact aeroallergen grouping used for Figure 2
# -----------------------------------------------------------------------------
# Names are based on the study grouping used in the manuscript and the labels
# present in the original analysis scripts. Names are normalised to lower case
# and extra spaces are removed before matching.

allergen_lookup <- tibble::tribble(
  ~allergen_key,                         ~allergen_group,

  # House dust mite
  "dermatophagoides farinae",            "Dust mite",
  "dermatophagoides pteronyssinus",      "Dust mite",
  "blomia tropicalis",                   "Dust mite",

  # Mould
  "alternaria alternata",                "Mould",
  "aspergillus fumigatus",               "Mould",
  "cladosporium herbarum",               "Mould",
  "mould mix 1",                         "Mould",

  # Grass pollen
  "bahia grass",                         "Grass pollen",
  "bahia grass (paspalum)",              "Grass pollen",
  "bermuda (couch) grass",               "Grass pollen",
  "johnson grass",                       "Grass pollen",
  "johnson grass (sorghum)",             "Grass pollen",
  "perennial rye grass",                 "Grass pollen",
  "grass mix 2",                         "Grass pollen",

  # Weed pollen
  "common ragweed",                      "Weed pollen",
  "weed mix 1",                          "Weed pollen",

  # Tree pollen
  "eucalyptus",                          "Tree pollen",
  "gum tree",                            "Tree pollen",
  "acacia",                              "Tree pollen",
  "tree mix 7",                          "Tree pollen",

  # Animal dander / animal-derived allergens
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


# -----------------------------------------------------------------------------
# 6. Prepare analysis variables
# -----------------------------------------------------------------------------
#
# IMPORTANT VALUE HANDLING
# Some laboratory results may be stored as strings such as "<0.35".
# Simply removing non-numeric characters would turn "<0.35" into 0.35 and could
# incorrectly classify it as positive. We therefore retain an explicit
# below-detection flag.
#
# For Panel A:
#   "<0.35" is classified as negative.
#
# For Panel B:
#   only sensitised tests are retained, so below-threshold results are excluded.

df <- df_raw %>%
  transmute(
    age = suppressWarnings(as.numeric(`Age.in.Years`)),
    sex = as.character(`Sex.at.birth`),
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
    !is.na(age),
    age >= 12,
    sex %in% c("Female", "Male"),
    !is.na(allergen_group),
    !is.na(sige)
  ) %>%
  mutate(
    sex = factor(sex, levels = c("Female", "Male")),

    # A result reported as "<0.35" is negative even though parse_number()
    # extracts the numeric boundary 0.35.
    sensitised = if_else(
      below_limit,
      0L,
      as.integer(sige >= 0.35)
    )
  )


# -----------------------------------------------------------------------------
# 7. Audit the grouping before modelling
# -----------------------------------------------------------------------------
# These checks are intentionally retained in the GitHub script because they make
# the analysis auditable and help detect unexpected allergen labels.

unmapped_allergens <- df_raw %>%
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
  file.path("outputs", "tables", "Figure2_unmapped_allergen_labels.csv")
)

group_counts <- df %>%
  count(allergen_group, sex, name = "n_test_records")

print(group_counts)

write_csv(
  group_counts,
  file.path("outputs", "tables", "Figure2_test_record_counts.csv")
)


# =============================================================================
# PANEL A — AGE- AND SEX-SPECIFIC SENSITISATION FREQUENCY
# =============================================================================


# -----------------------------------------------------------------------------
# 8. Fit separate binomial GAMs by aeroallergen group
# -----------------------------------------------------------------------------
#
# Model:
#   sensitised ~ sex + s(age, by = sex, k = 10)
#
# The factor-by smooth gives a separate nonlinear age curve for females and
# males. The sex main effect is retained because factor-by smooths are centred.

fit_panel_a <- function(dat) {

  group_name <- unique(dat$allergen_group)

  model <- bam(
    sensitised ~
      sex +
      s(age, by = sex, k = 10),
    data = dat,
    family = binomial(link = "logit"),
    method = "fREML",
    discrete = TRUE
  )

  list(
    group = group_name,
    model = model,
    data = dat
  )
}


panel_a_models <- df %>%
  group_by(allergen_group) %>%
  group_split() %>%
  map(fit_panel_a)

names(panel_a_models) <- map_chr(panel_a_models, "group")


# -----------------------------------------------------------------------------
# 9. Model diagnostics for Panel A
# -----------------------------------------------------------------------------

for (g in names(panel_a_models)) {

  model <- panel_a_models[[g]]$model

  sink(
    file.path(
      "outputs",
      "model_diagnostics",
      paste0("Figure2A_", gsub(" ", "_", g), "_summary.txt")
    )
  )

  cat("Aeroallergen group:", g, "\n\n")
  print(summary(model))
  cat("\n\nGAM CHECK\n")
  gam.check(model)

  sink()
}


# -----------------------------------------------------------------------------
# 10. Generate predictions and 95% CIs for Panel A
# -----------------------------------------------------------------------------

panel_a_predictions <- map_dfr(
  panel_a_models,
  function(x) {

    dat <- x$data

    newdata <- expand_grid(
      age = seq(
        min(dat$age, na.rm = TRUE),
        max(dat$age, na.rm = TRUE),
        length.out = 200
      ),
      sex = levels(dat$sex)
    )

    pred <- predict(
      x$model,
      newdata = newdata,
      type = "link",
      se.fit = TRUE
    )

    newdata %>%
      mutate(
        predicted_probability = plogis(pred$fit),
        lower_95 = plogis(pred$fit - 1.96 * pred$se.fit),
        upper_95 = plogis(pred$fit + 1.96 * pred$se.fit),
        allergen_group = x$group
      )
  }
)


panel_order <- c(
  "Animal dander",
  "Cockroach",
  "Dust mite",
  "Grass pollen",
  "Mould",
  "Tree pollen",
  "Weed pollen"
)

panel_a_predictions <- panel_a_predictions %>%
  mutate(
    allergen_group = factor(allergen_group, levels = panel_order)
  )

write_csv(
  panel_a_predictions,
  file.path("outputs", "tables", "Figure2A_GAM_predictions.csv")
)


# -----------------------------------------------------------------------------
# 11. Plot Figure 2A
# -----------------------------------------------------------------------------

p2a <- ggplot(
  panel_a_predictions,
  aes(
    x = age,
    y = predicted_probability,
    colour = sex,
    fill = sex
  )
) +
  geom_ribbon(
    aes(
      ymin = lower_95,
      ymax = upper_95
    ),
    alpha = 0.20,
    colour = NA
  ) +
  geom_line(linewidth = 1.1) +
  facet_wrap(
    ~ allergen_group,
    ncol = 3,
    scales = "fixed"
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    x = "Age (years)",
    y = "Sensitisation rate",
    colour = "Sex",
    fill = "Sex"
  ) +
  theme_bw(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold"),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold")
  )

print(p2a)

ggsave(
  file.path(
    "outputs",
    "figures",
    "Figure2A_age_sex_sensitisation_frequency_GAM.png"
  ),
  p2a,
  width = 11,
  height = 7,
  dpi = 600
)


# =============================================================================
# PANEL B — AGE- AND SEX-SPECIFIC sIgE LEVELS
# =============================================================================


# -----------------------------------------------------------------------------
# 12. Keep sensitised tests only and log-transform sIgE
# -----------------------------------------------------------------------------

df_positive <- df %>%
  filter(sensitised == 1) %>%
  mutate(
    log_sige = log(sige)
  )


# -----------------------------------------------------------------------------
# 13. Fit separate Gaussian GAMs by aeroallergen group
# -----------------------------------------------------------------------------
#
# Model:
#   log(sIgE) ~ sex + s(age, by = sex, k = 10)
#
# Predictions are made on the log scale and then exponentiated.

fit_panel_b <- function(dat) {

  group_name <- unique(dat$allergen_group)

  model <- bam(
    log_sige ~
      sex +
      s(age, by = sex, k = 10),
    data = dat,
    family = gaussian(),
    method = "fREML",
    discrete = TRUE
  )

  list(
    group = group_name,
    model = model,
    data = dat
  )
}


panel_b_models <- df_positive %>%
  group_by(allergen_group) %>%
  group_split() %>%
  map(fit_panel_b)

names(panel_b_models) <- map_chr(panel_b_models, "group")


# -----------------------------------------------------------------------------
# 14. Model diagnostics for Panel B
# -----------------------------------------------------------------------------

for (g in names(panel_b_models)) {

  model <- panel_b_models[[g]]$model

  sink(
    file.path(
      "outputs",
      "model_diagnostics",
      paste0("Figure2B_", gsub(" ", "_", g), "_summary.txt")
    )
  )

  cat("Aeroallergen group:", g, "\n\n")
  print(summary(model))
  cat("\n\nGAM CHECK\n")
  gam.check(model)

  sink()
}


# -----------------------------------------------------------------------------
# 15. Generate predictions and back-transform to kU/L
# -----------------------------------------------------------------------------

panel_b_predictions <- map_dfr(
  panel_b_models,
  function(x) {

    dat <- x$data

    newdata <- expand_grid(
      age = seq(
        min(dat$age, na.rm = TRUE),
        max(dat$age, na.rm = TRUE),
        length.out = 200
      ),
      sex = levels(dat$sex)
    )

    pred <- predict(
      x$model,
      newdata = newdata,
      type = "link",
      se.fit = TRUE
    )

    newdata %>%
      mutate(
        predicted_sige = exp(pred$fit),
        lower_95 = exp(pred$fit - 1.96 * pred$se.fit),
        upper_95 = exp(pred$fit + 1.96 * pred$se.fit),
        allergen_group = x$group
      )
  }
)

panel_b_predictions <- panel_b_predictions %>%
  mutate(
    allergen_group = factor(allergen_group, levels = panel_order)
  )

write_csv(
  panel_b_predictions,
  file.path("outputs", "tables", "Figure2B_GAM_predictions.csv")
)


# -----------------------------------------------------------------------------
# 16. Plot Figure 2B
# -----------------------------------------------------------------------------

p2b <- ggplot(
  panel_b_predictions,
  aes(
    x = age,
    y = predicted_sige,
    colour = sex,
    fill = sex
  )
) +
  geom_ribbon(
    aes(
      ymin = lower_95,
      ymax = upper_95
    ),
    alpha = 0.20,
    colour = NA
  ) +
  geom_line(linewidth = 1.1) +
  facet_wrap(
    ~ allergen_group,
    ncol = 3
  ) +
  scale_y_log10() +
  labs(
    x = "Age (years)",
    y = "Predicted sIgE (kU/L)",
    colour = "Sex",
    fill = "Sex"
  ) +
  theme_bw(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold"),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold")
  )

print(p2b)

ggsave(
  file.path(
    "outputs",
    "figures",
    "Figure2B_age_sex_sIgE_levels_GAM.png"
  ),
  p2b,
  width = 11,
  height = 7,
  dpi = 600
)


# -----------------------------------------------------------------------------
# 17. Save session information for reproducibility
# -----------------------------------------------------------------------------

capture.output(
  sessionInfo(),
  file = file.path("outputs", "Figure2_sessionInfo.txt")
)
