# VIATAR_QLD

## Overview

This repository contains the R analysis code supporting the Queensland component of the VIATAR project and the associated manuscript:

**“Real-world specific IgE data reveal geospatial variations in aeroallergen sensitization profiles.”**

The study uses real-world serum-specific immunoglobulin E (sIgE) data collected through pathology laboratories to investigate demographic and geospatial variation in aeroallergen sensitisation across Queensland, Australia.

Analyses examine variation in aeroallergen sensitisation according to allergen type, age, sex, urbanicity, coastal proximity, and geographic location, with the aim of improving understanding of regional sensitisation patterns and informing geographically relevant allergy testing, diagnosis, and management strategies.


## Study design

This study is a retrospective cross-sectional analysis of de-identified diagnostic sIgE data from individuals aged ≥12 years in Queensland, Australia, collected between 2019 and 2025.

The Queensland dataset included 744,946 eligible sIgE test records from 160,219 individuals.

Aeroallergen sensitisation was defined as:

**sIgE ≥0.35 kU/L**

Most analyses were conducted at the individual sIgE test-record level, meaning that individuals could contribute multiple test records across different allergens.

Participant-level aggregation was used only for analyses requiring unique individuals, such as the assessment of sensitisation combinations across aeroallergen groups.


## Aeroallergen groups

Aeroallergens were grouped into the following major categories:

- House dust mites
- Grass pollens
- Weed pollens
- Tree pollens
- Moulds
- Animal dander
- Cockroach


## Analysis scripts

Data preparation
R/01_data_preparation.R
Contains data-cleaning and preparation steps used to standardise sIgE values, demographic variables, allergen names, and aeroallergen-group classifications before statistical analysis.
Figure 1 — Sensitisation frequency and sIgE levels
R/02_figure1_sensitisation_frequency_sige_levels.R
Examines overall differences in sensitisation frequency and sIgE concentrations across major aeroallergen groups.
The analysis includes:
- Sensitisation frequency across aeroallergen groups
- Distribution of sIgE concentrations among sensitised tests
- Kruskal–Wallis testing
- Dunn post-hoc comparisons with Benjamini–Hochberg adjustment
- Bar plots and violin plots with embedded boxplots
Figure 2 — Age- and sex-specific patterns
R/03_figure2_age_sex_GAMs.R
Examines non-linear age- and sex-specific patterns in aeroallergen sensitisation frequency and sIgE concentrations using generalised additive models (GAMs).
For sensitisation frequency, separate binomial GAMs with a logit link are fitted for each aeroallergen group to estimate the probability of a positive sIgE test across age in females and males.
For sensitisation level, separate Gaussian GAMs are fitted to log-transformed sIgE concentrations among sensitised tests. Predictions and 95% confidence intervals are back-transformed to the original kU/L scale for presentation.
Figure 3 — Geographic variation in aeroallergen sensitisation
R/04_figure3_geospatial_EB_maps.R
Characterises postcode-level geographic variation in aeroallergen sensitisation across Queensland using Empirical Bayes smoothing.
Panel A estimates EB-smoothed sensitisation frequency for each aeroallergen group.
Panel B estimates EB-smoothed geometric mean sIgE concentrations among positive tests after smoothing postcode-level mean log-transformed sIgE values and back-transforming the estimates to the original scale.
A fixed prior weight of k = 20 is used for the displayed maps.
Table 1 — Sensitisation combinations
R/05_table1_sensitisation_combinations.R
Identifies participant-level patterns of sensitisation across the seven aeroallergen groups.
A participant is considered sensitised to an aeroallergen group if at least one sIgE result within that group is ≥0.35 kU/L.
The script calculates the frequency and percentage of unique sensitisation combinations among participants with at least one positive aeroallergen group.
Table 2 — Coastal proximity and aeroallergen sensitisation
R/06_table2_coastal_proximity.R
Examines variation in aeroallergen sensitisation according to residential coastal proximity:
- <15 km
- 15–30 km
- 30 km
Analyses are conducted at the sIgE test-record level.
The script reports:
- Sensitisation frequencies
- Relative percentage changes
- Risk differences
- Odds ratios
- Overall group comparisons
- Benjamini–Hochberg FDR correction
- Bonferroni-adjusted pairwise comparisons
Supplementary Table 1 — Sex-specific aeroallergen sensitisation
R/07_supp_table1_sex_specific_sensitisation.R
Compares aeroallergen sensitisation frequencies between female and male sIgE test records among individuals aged ≥12 years.
Female–male differences are assessed separately for individual allergens using chi-square or Fisher’s exact tests, with Benjamini–Hochberg correction for multiple comparisons.
Supplementary Table 2 — Age-specific aeroallergen sensitisation
R/08_supp_table2_age_specific_sensitisation.R
Estimates aeroallergen sensitisation frequency across five age groups:
- 12–17 years
- 18–29 years
- 30–44 years
- 45–64 years
- ≥65 years
Differences in sensitisation prevalence across age groups are assessed using chi-square tests, with Benjamini–Hochberg false discovery rate correction across allergens.
Supplementary Table 3 — Urban and regional sensitisation
R/09_supp_table3_urban_regional.R
Compares aeroallergen sensitisation frequencies between urban and regional Queensland at the sIgE test-record level.
The analysis reports:
- Number of tests
- Number sensitised
- Sensitisation percentage
- Risk differences with 95% confidence intervals
- Odds ratios with 95% confidence intervals
- Chi-square or Fisher’s exact tests
- Benjamini–Hochberg FDR-adjusted p-values
Risk differences are expressed as urban minus regional prevalence.
Statistical methods
Depending on the analysis, statistical methods include:
- Chi-square tests
- Fisher’s exact tests
- Two-sample proportion comparisons
- Benjamini–Hochberg false discovery rate correction
- Bonferroni correction for selected pairwise comparisons
- Kruskal–Wallis tests
- Dunn post-hoc comparisons
- Generalised additive models
- Empirical Bayes smoothing
- Risk differences with 95% confidence intervals
- Odds ratios with 95% confidence intervals
Continuous sIgE concentrations are highly right-skewed and are therefore log-transformed where appropriate for modelling and visualisation.
Geographic analyses
Geographic analyses use Queensland residential postcodes.
Spatial analyses include:
- Urban versus regional classification
- Coastal proximity
- Postcode-level geographic variation
- Empirical Bayes smoothing of sensitisation frequency
- Empirical Bayes smoothing of geometric mean sIgE concentration
Australian Statistical Geography Standard postal-area boundaries are used for postcode-level spatial analyses.
Data availability
The individual-level pathology data used in this study are not included in this public repository because of ethics, privacy, confidentiality, and data-governance requirements.
This repository contains analysis code only and does not contain identifiable or de-identified patient-level pathology data.
Users wishing to reproduce the analyses must provide appropriately authorised input data locally.
Local data files should be stored outside the public repository or in directories excluded through .gitignore.
Reproducibility
The analysis scripts use relative local file paths where possible.
Several scripts also generate audit outputs to assist with reproducibility, including:
- Allergen-mapping checks
- Numerical analysis outputs
- Model summaries
- Statistical-test results
- Validation against reported manuscript values
- R session information
These outputs can be used to verify that the analysis reproduces the values reported in the manuscript.
Software
Statistical analyses were primarily conducted using R and RStudio.
Key R packages used across the analyses include:
- dplyr
- tidyr
- stringr
- readr
- readxl
- ggplot2
- mgcv
- sf
- purrr
- RColorBrewer
- patchwork
Authors
Sajjad Chamani
Centre for Immunology and Infection Control

## Repository structure

```text
VIATAR_QLD/
├── README.md
├── .gitignore
└── R/
    ├── 01_data_preparation.R
    ├── 02_figure1_sensitisation_frequency_sige_levels.R
    ├── 03_figure2_age_sex_GAMs.R
    ├── 04_figure3_geospatial_EB_maps.R
    ├── 05_table1_sensitisation_combinations.R
    ├── 06_table2_coastal_proximity.R
    ├── 07_supp_table1_sex_specific_sensitisation.R
    ├── 08_supp_table2_age_specific_sensitisation.R
    └── 09_supp_table3_urban_regional.R


Queensland University of Technology
Brisbane, Queensland, Australia
Citation
Citation information will be added following publication of the associated manuscript.
