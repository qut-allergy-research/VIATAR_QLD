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

