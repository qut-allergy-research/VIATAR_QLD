# VIATAR Queensland – Synthetic Demonstration Dataset

## Overview

This directory contains `VIATAR_synthetic_26102026.xlsx`, a fully synthetic demonstration dataset prepared to illustrate the structure of allergen-specific immunoglobulin E (sIgE) test-event data used in the Queensland component of the VIATAR research project.

**All observations and values are fabricated.** No original study records, patient information, or clinical identifiers were copied or released.

## How the Synthetic Data Were Generated

- **Record structure:** The dataset contains 100 artificial allergen-specific test events representing 90 fictional patients. Some fictional patients have multiple records to demonstrate repeated testing and record linkage.
- **Identifiers:** Artificial patient identifiers and 64-character hexadecimal codes were assigned for demonstration purposes, without using original study identifiers.
- **Demographics:** Age, sex, referral information, Queensland postcodes, and collection/test dates were fabricated to provide plausible example records. Dates were assigned within 2019–2025.
- **Laboratory measurements:** Synthetic sIgE and total IgE concentrations were generated using right-skewed distributions. Allergen-group positivity frequencies and selected mean sIgE concentrations were calibrated toward predefined illustrative targets.
- **Quality checks:** The dataset includes 60 female and 40 male test events, with a median event-level age of 43 years (IQR: 23–60 years). These descriptive statistics relate exclusively to the synthetic dataset.

## Sensitisation Definition

For this demonstration dataset, sIgE values >0.35 kU/L were classified as positive and values ≤0.35 kU/L as negative.

This demonstration-specific threshold differs from the ≥0.35 kU/L criterion used in the original VIATAR manuscript.

Sensitisation percentages were calculated using allergen-group-specific test-event denominators rather than unique patient counts.

## Data Availability and Limitations

The synthetic dataset is provided solely for demonstration, educational purposes, and testing of data-processing workflows.

It does not contain actual clinical observations and cannot reproduce the statistical estimates presented in the associated manuscript.

The original individual-level pathology dataset is not publicly available due to ethical, privacy, confidentiality, and data-governance requirements.

## Use of Artificial Intelligence

Generative AI was used to assist with the creation and documentation of this fully synthetic demonstration dataset. All generated values are artificial and do not represent original clinical observations.

## Related Analysis Code

The associated R scripts and main project documentation are available in the [VIATAR_QLD repository](../).
