# Project Title

Stochastic Treatment Regimes in Climate-Health Research: Reassessing Malaria Risk Under Warming Scenarios in Colombia

## Description

Dataset and code shared to reproduce the paper's results "Stochastic Treatment Regimes in Climate-Health Research: Reassessing Malaria Risk Under Warming Scenarios in Colombia". 
The appropriate adjustment can be obtained with the file DAG_analysis.R.
The causal exposure-response curve of the effect of temperature on the incidence of malaria can be obtained with the file exposure_response_curve.py.
The regimenes analysis can be obtained with the file regimenes_temperature.R.
Note that the data on malaria cases were provided by the SIVIGILA and corresponded to anonymized data; for this reason, the data shared in the repository do not contain potentially identifying participant information.

## Dependencies

None

## Data Privacy and Anonymization

This dataset has been processed to ensure complete anonymization and contains no personally identifiable information (PII). All data has been:

- Aggregated at appropriate spatial/temporal scales
- Stripped of any individual identifiers
- Processed to remove direct or indirect identifying elements

The dataset is suitable for public sharing and complies with data privacy standards.

## Privacy Statement

This repository contains datasets that have been carefully processed to protect individual privacy:

### What is NOT included:
- Names, addresses, or contact information
- Individual-level identifiers
- Location data below 25 km resolution
- Timestamps more precise than monthly
- Any data that could be used to re-identify individuals

### Data Processing:
- Spatial aggregation to municipality
- Temporal aggregation to monthly averages

### Compliance:
This dataset meets requirements for public data sharing under applicable privacy regulations.

## Author

Juan David Gutiérrez  

## Versions

causal-curve Python package 1.0.6,
 DAGitty  R package 0.3-4,
 ggdag  R package 0.2.12,
 DAGitty  R package 0.3-4,
 haldensify R package 0.2.3,
 sl3 R package 1.4.3,
 tmle3shift  R package 0.2.1
