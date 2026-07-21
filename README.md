<h1>Association between cell composition PCs/epigenetic clocks and health outcomes </h1>
 
 This directory contains some of the R scripts used for Jonkman et al. to test whether cell composition PCs and/or epigenetic clocks are associated with health outcomes. These scripts are a fork of the scripts developed for Mavrommatis et al.
 
 [Jonkman et al. (current scripts)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12132128/) <br>
 [Mavrommatis et al. (original script source)](https://www.nature.com/articles/s41467-025-66106-y)

Only a subset of the original repo's scripts (scripts 1, 1a, 2, and 4) were run here. One additional script (Script 0) was added to the start of the pipeline.

**Script0_zhang_clock.r** <br>
Calculates the [Zhang clock](https://pmc.ncbi.nlm.nih.gov/articles/PMC6708158/) to match Mavrommatis et al.'s methods with our own.

**Script1_data_prep_25Aug2025.r** <br>
Includes data cleaning steps to calculate the time to censoring, impute missing values for the covariates using KNN, and calculates age acceleration residuals that additionally adjust for estimated white blood cell proportions and relatedness (kinship matrix).

**Script1a_disease_meta_data_25Aug2025.r** <br>
Calculates the descriptive statistics for each disease 

**Script2_Cox_disease_15Aug2025.r** <br>
Cox regerssion analyses for each clock/disease pairing. We consider diagnoses made after the blood draw and over the first 10 years of follow up to electronic health records. Death is treated as a censoring event. The code cycles through each disease outcome, runs the regression for each clock (adjusting for age, sex, deprivation, BMI, smoking, education and alcohol) and tests the regression assumptions (Schoenfeld residuals). Basic models adjust for age and sex, fully adjusted models include all covariates. Sex and smoking-stratified and interaction models are also considered.

**Script4_mortality_Cox_AUC_25Aug2025.r** <br>
As per scripts 2 and 3 but, instead of disease outcomes, we consider all-cause mortality.

**Script5_analyse_results.Rmd** <br>
Analyses the results obtained from Scripts 1-4 and makes figures for the manuscript. This script has also been knitted into an [HTML-file](https://thomasjonkmanlumc.github.io/FORK_Epigenetic_Clocks_and_Disease_GS/Knitted%20scripts/Script5_analyse_results.html) to allow closer inspection of which data and scripts generated which figure.

