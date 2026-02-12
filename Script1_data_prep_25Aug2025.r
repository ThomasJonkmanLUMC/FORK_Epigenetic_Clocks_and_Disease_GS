### Load Required Packages ###
packages <- c("dplyr", "coxme", "stringr")
lapply(packages, function(x) {
  if (!requireNamespace(x, quietly = TRUE)) {
    install.packages(x)
  }
  library(x, character.only = TRUE)
})

## if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
## BiocManager::install("impute")
library(impute)


### load in Epigenetic Clocks ###
### NB: this file includes the super-accurate clock developed by Zhang et al. (NOT Zhang_10). --Thomas
clock_data <- readRDS("file_path/appended_clock_data.rds")

### load in disease, Covariate and Mortality data  ###
disease_data <- read.csv("/Cluster_Filespace/Marioni_Group/Riccardo/Christos_17Oct2024/Data/2024-03-29_diseases.csv")
cov_data <- read.csv("/Cluster_Filespace/Marioni_Group/Riccardo/Christos_17Oct2024/Data/2024-04-19_covariates.csv")
mortality_data <- read.csv("/Cluster_Filespace/Marioni_Group/Riccardo/Christos_17Oct2024/Data/2024-04-19_deaths.csv")

### Load in file linking ids to meth dis ###
GSK_MAP <- readRDS("/Cluster_Filespace/Marioni_Group/Riccardo/Christos_17Oct2024/Data/GS20k_Targets_18869.rds")
names(GSK_MAP)[1] <- "DNAm_ID"

### load in appt records ###
appt_rds <- readRDS("/Cluster_Filespace/Marioni_Group/Riccardo/Christos_17Oct2024/Data/202-04-26_gs_appt.rds")


## Filter Clocks to relevant predictors ##
names(clock_data)[1] <- "DNAm_ID"
filt_clock_data<- select(clock_data, c(DNAm_ID, Horvathv1, Hannum, Lin, PhenoAge, YingCausAge, YingDamAge,
									   YingAdaptAge, Horvathv2, Zhang_10, DunedinPoAm38, DunedinPACE, 
									   DNAmGrimAge, DNAmGrimAge.1, DNAmTL, Zhang_Acc))


### merge appt data to GSK ###
appt_GSK <- merge(appt_rds, GSK_MAP[,c("DNAm_ID","Sample_Name")], by.x= "id", by.y= "Sample_Name")

merged_clock_data_1 <- merge(filt_clock_data,appt_GSK, by = "DNAm_ID")

### merge merged_clock_data with mortality data ### # if columns are same, no need for by.x or by.y just by = id ##
merged_clock_data_2<- merge(merged_clock_data_1, mortality_data[,2:3], by.x= "id", by.y= "id", all.x= TRUE )

#Check NA= alive (not in mortality data)# 
merged_clock_data_2_alive <- merged_clock_data_2 %>% 
  filter(!is.na(dod_ym))

### merge in covariates cov ###
merged_clock_data_3 <- merge(merged_clock_data_2, cov_data[,2:11], by.x= "id", by.y= "id")

### Calculate Time to censor: ###

### Substr APPT dates ###
merged_clock_data_3$GS_yoa <- substr(merged_clock_data_3$ym, 1,4) %>% as.numeric()
merged_clock_data_3$GS_moa <- substr(merged_clock_data_3$ym, 5,6) %>% as.numeric()

### substr YOD dates ### 
merged_clock_data_3$GS_yod <- substr(merged_clock_data_3$dod_ym, 1,4) %>% as.numeric()
merged_clock_data_3$GS_mod <- substr(merged_clock_data_3$dod_ym, 5,6) %>% as.numeric()

# recheck Check NA= as people have not died# 
merged_clock_data_3_alive <- merged_clock_data_3 %>% 
  filter(!is.na(GS_yod))

### new column : 202204- GS APPT ###
merged_clock_data_3$cutoff_minus_GSAPPT <- (2022 - merged_clock_data_3$GS_yoa) + ((04- merged_clock_data_3$GS_moa)/12) %>% as.numeric()

### new column : dod- GSappt ### 
merged_clock_data_3$DOD_minus_GSAPPT <- (merged_clock_data_3$GS_yod - merged_clock_data_3$GS_yoa) + ((merged_clock_data_3$GS_mod - merged_clock_data_3$GS_moa)/12) %>% as.numeric()

## New column T_Censor ##
merged_clock_data_3$t_censor <- ifelse(!is.na(merged_clock_data_3$dod_ym), merged_clock_data_3$DOD_minus_GSAPPT,  merged_clock_data_3$cutoff_minus_GSAPPT)

##  Create new DF trimmed  ###
merged_clock_data_4 <- subset(merged_clock_data_3, select = - c(ym, dod_ym, GS_yoa, GS_moa, GS_yod, GS_mod, cutoff_minus_GSAPPT, DOD_minus_GSAPPT))

### read in and merge WBC proportions ###
### Using EpiDISH as it matches our paper --Thomas ###
WBC_epidish <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/Standardised_Data/EpiDish_WBCs_25Aug2025_REM.rds")
WBC_epidish <- as.data.frame(WBC_epidish)
WBC_epidish <- WBC_epidish[,c("Neu", "Eos", "Baso", "Mono", "Bnv", "Bmem", "CD4Tnv", "CD4Tmem", "CD8Tnv", "CD8Tmem", "Treg", "NK")]
WBC_epidish <- WBC_epidish*100

### Project PC loadings onto WBC proportions --Thomas ###
### Please change the load command to the proper file path. --Thomas ###
load("/path_to_file/PC projection.rda")
WBC_sc <- scale(WBC_epidish, center = pca.trunc$center, scale = pca.trunc$scale)
WBC_PCA <- WBC_sc %*% pca.trunc$rotation
WBC_PCA$DNAm_ID <- row.names(WBC_epidish)

### Sanity check: print correlation between PCs and individual WBC fractions and PCs. --Thomas ###
PC_cor <- cor(WBC_epidish, WBC_PCA)
print(PC_cor)

### Merge WBCs ###
merged_clock_data_5 <- merge(merged_clock_data_4, WBC_PCA, by="DNAm_ID", all.x=T)

### read in the kinship matrix ###
kinship <- readRDS("/Cluster_Filespace/Marioni_Group/Josie/Proteins/ewas/input_files/kinship_matrix_using_fixed_2022-01-28-pedigree.rds")

### regress out kinship from each clock ### 
### Note that WBC props are no longer regressed out here! These will be added as covariates to the cox-model instead --Thomas ###
merged_clock_data_5 <- merged_clock_data_4
for (i in 3:16){
merged_clock_data_5[,i] <- resid(lmekin(merged_clock_data_5[,i] ~ age + (1|merged_clock_data_5$id), 
								varlist = kinship*2, data=merged_clock_data_5, na.action="na.exclude"))
}

### Do not remove WBC PCs --Thomas ###
merged_clock_data_6 <- merged_clock_data_5[,-c(22,25)]



### KNN imputation for missing covariates ###

covs <- as.data.frame(merged_clock_data_6[,c("age","sex","rank","bmi","years","pack_years","units")])
covs$sex <- ifelse(covs$sex=="F", 1, 0)

covs$sex <- scale(covs$sex)
covs$age <- scale(covs$age)
covs$rank <- scale(covs$rank)
covs$bmi <- scale(log(covs$bmi))
covs$years <- scale(covs$years)
covs$pack_years <- scale(log(1 + covs$pack_years))
covs$units <- scale(log(1 + covs$units))

set.seed(1.234)
covs_imp <- as.data.frame(impute.knn(as.matrix(covs), rowmax=0.6)$data)

table(rownames(covs_imp) == rownames(merged_clock_data_6))

### merge imputed covariates with clocks ###
merged_clock_data_7 <- cbind(merged_clock_data_6[,c(1:18,24)], covs_imp[,3:7])

### binary code for sex ###
merged_clock_data_7$sex <- ifelse(merged_clock_data_7$sex=="F", 1, 0)

# Please change file path. --Thomas
saveRDS(merged_clock_data_7, file="/file_path/Merged_clock_data.RDS")
saveRDS(PC_cor, file = "/file_path/WBC_PC_correlations.RDS")
