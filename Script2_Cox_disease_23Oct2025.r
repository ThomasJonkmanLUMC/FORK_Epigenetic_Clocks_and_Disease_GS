# Load Required Packages
packages <- c("dplyr", "survival")
lapply(packages, function(x) {
  if (!requireNamespace(x, quietly = TRUE)) {
    install.packages(x)
  }
  library(x, character.only = TRUE)
})

# Read in clock data
# NB: this data also contains the WBC PCs --Thomas
# Please change file path. --Thomas
clock_data <- readRDS("/file_path/Merged_clock_data.RDS")

# merge in smoking status for downstream subsetting
cov_data <- read.csv("/Cluster_Filespace/Marioni_Group/Riccardo/Christos_17Oct2024/Data/2024-04-19_covariates.csv")

# create binary smoking (NB 633 missing pack year info)
cov_data$smoke <- ifelse(cov_data$pack_years == 0, 0, 1)

clock_data <- merge(clock_data,cov_data[,c("id", "smoke")], by="id")

# List of clocks to loop through
# NB: includes the super-accurate clock developed by Zhang et al.
clocks <- c("Hannum", "Horvathv1", "Zhang_Acc", "PhenoAge", "DNAmGrimAge.1", "DunedinPACE")

# Read in diseases
diseases <- read.csv("/Cluster_Filespace/Marioni_Group/Riccardo/Christos_17Oct2024/Data/2024-03-29_diseases.csv", stringsAsFactors = FALSE)
diseases <- diseases[-which(diseases$Source == "Primary_Care" & diseases$GP_Consent==0),]

# Get unique diseases
dis_list <- unique(diseases$Disease)

# Check number of unique diseases
cat("Number of unique diseases:", length(dis_list), "\n")

# Create storage lists for results
# The following models will be run: clocks alone, clocks adjusted for WBC PCs (11 PCs, including one more each time), and WBC PCs alone (all 11). --Thomas


#Result objects for Cox-regression models.
res.cox.clocks <- setNames(vector("list", length(clocks)), clocks)

res.cox.clocks.pc1 <- setNames(vector("list", length(clocks)), clocks)
res.cox.clocks.pc1.2 <- setNames(vector("list", length(clocks)), clocks)
res.cox.clocks.pc1.3 <- setNames(vector("list", length(clocks)), clocks)
res.cox.clocks.pc1.4 <- setNames(vector("list", length(clocks)), clocks)
res.cox.clocks.pc1.5 <- setNames(vector("list", length(clocks)), clocks)
res.cox.clocks.pc1.6 <- setNames(vector("list", length(clocks)), clocks)
res.cox.clocks.pc1.7 <- setNames(vector("list", length(clocks)), clocks)
res.cox.clocks.pc1.8 <- setNames(vector("list", length(clocks)), clocks)
res.cox.clocks.pc1.9 <- setNames(vector("list", length(clocks)), clocks)
res.cox.clocks.pc1.10 <- setNames(vector("list", length(clocks)), clocks)
res.cox.clocks.pc1.11 <- setNames(vector("list", length(clocks)), clocks)

res.cox.pc1 <- setNames(vector("list", length(clocks)), clocks)
res.cox.pc2 <- setNames(vector("list", length(clocks)), clocks)
res.cox.pc3 <- setNames(vector("list", length(clocks)), clocks)
res.cox.pc4 <- setNames(vector("list", length(clocks)), clocks)
res.cox.pc5 <- setNames(vector("list", length(clocks)), clocks)
res.cox.pc6 <- setNames(vector("list", length(clocks)), clocks)
res.cox.pc7 <- setNames(vector("list", length(clocks)), clocks)
res.cox.pc8 <- setNames(vector("list", length(clocks)), clocks)
res.cox.pc9 <- setNames(vector("list", length(clocks)), clocks)
res.cox.pc10 <- setNames(vector("list", length(clocks)), clocks)
res.cox.pc11 <- setNames(vector("list", length(clocks)), clocks)


#Result objects for proportional hazards assumption tests.
res.cox.clocks.zph.local <- setNames(vector("list", length(clocks)), clocks)
res.cox.clocks.zph.global <- setNames(vector("list", length(clocks)), clocks)

res.cox.clocks.pc1.zph.local <- setNames(vector("list", length(clocks)), clocks)
res.cox.clocks.pc1.zph.global <- setNames(vector("list", length(clocks)), clocks)
res.cox.clocks.pc1.2.zph.local <- setNames(vector("list", length(clocks)), clocks)
res.cox.clocks.pc1.2.zph.global <- setNames(vector("list", length(clocks)), clocks)
res.cox.clocks.pc1.3.zph.local <- setNames(vector("list", length(clocks)), clocks)
res.cox.clocks.pc1.3.zph.global <- setNames(vector("list", length(clocks)), clocks)
res.cox.clocks.pc1.4.zph.local <- setNames(vector("list", length(clocks)), clocks)
res.cox.clocks.pc1.4.zph.global <- setNames(vector("list", length(clocks)), clocks)
res.cox.clocks.pc1.5.zph.local <- setNames(vector("list", length(clocks)), clocks)
res.cox.clocks.pc1.5.zph.global <- setNames(vector("list", length(clocks)), clocks)
res.cox.clocks.pc1.6.zph.local <- setNames(vector("list", length(clocks)), clocks)
res.cox.clocks.pc1.6.zph.global <- setNames(vector("list", length(clocks)), clocks)
res.cox.clocks.pc1.7.zph.local <- setNames(vector("list", length(clocks)), clocks)
res.cox.clocks.pc1.7.zph.global <- setNames(vector("list", length(clocks)), clocks)
res.cox.clocks.pc1.8.zph.local <- setNames(vector("list", length(clocks)), clocks)
res.cox.clocks.pc1.8.zph.global <- setNames(vector("list", length(clocks)), clocks)
res.cox.clocks.pc1.9.zph.local <- setNames(vector("list", length(clocks)), clocks)
res.cox.clocks.pc1.9.zph.global <- setNames(vector("list", length(clocks)), clocks)
res.cox.clocks.pc1.10.zph.local <- setNames(vector("list", length(clocks)), clocks)
res.cox.clocks.pc1.10.zph.global <- setNames(vector("list", length(clocks)), clocks)
res.cox.clocks.pc1.11.zph.local <- setNames(vector("list", length(clocks)), clocks)
res.cox.clocks.pc1.11.zph.global <- setNames(vector("list", length(clocks)), clocks)

res.cox.pc1.zph.local <- setNames(vector("list", length(clocks)), clocks)
res.cox.pc1.zph.global <- setNames(vector("list", length(clocks)), clocks)
res.cox.pc2.zph.local <- setNames(vector("list", length(clocks)), clocks)
res.cox.pc2.zph.global <- setNames(vector("list", length(clocks)), clocks)
res.cox.pc3.zph.local <- setNames(vector("list", length(clocks)), clocks)
res.cox.pc3.zph.global <- setNames(vector("list", length(clocks)), clocks)
res.cox.pc4.zph.local <- setNames(vector("list", length(clocks)), clocks)
res.cox.pc4.zph.global <- setNames(vector("list", length(clocks)), clocks)
res.cox.pc5.zph.local <- setNames(vector("list", length(clocks)), clocks)
res.cox.pc5.zph.global <- setNames(vector("list", length(clocks)), clocks)
res.cox.pc6.zph.local <- setNames(vector("list", length(clocks)), clocks)
res.cox.pc6.zph.global <- setNames(vector("list", length(clocks)), clocks)
res.cox.pc7.zph.local <- setNames(vector("list", length(clocks)), clocks)
res.cox.pc7.zph.global <- setNames(vector("list", length(clocks)), clocks)
res.cox.pc8.zph.local <- setNames(vector("list", length(clocks)), clocks)
res.cox.pc8.zph.global <- setNames(vector("list", length(clocks)), clocks)
res.cox.pc9.zph.local <- setNames(vector("list", length(clocks)), clocks)
res.cox.pc9.zph.global <- setNames(vector("list", length(clocks)), clocks)
res.cox.pc10.zph.local <- setNames(vector("list", length(clocks)), clocks)
res.cox.pc10.zph.global <- setNames(vector("list", length(clocks)), clocks)
res.cox.pc11.zph.local <- setNames(vector("list", length(clocks)), clocks)
res.cox.pc11.zph.global <- setNames(vector("list", length(clocks)), clocks)


#Result objects for LRT of WBC effect for clock models.
res.an.clocks <- setNames(vector("list", length(clocks)), clocks)


# Disease summary table
disease_summary <- list()

# Loop through each disease
# Currently looping through only 3 diseases as a test. --Thomas
for (disease in dis_list[1:3]) {
  disease_data <- diseases %>% filter(Disease == disease)
  cat(sprintf("Processing disease %s\n", disease))
  
  # Calculate time to incident
  disease_data <- disease_data %>%
    mutate(yod = as.numeric(substr(dt1_ym, 1, 4)),
           mod = as.numeric(substr(dt1_ym, 5, 6)),
           yoa = as.numeric(substr(gs_appt, 1, 4)),
           moa = as.numeric(substr(gs_appt, 5, 6)),
           t_disease = (yod - yoa) + ((mod - moa) / 12)) %>%
    select(-mod, -yoa, -moa)  # Remove excess columns

  # Merge with clock data
  merged_data <- merge(disease_data, clock_data, by = "id", all.y = TRUE)
  
  # Filter for incident cases
  merged_data_2 <- merged_data %>% filter(incident != 0 | is.na(incident))
  merged_data_2 <- merged_data_2 %>% filter(t_disease > 0 | is.na(t_disease))
  merged_data_2 <- merged_data_2 %>%
    mutate(DP = ifelse(t_disease > 0, 1, 0)) %>%
    filter(is.na(yod) | yod > 1980) 
 
  # Order data and remove duplicates
  merged_data_ord <- merged_data_2[order(merged_data_2$dt1_ym), ]
  merged_data_ord <- merged_data_ord[!duplicated(na.omit(merged_data_ord$id)), ]
  
  # Create event column and time to event
  merged_data_ord <- merged_data_ord %>%
    mutate(event = ifelse(is.na(incident), 0, 1),
           t_event = ifelse(event == 1, t_disease, t_censor))

  # Trim data
  final_data <- merged_data_ord %>% select(-yod, -dt1_ym, -gs_appt)

  # Skip disease if no events
  if (sum(final_data$event == 1) == 0) next

  # Gender filtering
  prop_F <-  length(which(final_data$sex == 1 & final_data$event ==1))/ length(which(final_data$event ==1 ))
  if (prop_F > 0.9) {
    final_data <- final_data %>% filter(sex == 1)
  } else if (prop_F < 0.1) {
    final_data <- final_data %>% filter(sex == 0)
  }


  # Filter to 10 year analyses
  final_data$event[final_data$t_event > 10 & final_data$event == 1] <- 0
  final_data$t_event[final_data$t_event > 10] <- 10 

  # Count events
  n_event <- sum(final_data$event == 1)
  final_data$n_event <- n_event

  # Skip if n_event < 30
  if (n_event < 30) next

  ## Run Cox regression model
  if (prop_F > 0.1 & prop_F < 0.9) {
    # Case 1: prop_F between 0.1 and 0.9
    for (clock in clocks) {

      model_formula_clock <- as.formula(paste0("Surv(t_event, event) ~ scale(", clock, ") + age + bmi + years + pack_years + units + rank + sex"))
      
      model_formula_clock.pc1 <- as.formula(paste0("Surv(t_event, event) ~ scale(", clock, ") + age + bmi + years + pack_years + units + rank + sex + PC1"))
      model_formula_clock.pc1.2 <- as.formula(paste0("Surv(t_event, event) ~ scale(", clock, ") + age + bmi + years + pack_years + units + rank + sex + PC1 + PC2"))
      model_formula_clock.pc1.3 <- as.formula(paste0("Surv(t_event, event) ~ scale(", clock, ") + age + bmi + years + pack_years + units + rank + sex + PC1 + PC2 + PC3"))
      model_formula_clock.pc1.4 <- as.formula(paste0("Surv(t_event, event) ~ scale(", clock, ") + age + bmi + years + pack_years + units + rank + sex + PC1 + PC2 + PC3 + PC4"))
      model_formula_clock.pc1.5 <- as.formula(paste0("Surv(t_event, event) ~ scale(", clock, ") + age + bmi + years + pack_years + units + rank + sex + PC1 + PC2 + PC3 + PC4 + PC5"))
      model_formula_clock.pc1.6 <- as.formula(paste0("Surv(t_event, event) ~ scale(", clock, ") + age + bmi + years + pack_years + units + rank + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6"))
      model_formula_clock.pc1.7 <- as.formula(paste0("Surv(t_event, event) ~ scale(", clock, ") + age + bmi + years + pack_years + units + rank + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7"))
      model_formula_clock.pc1.8 <- as.formula(paste0("Surv(t_event, event) ~ scale(", clock, ") + age + bmi + years + pack_years + units + rank + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8"))
      model_formula_clock.pc1.9 <- as.formula(paste0("Surv(t_event, event) ~ scale(", clock, ") + age + bmi + years + pack_years + units + rank + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9"))
      model_formula_clock.pc1.10 <- as.formula(paste0("Surv(t_event, event) ~ scale(", clock, ") + age + bmi + years + pack_years + units + rank + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))
      model_formula_clock.pc1.11 <- as.formula(paste0("Surv(t_event, event) ~ scale(", clock, ") + age + bmi + years + pack_years + units + rank + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11"))
      
      model_formula_pc1 <- as.formula(paste0("Surv(t_event, event) ~ PC1 + age + bmi + years + pack_years + units + rank + sex"))
      model_formula_pc2 <- as.formula(paste0("Surv(t_event, event) ~ PC2 + age + bmi + years + pack_years + units + rank + sex"))
      model_formula_pc3 <- as.formula(paste0("Surv(t_event, event) ~ PC3 + age + bmi + years + pack_years + units + rank + sex"))
      model_formula_pc4 <- as.formula(paste0("Surv(t_event, event) ~ PC4 + age + bmi + years + pack_years + units + rank + sex"))
      model_formula_pc5 <- as.formula(paste0("Surv(t_event, event) ~ PC5 + age + bmi + years + pack_years + units + rank + sex"))
      model_formula_pc6 <- as.formula(paste0("Surv(t_event, event) ~ PC6 + age + bmi + years + pack_years + units + rank + sex"))
      model_formula_pc7 <- as.formula(paste0("Surv(t_event, event) ~ PC7 + age + bmi + years + pack_years + units + rank + sex"))
      model_formula_pc8 <- as.formula(paste0("Surv(t_event, event) ~ PC8 + age + bmi + years + pack_years + units + rank + sex"))
      model_formula_pc9 <- as.formula(paste0("Surv(t_event, event) ~ PC9 + age + bmi + years + pack_years + units + rank + sex"))
      model_formula_pc10 <- as.formula(paste0("Surv(t_event, event) ~ PC10 + age + bmi + years + pack_years + units + rank + sex"))
      model_formula_pc11 <- as.formula(paste0("Surv(t_event, event) ~ PC11 + age + bmi + years + pack_years + units + rank + sex"))
      
      
      
      res.cox.clocks[[clock]][[disease]] <- coxph(model_formula_clock, data = final_data)
      
      res.cox.clocks.pc1[[clock]][[disease]] <- coxph(model_formula_clock.pc1, data = final_data)
      res.cox.clocks.pc1.2[[clock]][[disease]] <- coxph(model_formula_clock.pc1.2, data = final_data)
      res.cox.clocks.pc1.3[[clock]][[disease]] <- coxph(model_formula_clock.pc1.3, data = final_data)
      res.cox.clocks.pc1.4[[clock]][[disease]] <- coxph(model_formula_clock.pc1.4, data = final_data)
      res.cox.clocks.pc1.5[[clock]][[disease]] <- coxph(model_formula_clock.pc1.5, data = final_data)
      res.cox.clocks.pc1.6[[clock]][[disease]] <- coxph(model_formula_clock.pc1.6, data = final_data)
      res.cox.clocks.pc1.7[[clock]][[disease]] <- coxph(model_formula_clock.pc1.7, data = final_data)
      res.cox.clocks.pc1.8[[clock]][[disease]] <- coxph(model_formula_clock.pc1.8, data = final_data)
      res.cox.clocks.pc1.9[[clock]][[disease]] <- coxph(model_formula_clock.pc1.9, data = final_data)
      res.cox.clocks.pc1.10[[clock]][[disease]] <- coxph(model_formula_clock.pc1.10, data = final_data)
      res.cox.clocks.pc1.11[[clock]][[disease]] <- coxph(model_formula_clock.pc1.11, data = final_data)
      
      res.cox.pc1[[clock]][[disease]] <- coxph(model_formula_pc1, data = final_data)
      res.cox.pc2[[clock]][[disease]] <- coxph(model_formula_pc2, data = final_data)
      res.cox.pc3[[clock]][[disease]] <- coxph(model_formula_pc3, data = final_data)
      res.cox.pc4[[clock]][[disease]] <- coxph(model_formula_pc4, data = final_data)
      res.cox.pc5[[clock]][[disease]] <- coxph(model_formula_pc5, data = final_data)
      res.cox.pc6[[clock]][[disease]] <- coxph(model_formula_pc6, data = final_data)
      res.cox.pc7[[clock]][[disease]] <- coxph(model_formula_pc7, data = final_data)
      res.cox.pc8[[clock]][[disease]] <- coxph(model_formula_pc8, data = final_data)
      res.cox.pc9[[clock]][[disease]] <- coxph(model_formula_pc9, data = final_data)
      res.cox.pc10[[clock]][[disease]] <- coxph(model_formula_pc10, data = final_data)
      res.cox.pc11[[clock]][[disease]] <- coxph(model_formula_pc11, data = final_data)
      
      #Calculate the level of evidence that WBC composition affects the clock-disease association. --Thomas
      fit1 <- coxph(model_formula_clock, data = final_data)
      fit2 <- coxph(model_formula_clock.pc1.11, data = final_data)
      res.an.clocks[[clock]][[disease]] <- anova(fit1, fit2, test = "LRT")$`Pr(>|Chi|)`[2]
      
    }
  } else if (prop_F <= 0.1) {
    # Case 2: prop_F less than 0.1
    for (clock in clocks) {

      model_formula_clock <- as.formula(paste0("Surv(t_event, event) ~ scale(", clock, ") + age + bmi + years + pack_years + units + rank"))
      
      model_formula_clock.pc1 <- as.formula(paste0("Surv(t_event, event) ~ scale(", clock, ") + age + bmi + years + pack_years + units + rank + PC1"))
      model_formula_clock.pc1.2 <- as.formula(paste0("Surv(t_event, event) ~ scale(", clock, ") + age + bmi + years + pack_years + units + rank + PC1 + PC2"))
      model_formula_clock.pc1.3 <- as.formula(paste0("Surv(t_event, event) ~ scale(", clock, ") + age + bmi + years + pack_years + units + rank + PC1 + PC2 + PC3"))
      model_formula_clock.pc1.4 <- as.formula(paste0("Surv(t_event, event) ~ scale(", clock, ") + age + bmi + years + pack_years + units + rank + PC1 + PC2 + PC3 + PC4"))
      model_formula_clock.pc1.5 <- as.formula(paste0("Surv(t_event, event) ~ scale(", clock, ") + age + bmi + years + pack_years + units + rank + PC1 + PC2 + PC3 + PC4 + PC5"))
      model_formula_clock.pc1.6 <- as.formula(paste0("Surv(t_event, event) ~ scale(", clock, ") + age + bmi + years + pack_years + units + rank + PC1 + PC2 + PC3 + PC4 + PC5 + PC6"))
      model_formula_clock.pc1.7 <- as.formula(paste0("Surv(t_event, event) ~ scale(", clock, ") + age + bmi + years + pack_years + units + rank + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7"))
      model_formula_clock.pc1.8 <- as.formula(paste0("Surv(t_event, event) ~ scale(", clock, ") + age + bmi + years + pack_years + units + rank + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8"))
      model_formula_clock.pc1.9 <- as.formula(paste0("Surv(t_event, event) ~ scale(", clock, ") + age + bmi + years + pack_years + units + rank + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9"))
      model_formula_clock.pc1.10 <- as.formula(paste0("Surv(t_event, event) ~ scale(", clock, ") + age + bmi + years + pack_years + units + rank + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))
      model_formula_clock.pc1.11 <- as.formula(paste0("Surv(t_event, event) ~ scale(", clock, ") + age + bmi + years + pack_years + units + rank + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11"))
      
      model_formula_pc1 <- as.formula(paste0("Surv(t_event, event) ~ PC1 + age + bmi + years + pack_years + units + rank"))
      model_formula_pc2 <- as.formula(paste0("Surv(t_event, event) ~ PC2 + age + bmi + years + pack_years + units + rank"))
      model_formula_pc3 <- as.formula(paste0("Surv(t_event, event) ~ PC3 + age + bmi + years + pack_years + units + rank"))
      model_formula_pc4 <- as.formula(paste0("Surv(t_event, event) ~ PC4 + age + bmi + years + pack_years + units + rank"))
      model_formula_pc5 <- as.formula(paste0("Surv(t_event, event) ~ PC5 + age + bmi + years + pack_years + units + rank"))
      model_formula_pc6 <- as.formula(paste0("Surv(t_event, event) ~ PC6 + age + bmi + years + pack_years + units + rank"))
      model_formula_pc7 <- as.formula(paste0("Surv(t_event, event) ~ PC7 + age + bmi + years + pack_years + units + rank"))
      model_formula_pc8 <- as.formula(paste0("Surv(t_event, event) ~ PC8 + age + bmi + years + pack_years + units + rank"))
      model_formula_pc9 <- as.formula(paste0("Surv(t_event, event) ~ PC9 + age + bmi + years + pack_years + units + rank"))
      model_formula_pc10 <- as.formula(paste0("Surv(t_event, event) ~ PC10 + age + bmi + years + pack_years + units + rank"))
      model_formula_pc11 <- as.formula(paste0("Surv(t_event, event) ~ PC11 + age + bmi + years + pack_years + units + rank"))

      
      
      res.cox.clocks[[clock]][[disease]] <- coxph(model_formula_clock, data = final_data[final_data$sex==0,])
      
      res.cox.clocks.pc1[[clock]][[disease]] <- coxph(model_formula_clock.pc1, data = final_data[final_data$sex==0,])
      res.cox.clocks.pc1.2[[clock]][[disease]] <- coxph(model_formula_clock.pc1.2, data = final_data[final_data$sex==0,])
      res.cox.clocks.pc1.3[[clock]][[disease]] <- coxph(model_formula_clock.pc1.3, data = final_data[final_data$sex==0,])
      res.cox.clocks.pc1.4[[clock]][[disease]] <- coxph(model_formula_clock.pc1.4, data = final_data[final_data$sex==0,])
      res.cox.clocks.pc1.5[[clock]][[disease]] <- coxph(model_formula_clock.pc1.5, data = final_data[final_data$sex==0,])
      res.cox.clocks.pc1.6[[clock]][[disease]] <- coxph(model_formula_clock.pc1.6, data = final_data[final_data$sex==0,])
      res.cox.clocks.pc1.7[[clock]][[disease]] <- coxph(model_formula_clock.pc1.7, data = final_data[final_data$sex==0,])
      res.cox.clocks.pc1.8[[clock]][[disease]] <- coxph(model_formula_clock.pc1.8, data = final_data[final_data$sex==0,])
      res.cox.clocks.pc1.9[[clock]][[disease]] <- coxph(model_formula_clock.pc1.9, data = final_data[final_data$sex==0,])
      res.cox.clocks.pc1.10[[clock]][[disease]] <- coxph(model_formula_clock.pc1.10, data = final_data[final_data$sex==0,])
      res.cox.clocks.pc1.11[[clock]][[disease]] <- coxph(model_formula_clock.pc1.11, data = final_data[final_data$sex==0,])
      
      res.cox.pc1[[clock]][[disease]] <- coxph(model_formula_pc1, data = final_data[final_data$sex==0,])
      res.cox.pc2[[clock]][[disease]] <- coxph(model_formula_pc2, data = final_data[final_data$sex==0,])
      res.cox.pc3[[clock]][[disease]] <- coxph(model_formula_pc3, data = final_data[final_data$sex==0,])
      res.cox.pc4[[clock]][[disease]] <- coxph(model_formula_pc4, data = final_data[final_data$sex==0,])
      res.cox.pc5[[clock]][[disease]] <- coxph(model_formula_pc5, data = final_data[final_data$sex==0,])
      res.cox.pc6[[clock]][[disease]] <- coxph(model_formula_pc6, data = final_data[final_data$sex==0,])
      res.cox.pc7[[clock]][[disease]] <- coxph(model_formula_pc7, data = final_data[final_data$sex==0,])
      res.cox.pc8[[clock]][[disease]] <- coxph(model_formula_pc8, data = final_data[final_data$sex==0,])
      res.cox.pc9[[clock]][[disease]] <- coxph(model_formula_pc9, data = final_data[final_data$sex==0,])
      res.cox.pc10[[clock]][[disease]] <- coxph(model_formula_pc10, data = final_data[final_data$sex==0,])
      res.cox.pc11[[clock]][[disease]] <- coxph(model_formula_pc11, data = final_data[final_data$sex==0,])

    }
  } else {
    # Case 3: prop_F greater than or equal to 0.9
    for (clock in clocks) {
      
      model_formula_clock <- as.formula(paste0("Surv(t_event, event) ~ scale(", clock, ") + age + bmi + years + pack_years + units + rank"))
      
      model_formula_clock.pc1 <- as.formula(paste0("Surv(t_event, event) ~ scale(", clock, ") + age + bmi + years + pack_years + units + rank + PC1"))
      model_formula_clock.pc1.2 <- as.formula(paste0("Surv(t_event, event) ~ scale(", clock, ") + age + bmi + years + pack_years + units + rank + PC1 + PC2"))
      model_formula_clock.pc1.3 <- as.formula(paste0("Surv(t_event, event) ~ scale(", clock, ") + age + bmi + years + pack_years + units + rank + PC1 + PC2 + PC3"))
      model_formula_clock.pc1.4 <- as.formula(paste0("Surv(t_event, event) ~ scale(", clock, ") + age + bmi + years + pack_years + units + rank + PC1 + PC2 + PC3 + PC4"))
      model_formula_clock.pc1.5 <- as.formula(paste0("Surv(t_event, event) ~ scale(", clock, ") + age + bmi + years + pack_years + units + rank + PC1 + PC2 + PC3 + PC4 + PC5"))
      model_formula_clock.pc1.6 <- as.formula(paste0("Surv(t_event, event) ~ scale(", clock, ") + age + bmi + years + pack_years + units + rank + PC1 + PC2 + PC3 + PC4 + PC5 + PC6"))
      model_formula_clock.pc1.7 <- as.formula(paste0("Surv(t_event, event) ~ scale(", clock, ") + age + bmi + years + pack_years + units + rank + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7"))
      model_formula_clock.pc1.8 <- as.formula(paste0("Surv(t_event, event) ~ scale(", clock, ") + age + bmi + years + pack_years + units + rank + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8"))
      model_formula_clock.pc1.9 <- as.formula(paste0("Surv(t_event, event) ~ scale(", clock, ") + age + bmi + years + pack_years + units + rank + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9"))
      model_formula_clock.pc1.10 <- as.formula(paste0("Surv(t_event, event) ~ scale(", clock, ") + age + bmi + years + pack_years + units + rank + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))
      model_formula_clock.pc1.11 <- as.formula(paste0("Surv(t_event, event) ~ scale(", clock, ") + age + bmi + years + pack_years + units + rank + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11"))
      
      model_formula_pc1 <- as.formula(paste0("Surv(t_event, event) ~ PC1 + age + bmi + years + pack_years + units + rank"))
      model_formula_pc2 <- as.formula(paste0("Surv(t_event, event) ~ PC2 + age + bmi + years + pack_years + units + rank"))
      model_formula_pc3 <- as.formula(paste0("Surv(t_event, event) ~ PC3 + age + bmi + years + pack_years + units + rank"))
      model_formula_pc4 <- as.formula(paste0("Surv(t_event, event) ~ PC4 + age + bmi + years + pack_years + units + rank"))
      model_formula_pc5 <- as.formula(paste0("Surv(t_event, event) ~ PC5 + age + bmi + years + pack_years + units + rank"))
      model_formula_pc6 <- as.formula(paste0("Surv(t_event, event) ~ PC6 + age + bmi + years + pack_years + units + rank"))
      model_formula_pc7 <- as.formula(paste0("Surv(t_event, event) ~ PC7 + age + bmi + years + pack_years + units + rank"))
      model_formula_pc8 <- as.formula(paste0("Surv(t_event, event) ~ PC8 + age + bmi + years + pack_years + units + rank"))
      model_formula_pc9 <- as.formula(paste0("Surv(t_event, event) ~ PC9 + age + bmi + years + pack_years + units + rank"))
      model_formula_pc10 <- as.formula(paste0("Surv(t_event, event) ~ PC10 + age + bmi + years + pack_years + units + rank"))
      model_formula_pc11 <- as.formula(paste0("Surv(t_event, event) ~ PC11 + age + bmi + years + pack_years + units + rank"))
      
      
      res.cox.clocks[[clock]][[disease]] <- coxph(model_formula_clock, data = final_data[final_data$sex==1,])
      
      res.cox.clocks.pc1[[clock]][[disease]] <- coxph(model_formula_clock.pc1, data = final_data[final_data$sex==1,])
      res.cox.clocks.pc1.2[[clock]][[disease]] <- coxph(model_formula_clock.pc1.2, data = final_data[final_data$sex==1,])
      res.cox.clocks.pc1.3[[clock]][[disease]] <- coxph(model_formula_clock.pc1.3, data = final_data[final_data$sex==1,])
      res.cox.clocks.pc1.4[[clock]][[disease]] <- coxph(model_formula_clock.pc1.4, data = final_data[final_data$sex==1,])
      res.cox.clocks.pc1.5[[clock]][[disease]] <- coxph(model_formula_clock.pc1.5, data = final_data[final_data$sex==1,])
      res.cox.clocks.pc1.6[[clock]][[disease]] <- coxph(model_formula_clock.pc1.6, data = final_data[final_data$sex==1,])
      res.cox.clocks.pc1.7[[clock]][[disease]] <- coxph(model_formula_clock.pc1.7, data = final_data[final_data$sex==1,])
      res.cox.clocks.pc1.8[[clock]][[disease]] <- coxph(model_formula_clock.pc1.8, data = final_data[final_data$sex==1,])
      res.cox.clocks.pc1.9[[clock]][[disease]] <- coxph(model_formula_clock.pc1.9, data = final_data[final_data$sex==1,])
      res.cox.clocks.pc1.10[[clock]][[disease]] <- coxph(model_formula_clock.pc1.10, data = final_data[final_data$sex==1,])
      res.cox.clocks.pc1.11[[clock]][[disease]] <- coxph(model_formula_clock.pc1.11, data = final_data[final_data$sex==1,])
      
      res.cox.pc1[[clock]][[disease]] <- coxph(model_formula_pc1, data = final_data[final_data$sex==1,])
      res.cox.pc2[[clock]][[disease]] <- coxph(model_formula_pc2, data = final_data[final_data$sex==1,])
      res.cox.pc3[[clock]][[disease]] <- coxph(model_formula_pc3, data = final_data[final_data$sex==1,])
      res.cox.pc4[[clock]][[disease]] <- coxph(model_formula_pc4, data = final_data[final_data$sex==1,])
      res.cox.pc5[[clock]][[disease]] <- coxph(model_formula_pc5, data = final_data[final_data$sex==1,])
      res.cox.pc6[[clock]][[disease]] <- coxph(model_formula_pc6, data = final_data[final_data$sex==1,])
      res.cox.pc7[[clock]][[disease]] <- coxph(model_formula_pc7, data = final_data[final_data$sex==1,])
      res.cox.pc8[[clock]][[disease]] <- coxph(model_formula_pc8, data = final_data[final_data$sex==1,])
      res.cox.pc9[[clock]][[disease]] <- coxph(model_formula_pc9, data = final_data[final_data$sex==1,])
      res.cox.pc10[[clock]][[disease]] <- coxph(model_formula_pc10, data = final_data[final_data$sex==1,])
      res.cox.pc11[[clock]][[disease]] <- coxph(model_formula_pc11, data = final_data[final_data$sex==1,])

    }
  }

cat(sprintf("Running zph for %s\n", disease))

for (clock in clocks) {

  res.cox.clocks.zph.local[[clock]][[disease]]         <- cox.zph(res.cox.clocks[[clock]][[disease]])$table[1, ]
  res.cox.clocks.zph.global[[clock]][[disease]]        <- cox.zph(res.cox.clocks[[clock]][[disease]])$table["GLOBAL", ]
  
  res.cox.clocks.pc1.zph.local[[clock]][[disease]]     <- cox.zph(res.cox.clocks.pc1[[clock]][[disease]])$table[1, ]
  res.cox.clocks.pc1.zph.global[[clock]][[disease]]    <- cox.zph(res.cox.clocks.pc1[[clock]][[disease]])$table["GLOBAL", ]
  res.cox.clocks.pc1.2.zph.local[[clock]][[disease]]   <- cox.zph(res.cox.clocks.pc1.2[[clock]][[disease]])$table[1, ]
  res.cox.clocks.pc1.2.zph.global[[clock]][[disease]]  <- cox.zph(res.cox.clocks.pc1.2[[clock]][[disease]])$table["GLOBAL", ]
  res.cox.clocks.pc1.3.zph.local[[clock]][[disease]]   <- cox.zph(res.cox.clocks.pc1.3[[clock]][[disease]])$table[1, ]
  res.cox.clocks.pc1.3.zph.global[[clock]][[disease]]  <- cox.zph(res.cox.clocks.pc1.3[[clock]][[disease]])$table["GLOBAL", ]
  res.cox.clocks.pc1.4.zph.local[[clock]][[disease]]   <- cox.zph(res.cox.clocks.pc1.4[[clock]][[disease]])$table[1, ]
  res.cox.clocks.pc1.4.zph.global[[clock]][[disease]]  <- cox.zph(res.cox.clocks.pc1.4[[clock]][[disease]])$table["GLOBAL", ]
  res.cox.clocks.pc1.5.zph.local[[clock]][[disease]]   <- cox.zph(res.cox.clocks.pc1.5[[clock]][[disease]])$table[1, ]
  res.cox.clocks.pc1.5.zph.global[[clock]][[disease]]  <- cox.zph(res.cox.clocks.pc1.5[[clock]][[disease]])$table["GLOBAL", ]
  res.cox.clocks.pc1.6.zph.local[[clock]][[disease]]   <- cox.zph(res.cox.clocks.pc1.6[[clock]][[disease]])$table[1, ]
  res.cox.clocks.pc1.6.zph.global[[clock]][[disease]]  <- cox.zph(res.cox.clocks.pc1.6[[clock]][[disease]])$table["GLOBAL", ]
  res.cox.clocks.pc1.7.zph.local[[clock]][[disease]]   <- cox.zph(res.cox.clocks.pc1.7[[clock]][[disease]])$table[1, ]
  res.cox.clocks.pc1.7.zph.global[[clock]][[disease]]  <- cox.zph(res.cox.clocks.pc1.7[[clock]][[disease]])$table["GLOBAL", ]
  res.cox.clocks.pc1.8.zph.local[[clock]][[disease]]   <- cox.zph(res.cox.clocks.pc1.8[[clock]][[disease]])$table[1, ]
  res.cox.clocks.pc1.8.zph.global[[clock]][[disease]]  <- cox.zph(res.cox.clocks.pc1.8[[clock]][[disease]])$table["GLOBAL", ]
  res.cox.clocks.pc1.9.zph.local[[clock]][[disease]]   <- cox.zph(res.cox.clocks.pc1.9[[clock]][[disease]])$table[1, ]
  res.cox.clocks.pc1.9.zph.global[[clock]][[disease]]  <- cox.zph(res.cox.clocks.pc1.9[[clock]][[disease]])$table["GLOBAL", ]
  res.cox.clocks.pc1.10.zph.local[[clock]][[disease]]  <- cox.zph(res.cox.clocks.pc1.10[[clock]][[disease]])$table[1, ]
  res.cox.clocks.pc1.10.zph.global[[clock]][[disease]] <- cox.zph(res.cox.clocks.pc1.10[[clock]][[disease]])$table["GLOBAL", ]
  res.cox.clocks.pc1.11.zph.local[[clock]][[disease]]  <- cox.zph(res.cox.clocks.pc1.11[[clock]][[disease]])$table[1, ]
  res.cox.clocks.pc1.11.zph.global[[clock]][[disease]] <- cox.zph(res.cox.clocks.pc1.11[[clock]][[disease]])$table["GLOBAL", ]
  
  res.cox.pc1.zph.local[[clock]][[disease]]            <- cox.zph(res.cox.pc1[[clock]][[disease]])$table[1, ]
  res.cox.pc1.zph.global[[clock]][[disease]]           <- cox.zph(res.cox.pc1[[clock]][[disease]])$table["GLOBAL", ]
  res.cox.pc2.zph.local[[clock]][[disease]]            <- cox.zph(res.cox.pc2[[clock]][[disease]])$table[1, ]
  res.cox.pc2.zph.global[[clock]][[disease]]           <- cox.zph(res.cox.pc2[[clock]][[disease]])$table["GLOBAL", ]
  res.cox.pc3.zph.local[[clock]][[disease]]            <- cox.zph(res.cox.pc3[[clock]][[disease]])$table[1, ]
  res.cox.pc3.zph.global[[clock]][[disease]]           <- cox.zph(res.cox.pc3[[clock]][[disease]])$table["GLOBAL", ]
  res.cox.pc4.zph.local[[clock]][[disease]]            <- cox.zph(res.cox.pc4[[clock]][[disease]])$table[1, ]
  res.cox.pc4.zph.global[[clock]][[disease]]           <- cox.zph(res.cox.pc4[[clock]][[disease]])$table["GLOBAL", ]
  res.cox.pc5.zph.local[[clock]][[disease]]            <- cox.zph(res.cox.pc5[[clock]][[disease]])$table[1, ]
  res.cox.pc5.zph.global[[clock]][[disease]]           <- cox.zph(res.cox.pc5[[clock]][[disease]])$table["GLOBAL", ]
  res.cox.pc6.zph.local[[clock]][[disease]]            <- cox.zph(res.cox.pc6[[clock]][[disease]])$table[1, ]
  res.cox.pc6.zph.global[[clock]][[disease]]           <- cox.zph(res.cox.pc6[[clock]][[disease]])$table["GLOBAL", ]
  res.cox.pc7.zph.local[[clock]][[disease]]            <- cox.zph(res.cox.pc7[[clock]][[disease]])$table[1, ]
  res.cox.pc7.zph.global[[clock]][[disease]]           <- cox.zph(res.cox.pc7[[clock]][[disease]])$table["GLOBAL", ]
  res.cox.pc8.zph.local[[clock]][[disease]]            <- cox.zph(res.cox.pc8[[clock]][[disease]])$table[1, ]
  res.cox.pc8.zph.global[[clock]][[disease]]           <- cox.zph(res.cox.pc8[[clock]][[disease]])$table["GLOBAL", ]
  res.cox.pc9.zph.local[[clock]][[disease]]            <- cox.zph(res.cox.pc9[[clock]][[disease]])$table[1, ]
  res.cox.pc9.zph.global[[clock]][[disease]]           <- cox.zph(res.cox.pc9[[clock]][[disease]])$table["GLOBAL", ]
  res.cox.pc10.zph.local[[clock]][[disease]]           <- cox.zph(res.cox.pc10[[clock]][[disease]])$table[1, ]
  res.cox.pc10.zph.global[[clock]][[disease]]          <- cox.zph(res.cox.pc10[[clock]][[disease]])$table["GLOBAL", ]
  res.cox.pc11.zph.local[[clock]][[disease]]           <- cox.zph(res.cox.pc11[[clock]][[disease]])$table[1, ]
  res.cox.pc11.zph.global[[clock]][[disease]]          <- cox.zph(res.cox.pc11[[clock]][[disease]])$table["GLOBAL", ]
  
}

  # Cleanup
  rm(final_data)
  gc() 
}


process_model_results <- function(model_list, zph_local, zph_global, tag = "model") {
  
  # --- Extract local p-values
  local_p <- lapply(names(zph_local), function(clock) {
    tmp <- zph_local[[clock]]
    df <- data.frame(do.call('rbind', tmp))
    df$Clock <- clock
    df$disease <- rownames(df)
    df
  })
  local_p <- do.call(rbind, local_p)
  names(local_p)[1:3] <- c("chisq_local", "df_local", "local_p")
  
  # --- Extract global p-values
  global_p <- lapply(names(zph_global), function(clock) {
    tmp <- zph_global[[clock]]
    df <- data.frame(do.call('rbind', tmp))
    df$Clock <- clock
    df$disease <- rownames(df)
    df
  })
  global_p <- do.call(rbind, global_p)
  names(global_p)[1:3] <- c("chisq_global", "df_global", "global_p")
  
  # --- Extract summary output from Cox models
  extract_coefficients <- function(models_list) {
    results <- list()
    for (disease in names(models_list)) {
      disease_coefficients <- list()
      for (clock in names(models_list[[disease]])) {
        model <- models_list[[disease]][[clock]]
        coeff <- summary(model)$coefficients[1, ]
        n_event <- summary(model)$nevent
        n_total <- summary(model)$n
        disease_coefficients[[clock]] <- c(coeff, n_event, n_total)
      }
      results[[disease]] <- do.call(rbind, disease_coefficients)
    }
    final_results <- do.call(rbind, lapply(names(results), function(disease) {
      cbind(Disease = disease, results[[disease]])
    }))
    return(as.data.frame(final_results))
  }
  
  out <- extract_coefficients(model_list)
  names(out)[c(1,3,4,5,6,7,8)] <- c("Clock","HR", "SE","Z", "P","N_cases","N_total")
  out$disease <- sub("\\..*", "", rownames(out))
  
  out <- out %>%
    mutate(across(c("HR", "coef", "SE", "P", "Z", "N_cases", "N_total"), as.numeric))
  
  out$LCI <- exp(out$coef - 1.96*out$SE)
  out$UCI <- exp(out$coef + 1.96*out$SE)
  
  out1 <- out[,c("Clock","disease","N_cases","N_total","HR","LCI","UCI","P","Z")]
  out1$df <- 1
  
  # --- Merge all results
  out2 <- merge(out1, local_p, by=c("Clock", "disease")) 
  out3 <- merge(out2, global_p, by=c("Clock", "disease")) 

  
  # --- Optional: add tag
  out3$model <- tag
  return(out3)
}




results_list <- list(

  clocks        = list(model=res.cox.clocks, zph_local=res.cox.clocks.zph.local, zph_global=res.cox.clocks.zph.global),
  
  clocks.pc1    = list(model=res.cox.clocks.pc1, zph_local=res.cox.clocks.pc1.zph.local, zph_global=res.cox.clocks.pc1.zph.global),
  clocks.pc1.2  = list(model=res.cox.clocks.pc1.2, zph_local=res.cox.clocks.pc1.2.zph.local, zph_global=res.cox.clocks.pc1.2.zph.global),
  clocks.pc1.3  = list(model=res.cox.clocks.pc1.3, zph_local=res.cox.clocks.pc1.3.zph.local, zph_global=res.cox.clocks.pc1.3.zph.global),
  clocks.pc1.4  = list(model=res.cox.clocks.pc1.4, zph_local=res.cox.clocks.pc1.4.zph.local, zph_global=res.cox.clocks.pc1.4.zph.global),
  clocks.pc1.5  = list(model=res.cox.clocks.pc1.5, zph_local=res.cox.clocks.pc1.5.zph.local, zph_global=res.cox.clocks.pc1.5.zph.global),
  clocks.pc1.6  = list(model=res.cox.clocks.pc1.6, zph_local=res.cox.clocks.pc1.6.zph.local, zph_global=res.cox.clocks.pc1.6.zph.global),
  clocks.pc1.7  = list(model=res.cox.clocks.pc1.7, zph_local=res.cox.clocks.pc1.7.zph.local, zph_global=res.cox.clocks.pc1.7.zph.global),
  clocks.pc1.8  = list(model=res.cox.clocks.pc1.8, zph_local=res.cox.clocks.pc1.8.zph.local, zph_global=res.cox.clocks.pc1.8.zph.global),
  clocks.pc1.9  = list(model=res.cox.clocks.pc1.9, zph_local=res.cox.clocks.pc1.9.zph.local, zph_global=res.cox.clocks.pc1.9.zph.global),
  clocks.pc1.10 = list(model=res.cox.clocks.pc1.10, zph_local=res.cox.clocks.pc1.10.zph.local, zph_global=res.cox.clocks.pc1.10.zph.global),
  clocks.pc1.11 = list(model=res.cox.clocks.pc1.11, zph_local=res.cox.clocks.pc1.11.zph.local, zph_global=res.cox.clocks.pc1.11.zph.global),
  
  pc1           = list(model=res.cox.pc1, zph_local=res.cox.pc1.zph.local, zph_global=res.cox.pc1.zph.global),
  pc2           = list(model=res.cox.pc2, zph_local=res.cox.pc2.zph.local, zph_global=res.cox.pc2.zph.global),
  pc3           = list(model=res.cox.pc3, zph_local=res.cox.pc3.zph.local, zph_global=res.cox.pc3.zph.global),
  pc4           = list(model=res.cox.pc4, zph_local=res.cox.pc4.zph.local, zph_global=res.cox.pc4.zph.global),
  pc5           = list(model=res.cox.pc5, zph_local=res.cox.pc5.zph.local, zph_global=res.cox.pc5.zph.global),
  pc6           = list(model=res.cox.pc6, zph_local=res.cox.pc6.zph.local, zph_global=res.cox.pc6.zph.global),
  pc7           = list(model=res.cox.pc7, zph_local=res.cox.pc7.zph.local, zph_global=res.cox.pc7.zph.global),
  pc8           = list(model=res.cox.pc8, zph_local=res.cox.pc8.zph.local, zph_global=res.cox.pc8.zph.global),
  pc9           = list(model=res.cox.pc9, zph_local=res.cox.pc9.zph.local, zph_global=res.cox.pc9.zph.global),
  pc10          = list(model=res.cox.pc10, zph_local=res.cox.pc10.zph.local, zph_global=res.cox.pc10.zph.global),
  pc11          = list(model=res.cox.pc11, zph_local=res.cox.pc11.zph.local, zph_global=res.cox.pc11.zph.global)
)


all_results <- do.call(rbind, lapply(names(results_list), function(tag) {
  model_info <- results_list[[tag]]
  process_model_results(model_info$model, model_info$zph_local, model_info$zph_global, tag)
}))


# Please change file path. --Thomas
save(all_results, file = "/file_path/Disease_Results.rda")
save(res.an.clocks, file = "/file_path/clock_WBC_LRT.rda")