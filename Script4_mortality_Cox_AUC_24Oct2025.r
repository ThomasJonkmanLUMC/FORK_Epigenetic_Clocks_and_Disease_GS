### read in death data and QCd clock data ###
mortality_data<- read.csv("/Cluster_Filespace/Marioni_Group/Riccardo/Christos_17Oct2024/Data/2024-04-19_deaths.csv")
d <- readRDS(file="/Cluster_Filespace/Marioni_Group/Riccardo/Christos_17Oct2024/Data/Merged_clock_data_17Oct2024.RDS")

### create 10-year death phenotype and tte ###
temp <- which(d$id %in% mortality_data$id)
d$dead <- 0
d$dead[temp] <- 1
temp2 <- which(d$dead==1 & d$t_censor>10)  
d$dead[temp2] <- 0
d$t_censor[d$t_censor>10] <- 10

# List of clocks to loop through
clocks <- c("Hannum", "Horvathv1", "Zhang_Acc", "PhenoAge", "DNAmGrimAge.1", "DunedinPACE")

# objects to store results
# Just like Script2, we test the effect of clocks alone, clocks adjusted for WBC PCs, and WBC PCs alone.
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

library(survival)

# Prepare the cox regression models. --Thomas
model_formula_clocks <- as.formula(paste0("Surv(t_censor, dead) ~ scale(clock) + age + bmi + years + pack_years + units + rank + sex"))

model_formula_clocks.pc1 <- as.formula(paste0("Surv(t_censor, dead) ~ scale(clock) + age + bmi + years + pack_years + units + rank + sex + PC1"))
model_formula_clocks.pc1.2 <- as.formula(paste0("Surv(t_censor, dead) ~ scale(clock) + age + bmi + years + pack_years + units + rank + sex + PC1 + PC2"))
model_formula_clocks.pc1.3 <- as.formula(paste0("Surv(t_censor, dead) ~ scale(clock) + age + bmi + years + pack_years + units + rank + sex + PC1 + PC2 + PC3"))
model_formula_clocks.pc1.4 <- as.formula(paste0("Surv(t_censor, dead) ~ scale(clock) + age + bmi + years + pack_years + units + rank + sex + PC1 + PC2 + PC3 + PC4"))
model_formula_clocks.pc1.5 <- as.formula(paste0("Surv(t_censor, dead) ~ scale(clock) + age + bmi + years + pack_years + units + rank + sex + PC1 + PC2 + PC3 + PC4 + PC5"))
model_formula_clocks.pc1.6 <- as.formula(paste0("Surv(t_censor, dead) ~ scale(clock) + age + bmi + years + pack_years + units + rank + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6"))
model_formula_clocks.pc1.7 <- as.formula(paste0("Surv(t_censor, dead) ~ scale(clock) + age + bmi + years + pack_years + units + rank + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7"))
model_formula_clocks.pc1.8 <- as.formula(paste0("Surv(t_censor, dead) ~ scale(clock) + age + bmi + years + pack_years + units + rank + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8"))
model_formula_clocks.pc1.9 <- as.formula(paste0("Surv(t_censor, dead) ~ scale(clock) + age + bmi + years + pack_years + units + rank + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9"))
model_formula_clocks.pc1.10 <- as.formula(paste0("Surv(t_censor, dead) ~ scale(clock) + age + bmi + years + pack_years + units + rank + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))
model_formula_clocks.pc1.11 <- as.formula(paste0("Surv(t_censor, dead) ~ scale(clock) + age + bmi + years + pack_years + units + rank + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11"))

model_formula_pc1 <- as.formula(paste0("Surv(t_censor, dead) ~ PC1 + age + bmi + years + pack_years + units + rank + sex"))
model_formula_pc2 <- as.formula(paste0("Surv(t_censor, dead) ~ PC2 + age + bmi + years + pack_years + units + rank + sex"))
model_formula_pc3 <- as.formula(paste0("Surv(t_censor, dead) ~ PC3 + age + bmi + years + pack_years + units + rank + sex"))
model_formula_pc4 <- as.formula(paste0("Surv(t_censor, dead) ~ PC4 + age + bmi + years + pack_years + units + rank + sex"))
model_formula_pc5 <- as.formula(paste0("Surv(t_censor, dead) ~ PC5 + age + bmi + years + pack_years + units + rank + sex"))
model_formula_pc6 <- as.formula(paste0("Surv(t_censor, dead) ~ PC6 + age + bmi + years + pack_years + units + rank + sex"))
model_formula_pc7 <- as.formula(paste0("Surv(t_censor, dead) ~ PC7 + age + bmi + years + pack_years + units + rank + sex"))
model_formula_pc8 <- as.formula(paste0("Surv(t_censor, dead) ~ PC8 + age + bmi + years + pack_years + units + rank + sex"))
model_formula_pc9 <- as.formula(paste0("Surv(t_censor, dead) ~ PC9 + age + bmi + years + pack_years + units + rank + sex"))
model_formula_pc10 <- as.formula(paste0("Surv(t_censor, dead) ~ PC10 + age + bmi + years + pack_years + units + rank + sex"))
model_formula_pc11 <- as.formula(paste0("Surv(t_censor, dead) ~ PC11 + age + bmi + years + pack_years + units + rank + sex"))


# Run cox regression analyses for each clock.
CoxRegression <- function(res.cox, model_formula, res.cox.zph.local, res.cox.zph.global){
    
  for(clock in clocks){
        
    d$clock <<- d[,clock]
    res.cox[[clock]] <- coxph(model_formula, data = d)
    res.cox.zph.local[[clock]] <- cox.zph(res.cox[[clock]])$table[1,] 
    res.cox.zph.global[[clock]] <- cox.zph(res.cox[[clock]])$table["GLOBAL",]
  }
  
  ## prop hazards (local) ##
  local_p = list()
  
  for(clock in names(res.cox.zph.local)){
  local_p[[clock]] <- res.cox.zph.local[[clock]]
  }
  
  local_p = as.data.frame(do.call('rbind', local_p))
  local_p$Clock <- as.character(rownames(local_p))
  names(local_p)[1:3] <- c("chisq_local", "df_local", "local_p")
  
  
  ## prop hazards (global) ##
  global_p = list()
  
  for(clock in names(res.cox.zph.global)){
  global_p[[clock]] <- res.cox.zph.global[[clock]]
  }
  
  global_p = as.data.frame(do.call('rbind', global_p))
  global_p$Clock <- as.character(rownames(global_p))
  names(global_p)[1:3] <- c("chisq_global", "df_global", "global_p")
  
  
  ## Extract summary output from Cox models for plotting 
  extract_coefficients <- function(models_list) {
    # Initialize an empty list to store results
    results <- list()
    mort_coefficients <- list()
     
      # Iterate over each epigenetic clock model
      for (clock in names(models_list)) {
        model <- models_list[[clock]]
        
          # Extract the first row of coefficients
          coeff <- summary(model)$coefficients[1,]
  	      n_event <- summary(model)$nevent
  	      n_total <- summary(model)$n
          mort_coefficients[[clock]] <- c(coeff, n_event, n_total)
      }
      
      # Combine the coefficients into a data frame for the current disease
      results <- do.call(rbind, mort_coefficients)  
  
    return(as.data.frame(results))
  }
  
  out <- extract_coefficients(res.cox)
  names(out) <- c("logHR","HR", "SE","Z","P","N_cases","N_total")
  out$Clock <- rownames(out)
  
  out$LCI <- exp(out$logHR - 1.96*out$SE)
  out$UCI <- exp(out$logHR + 1.96*out$SE)
  
  out1 <- out[,c("Clock","N_cases","N_total","HR","LCI","UCI","Z","P")]
  
  out2 <- merge(out1, local_p, by=c("Clock")) 
  out3 <- merge(out2, global_p, by=c("Clock")) 
}


#Run Cox-regression models and store results. --Thomas
out.clocks <- CoxRegression(res.cox.clocks, model_formula_clocks, res.cox.clocks.zph.local, res.cox.clocks.zph.global)

out.clocks.pc1 <- CoxRegression(res.cox.clocks.pc1, model_formula_clocks.pc1, res.cox.clocks.pc1.zph.local, res.cox.clocks.pc1.zph.global)
out.clocks.pc1.2 <- CoxRegression(res.cox.clocks.pc1.2, model_formula_clocks.pc1.2, res.cox.clocks.pc1.2.zph.local, res.cox.clocks.pc1.2.zph.global)
out.clocks.pc1.3 <- CoxRegression(res.cox.clocks.pc1.3, model_formula_clocks.pc1.3, res.cox.clocks.pc1.3.zph.local, res.cox.clocks.pc1.3.zph.global)
out.clocks.pc1.4 <- CoxRegression(res.cox.clocks.pc1.4, model_formula_clocks.pc1.4, res.cox.clocks.pc1.4.zph.local, res.cox.clocks.pc1.4.zph.global)
out.clocks.pc1.5 <- CoxRegression(res.cox.clocks.pc1.5, model_formula_clocks.pc1.5, res.cox.clocks.pc1.5.zph.local, res.cox.clocks.pc1.5.zph.global)
out.clocks.pc1.6 <- CoxRegression(res.cox.clocks.pc1.6, model_formula_clocks.pc1.6, res.cox.clocks.pc1.6.zph.local, res.cox.clocks.pc1.6.zph.global)
out.clocks.pc1.7 <- CoxRegression(res.cox.clocks.pc1.7, model_formula_clocks.pc1.7, res.cox.clocks.pc1.7.zph.local, res.cox.clocks.pc1.7.zph.global)
out.clocks.pc1.8 <- CoxRegression(res.cox.clocks.pc1.8, model_formula_clocks.pc1.8, res.cox.clocks.pc1.8.zph.local, res.cox.clocks.pc1.8.zph.global)
out.clocks.pc1.9 <- CoxRegression(res.cox.clocks.pc1.9, model_formula_clocks.pc1.9, res.cox.clocks.pc1.9.zph.local, res.cox.clocks.pc1.9.zph.global)
out.clocks.pc1.10 <- CoxRegression(res.cox.clocks.pc1.10, model_formula_clocks.pc1.10, res.cox.clocks.pc1.10.zph.local, res.cox.clocks.pc1.10.zph.global)
out.clocks.pc1.11 <- CoxRegression(res.cox.clocks.pc1.11, model_formula_clocks.pc1.11, res.cox.clocks.pc1.11.zph.local, res.cox.clocks.pc1.11.zph.global)

out.pc1 <- CoxRegression(res.cox.pc1, model_formula_pc1, res.cox.pc1.zph.local, res.cox.pc1.zph.global)
out.pc2 <- CoxRegression(res.cox.pc2, model_formula_pc2, res.cox.pc2.zph.local, res.cox.pc2.zph.global)
out.pc3 <- CoxRegression(res.cox.pc3, model_formula_pc3, res.cox.pc3.zph.local, res.cox.pc3.zph.global)
out.pc4 <- CoxRegression(res.cox.pc4, model_formula_pc4, res.cox.pc4.zph.local, res.cox.pc4.zph.global)
out.pc5 <- CoxRegression(res.cox.pc5, model_formula_pc5, res.cox.pc5.zph.local, res.cox.pc5.zph.global)
out.pc6 <- CoxRegression(res.cox.pc6, model_formula_pc6, res.cox.pc6.zph.local, res.cox.pc6.zph.global)
out.pc7 <- CoxRegression(res.cox.pc7, model_formula_pc7, res.cox.pc7.zph.local, res.cox.pc7.zph.global)
out.pc8 <- CoxRegression(res.cox.pc8, model_formula_pc8, res.cox.pc8.zph.local, res.cox.pc8.zph.global)
out.pc9 <- CoxRegression(res.cox.pc9, model_formula_pc9, res.cox.pc9.zph.local, res.cox.pc9.zph.global)
out.pc10 <- CoxRegression(res.cox.pc10, model_formula_pc10, res.cox.pc10.zph.local, res.cox.pc10.zph.global)
out.pc11 <- CoxRegression(res.cox.pc11, model_formula_pc11, res.cox.pc11.zph.local, res.cox.pc11.zph.global)

out.pcs <- rbind(out.pc1[1,], out.pc2[1,], out.pc3[1,], out.pc4[1,], out.pc5[1,], out.pc6[1,], out.pc7[1,], out.pc8[1,], out.pc9[1,], out.pc10[1,], out.pc11[1,])
names(out.pcs)[1] <- "PC"
out.pcs$PC <- paste0("PC", 1:11)

#Calculate the level of evidence that WBC composition affects the clock-disease association. --Thomas
out.an.clocks <- setNames(vector("list", length(clocks)), clocks)
for(clock in clocks){
  
  d$clock <- d[,clock]
  fit1 <- coxph(model_formula_clocks, data = d)
  fit2 <- coxph(model_formula_clocks.pc1.11, data = d)
  out.an.clocks[[clock]] <- anova(fit1, fit2, test = "LRT")$`Pr(>|Chi|)`[2]
}

all_results <- list(
  clocks = out.clocks,
  
  clocks.pc1 = out.clocks.pc1,
  clocks.pc1.2 = out.clocks.pc1.2,
  clocks.pc1.3 = out.clocks.pc1.3,
  clocks.pc1.4 = out.clocks.pc1.4,
  clocks.pc1.5 = out.clocks.pc1.5,
  clocks.pc1.6 = out.clocks.pc1.6,
  clocks.pc1.7 = out.clocks.pc1.7,
  clocks.pc1.8 = out.clocks.pc1.8,
  clocks.pc1.9 = out.clocks.pc1.9,
  clocks.pc1.10 = out.clocks.pc1.10,
  clocks.pc1.11 = out.clocks.pc1.11,
  
  pcs = out.pcs,
  an.clocks = out.an.clocks
)

# Please change file path. --Thomas
save(all_results, file = "/file_path/Mortality_results.rda")