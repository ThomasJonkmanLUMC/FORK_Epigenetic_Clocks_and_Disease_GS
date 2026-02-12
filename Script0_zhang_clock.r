# .libPaths("~/RSC_BIOS/Users/tjonkman/cellcounts/test/Marioni_scripts/Packages")

# Tested on fresh R 4.3.1 install
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("methylclock")
library(methylclock)

#Load beta-values and calculate the super-accurate Zhang clock.
betas <- read.csv("path_to_beta_values")
zhang <- as.data.frame(DNAmAge(assay(betas), clocks = "EN"))

#Append the accurate Zhang-clock to the clock data.
clock_data <- read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/Standardised_Data/gs_clocks.csv")

#Please match samples appropriately if necessary before merging data.
clock_data$Zhang_Acc <- cbind(clock)

#Please plot DNAmAge against calendar age as a sanity check. Age information might be in a separate file.
cor(clock_data$Zhang_Acc, clock_data$Age)
library(ggplot2)
p <- ggplot(clock_data, aes(x = Age, y = Zhang_Acc)) +
  geom_abline(slope = 1, intercept = 0, colour = "red", linetype = "dashed") +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw()
p

#Save to disk (please change file path).
saveRDS(clock_data, file = "/file_path/appended_clock_data.rds")
save(p, file = "/file_path/Zhang_Acc_vs_calendar_age.rda")

# SessionInfo paste
# R version 4.3.1 (2023-06-16)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Rocky Linux 8.10 (Green Obsidian)
# 
# Matrix products: default
# BLAS/LAPACK: /usr/lib64/libopenblas-r0.3.15.so;  LAPACK version 3.9.0
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: Europe/Amsterdam
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggplot2_4.0.2               SummarizedExperiment_1.32.0 Biobase_2.62.0              GenomicRanges_1.54.1        GenomeInfoDb_1.38.8         IRanges_2.36.0              S4Vectors_0.40.2           
# [8] BiocGenerics_0.48.1         MatrixGenerics_1.14.0       matrixStats_1.5.0           methylclock_1.8.0           quadprog_1.5-8              devtools_2.4.6              usethis_3.2.1              
# [15] methylclockData_1.10.0      futile.logger_1.4.9        
# 
# loaded via a namespace (and not attached):
#   [1] splines_4.3.1                 later_1.4.5                   BiocIO_1.12.0                 bitops_1.0-9                  filelock_1.0.3                tibble_3.3.1                  ggpp_0.6.0                   
# [8] preprocessCore_1.64.0         graph_1.80.0                  xts_0.14.1                    XML_3.99-0.22                 lifecycle_1.0.5               httr2_1.2.2                   rstatix_0.7.3                
# [15] base64_2.0.2                  lattice_0.21-8                MASS_7.3-60                   OrganismDbi_1.44.0            scrime_1.3.7                  backports_1.5.0               magrittr_2.0.4               
# [22] limma_3.58.1                  minfi_1.48.0                  yaml_2.3.12                   remotes_2.5.0                 httpuv_1.6.16                 otel_0.2.0                    doRNG_1.8.6.3                
# [29] askpass_1.2.1                 sessioninfo_1.2.3             pkgbuild_1.4.8                RUnit_0.4.33.1                DBI_1.2.3                     RColorBrewer_1.1-3            abind_1.4-8                  
# [36] pkgload_1.5.0                 zlibbioc_1.48.2               purrr_1.2.1                   RCurl_1.98-1.17               rappdirs_0.3.4                AnnotationHubData_1.32.1      GenomeInfoDbData_1.2.11      
# [43] AnnotationForge_1.44.0        ggpmisc_0.6.3                 genefilter_1.84.0             MatrixModels_0.5-4            annotate_1.80.0               DelayedMatrixStats_1.24.0     codetools_0.2-19             
# [50] DelayedArray_0.28.0           xml2_1.5.2                    tidyselect_1.2.1              farver_2.1.2                  beanplot_1.3.1                BiocFileCache_2.10.2          dynamicTreeCut_1.63-1        
# [57] illuminaio_0.44.0             GenomicAlignments_1.38.2      jsonlite_2.0.0                multtest_2.58.0               ellipsis_0.3.2                Formula_1.2-5                 iterators_1.0.14             
# [64] survival_3.5-5                foreach_1.5.2                 tools_4.3.1                   progress_1.2.3                stringdist_0.9.17             Rcpp_1.1.1                    glue_1.8.0                   
# [71] gridExtra_2.3                 SparseArray_1.2.4             BiocBaseUtils_1.4.0           xfun_0.56                     ExperimentHubData_1.28.0      dplyr_1.2.0                   HDF5Array_1.30.1             
# [78] withr_3.0.2                   formatR_1.14                  BiocManager_1.30.27           fastmap_1.2.0                 rhdf5filters_1.14.1           openssl_2.3.4                 SparseM_1.84-2               
# [85] digest_0.6.39                 R6_2.6.1                      mime_0.13                     RPMM_1.25                     biomaRt_2.58.2                RSQLite_2.4.6                 tidyr_1.3.2                  
# [92] generics_0.1.4                PerformanceAnalytics_2.0.8    data.table_1.18.2.1           rtracklayer_1.62.0            prettyunits_1.2.0             httr_1.4.7                    S4Arrays_1.2.1               
# [99] pkgconfig_2.0.3               gtable_0.3.6                  blob_1.3.0                    siggenes_1.76.0               S7_0.2.1                      impute_1.76.0                 XVector_0.42.0               
# [106] htmltools_0.5.9               carData_3.0-6                 RBGL_1.78.0                   scales_1.4.0                  tidyverse_2.0.0               png_0.1-8                     knitr_1.51                   
# [113] lambda.r_1.2.4                rstudioapi_0.18.0             tzdb_0.5.0                    rjson_0.2.23                  nlme_3.1-162                  curl_7.0.0                    bumphunter_1.44.0            
# [120] biocViews_1.70.0              cachem_1.1.0                  zoo_1.8-15                    rhdf5_2.46.1                  stringr_1.6.0                 BiocVersion_3.18.1            parallel_4.3.1               
# [127] AnnotationDbi_1.64.1          restfulr_0.0.16               GEOquery_2.70.0               reshape_0.8.10                pillar_1.11.1                 grid_4.3.1                    vctrs_0.7.1                  
# [134] promises_1.5.0                ggpubr_0.6.2                  car_3.1-5                     dbplyr_2.5.1                  xtable_1.8-4                  cluster_2.1.4                 evaluate_1.0.5               
# [141] readr_2.1.6                   GenomicFeatures_1.54.4        locfit_1.5-9.12               cli_3.6.5                     compiler_4.3.1                futile.options_1.0.1          Rsamtools_2.18.0             
# [148] rngtools_1.5.2                rlang_1.1.7                   crayon_1.5.3                  ggsignif_0.6.4                nor1mix_1.3-3                 mclust_6.1.2                  plyr_1.8.9                   
# [155] fs_1.6.6                      stringi_1.8.7                 BiocParallel_1.36.0           BiocCheck_1.38.2              Biostrings_2.70.3             quantreg_6.1                  Matrix_1.6-5                 
# [162] ExperimentHub_2.10.0          hms_1.1.4                     sparseMatrixStats_1.14.0      bit64_4.6.0-1                 Rhdf5lib_1.24.2               statmod_1.5.1                 KEGGREST_1.42.0              
# [169] shiny_1.12.1                  interactiveDisplayBase_1.40.0 AnnotationHub_3.10.1          broom_1.0.12                  memoise_2.0.1                 bit_4.6.0                     polynom_1.4-1  