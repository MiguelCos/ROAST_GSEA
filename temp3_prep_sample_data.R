## Temp: prepare test dataset (CPTAC: CCRCC Tumor vs Healthy)  

# Reduced sample of 10 patients Tumor vs Healthy  

## Load packages -----

library(tidyverse)

## Load data 

load(file = here::here("CPTAC_CCRCC_Proteome_tabular_reduced.Rda"))

# check 
dim(selected_patdata)

str(selected_patdata)

## Create tabular txt file that simulates limma input ----

cptac_liminput1 <- dplyr::mutate(selected_patdata,
                                 ID = row.names(selected_patdata)) %>% 
                  dplyr::select(ID, 1:10)
row.names(cptac_liminput1) <- NULL

head(cptac_liminput1)
dim(cptac_liminput1)

write.table(cptac_liminput1,
            file = "roast_sample_input_cptac_ccrcc_reduced.txt",
            sep = "\t",
            row.names = FALSE)

## Load experimental design 

load(file  = "exp_design_CPTACCRCC_Proteomics_reduced.Rda")

# Check 
designCPTACCRCCProteomics_reduced

sample_idesign_cptac_ccrcc_reduced <- designCPTACCRCCProteomics_reduced

save(sample_idesign_cptac_ccrcc_reduced,
     file = "sample_idesign_cptac_ccrcc_reduced.Rda")

