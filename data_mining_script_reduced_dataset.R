#### MC_004: Renal Cancer public dataset re-analysis (Renal Cancer Grant) ####

## Data mining script ####

# Data downloaded from: https://cptac-data-portal.georgetown.edu/cptac/s/S050
# and other associated sources, such as:
# https://portal.gdc.cancer.gov/projects/CPTAC-3

## Load packages ####

library(tidyverse)
library(readxl)

## Load datasets ####

## * 1. Proteomics-related datasets ####

# This dataset contains log2 values of expression of each protein (row) per sample (column)
prot_tmt10 <- read_tsv(file = here::here("Test_data/CPTAC_ccRCC_data_2019/CPTAC3_Clear_Cell_Renal_Cell_Carcinoma_Proteome.tmt10.tsv"))

# This dataset contains the sample-to-TMT-channel-mapping   
samp_chnl <- read_excel(path = here::here("Test_data/CPTAC_ccRCC_data_2019/S044_CPTAC_CCRCC_Discovery_Cohort_Specimens_r1_Sept2018.xlsx"))  

## Data wrangling ####

## Proteomics quantitative data ####

# namesof the columns
nam <- names(prot_tmt10)
# select colunm names that have "Unshared"
xx <- str_subset(nam, "Unshared")

prot_quant <- select(prot_tmt10,
                     setdiff(nam,xx), -contains("QC"), -contains("NCI7"),
                     -c(Authority, Description, Organism, Chromosome, Locus,
                        NCBIGeneID))

names(prot_quant) <- str_sub(names(prot_quant), 1, 13)

length(names(prot_quant))

# * 1 Proteomics: Separate the columns according to tumor or not ----

sample_data <- filter(samp_chnl,
                      `Aliquot ID` %in% names(prot_quant),
                      duplicated(`Aliquot ID`) == FALSE)

selected_patients <- dplyr::filter(sample_data,
                                   ParticipantID %in% unique(ParticipantID)[c(20:25)]) %>%
                     filter(ParticipantID != "C3N-01180")

# Keep only patients with two samples (Tumor vs Normal) (Experimental design matrix for Limma)

sample_datafil <- group_by(sample_data, ParticipantID) %>%
                  filter(n() == 2) %>% arrange(ParticipantID) %>%
                  select(`Aliquot ID`, ParticipantID, Group) %>% ungroup()


sample_datafilsel <- group_by(selected_patients, ParticipantID) %>%
                     filter(n() == 2) %>% arrange(ParticipantID) %>%
                     select(`Aliquot ID`, ParticipantID, Group) %>% ungroup()


# Create matrix for Limma ----
tomat <- as.data.frame(prot_quant)

row.names(tomat) <- prot_quant$Gene

tomat <- select(tomat, -Gene)

tomatsli <- tomat[-c(1:3),]
      
data_protsel <- select(tomatsli,
                       sample_datafil$`Aliquot ID`) 

selected_patdata <- select(tomatsli,
                           sample_datafilsel$`Aliquot ID`) 



save(selected_patdata,
     file = "CPTAC_CCRCC_Proteome_tabular_reduced.Rda")

matrix_all <- as.matrix(data_protsel)

matrix_all[1:5,1:5]

class(matrix_all)


### Create design matrix and run Limma ----

#from sample_datafil

patientfct <- factor(sample_datafil$ParticipantID)
group <- factor(sample_datafil$Group, levels = c("Normal","Tumor"))

design <- model.matrix(~patientfct+group)

designCPTACCRCCProteomics <- design

save(designCPTACCRCCProteomics,
     file  = "exp_design_CPTACCRCC_Proteomics.Rda")


## Design matrix reducen data ----

patientfct2 <- factor(sample_datafilsel$ParticipantID)
group2 <- factor(sample_datafilsel$Group, levels = c("Normal","Tumor"))

design2 <- model.matrix(~patientfct2+group2)

designCPTACCRCCProteomics_reduced <- design2

save(designCPTACCRCCProteomics_reduced,
     file  = "exp_design_CPTACCRCC_Proteomics_reduced.Rda")



## Run Limma ----

library(limma)

fit <- lmFit(matrix_all, design = design)
fit <- eBayes(fit)

output_limma <- topTable(fit,
                         number = dim(fit$t)[1],
                         coef="groupTumor")

gene_id <- rownames(output_limma)

df_fitable <- cbind(gene_id, output_limma)


### Get filter Limma output by proteases ----

fitable_proteases <- dplyr::filter(df_fitable,
                                   gene_id %in% proteases$Gene)

write_tsv(x = fitable_proteases,
          path = here::here("limmaFit_output_filtered_by_Proteases_ccRCC_CPTAC_NormalvsTumor.tsv"))


diff_expr_proteases_modpval <- dplyr::filter(fitable_proteases, 
                                             P.Value <= 0.05,
                                             abs(logFC) >= log2(1.2)) %>% pull(gene_id)

diff_expr_proteases_adjpval <- dplyr::filter(fitable_proteases, 
                                             adj.P.Val <= 0.05,
                                             abs(logFC) >= log2(1.2)) %>% pull(gene_id) 

## Plot volcano ----  

### VOLCANO PROTEASES ####

library(ggrepel)

volcan_proteases_moderate <- ggplot(fitable_proteases,
                                    aes(x = logFC, y = -log10(P.Value))) + 
      geom_point() + geom_hline(yintercept = -log10(0.05), colour="#990000", linetype="dashed") + 
      geom_vline(xintercept = -log2(1.2), colour="#990000", linetype="dashed") + 
      geom_vline(xintercept = log2(1.2), colour="#990000", linetype="dashed") +
      xlim(c(-2.5, 2.5)) + ylim(c(0, 10)) +
      ylab("-Log10(Moderate p-value)") + 
      xlab("Log2(Fold-change)") +
      geom_text_repel(aes(label=ifelse(gene_id %in% diff_expr_proteases_modpval,
                                       as.character(gene_id),'')),
                      size = 3.2, fontface="bold")+
      ggtitle("Differentially expressed proteases (Limma). ccRCC_CPTAC-Proteomics",
              subtitle = "Up-regulated are more expressed in Tumor" )+
      theme(axis.text.x = element_text(hjust = 0.5),
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=1.5),
            axis.title=element_text(size=10,face="bold"))

volcan_proteases_moderate


volcan_proteases_adjustP <- ggplot(fitable_proteases,
                                   aes(x = logFC, y = -log10(adj.P.Val))) + 
      geom_point() +
      geom_hline(yintercept = -log10(0.05), colour="#990000", linetype="dashed") + 
      geom_vline(xintercept = -log2(1.2), colour="#990000", linetype="dashed") + 
      geom_vline(xintercept = log2(1.2), colour="#990000", linetype="dashed") +
      xlim(c(-2.5, 2.5)) + ylim(c(0, 10)) +
      ylab("-Log10(Adjusted p-value)") + 
      xlab("Log2(Fold-change)") +
      geom_text_repel(aes(label=ifelse(gene_id %in% diff_expr_proteases_adjpval,
                                       as.character(gene_id),'')),
                      size = 3.2, fontface="bold")+
      ggtitle("Differentially expressed proteases (Limma). ccRCC_CPTAC-Proteomics",
              subtitle = "Up-regulated are more expressed in Tumor" )+
      theme(axis.text.x = element_text(hjust = 0.5),
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=1.5),
            axis.title=element_text(size=10,face="bold"))


volcan_proteases_adjustP



# Define 



# Split dataset according to Tumor or Normal samples  

sample_tum <- filter(sample_datafil,
                     Group == "Tumor")

sample_nor <- filter(sample_datafil,
                     Group == "Normal")

# Get columns for Tumor samples 



length(sample_data$`Aliquot ID`)
duplicated(sample_data$`Aliquot ID`)

length(unique(sample_data$ParticipantID))

which(!duplicated(sample_data$ParticipantID))



