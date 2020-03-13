### Temp 5: generate MSstats/MSstatsTMT processed data ----

## Load required packages ----

library(MSstats)
library(MSstatsTMT)
library(tidyverse)

## Load required data ----

# Label-fre data Nr1:  -----

lf_evid <- read.delim(file = here::here("Test_data/Pancreatic_Cancer_BM3448_3470/Short_course/evidence_short_course_pc_001.txt"),
                      sep = "\t",
                      stringsAsFactors = FALSE)

lf_pgr <- read.delim(file = here::here("Test_data/Pancreatic_Cancer_BM3448_3470/Short_course/proteinGroups_short_course_pc_001.txt"),
                     sep = "\t",
                     stringsAsFactors = FALSE)

lf_annot <- read.csv(file = here::here("Test_data/Pancreatic_Cancer_BM3448_3470/Short_course/annotation_sc_pc_001.csv"))


## Label-free data Nr 2: ----

lf2_evid <- read.delim(file = here::here("Test_data/Pancreatic_Cancer_BM3448_3470/Long_course/evidence_long_course_pc_001.txt"),
                       sep = "\t",
                       stringsAsFactors = FALSE)

lf2_pgr <- read.delim(file = here::here("Test_data/Pancreatic_Cancer_BM3448_3470/Long_course/proteinGroups_long_course_pc_001.txt"),
                      sep = "\t",
                      stringsAsFactors = FALSE)

lf2_annot <- read.delim(file = here::here("Test_data/Pancreatic_Cancer_BM3448_3470/Long_course/annotation_lc_pc_001.csv"),
                        sep = ",",
                        stringsAsFactors = FALSE) %>%
             dplyr::filter(Condition != "unknown")


## TMT GBM data (with paired desing) ----

tmt_evid <- read.delim(file = here::here("Test_data/GBM_data/evidence_top_3pept_PIFfilter_0.75.txt"),
                       sep = "\t",
                       stringsAsFactors = FALSE)

tmt_pgr <- read.delim(file = here::here("Test_data/GBM_data/proteinGroups.txt"),
                      sep = "\t",
                      stringsAsFactors = FALSE)

tmt_annot <- read.delim(file = here::here("Test_data/GBM_data/annotation_w_paired_design.tsv"),
                        sep = "\t",
                        stringsAsFactors = FALSE)

## Process MSstats data (label free) ----  

## Label-free data with 4 groups ----

lf_msstats <- MaxQtoMSstatsFormat(evidence = lf_evid,
                                  annotation = lf_annot,
                                  proteinGroups = lf_pgr,
                                  proteinID = "Leading.razor.protein",
                                  useUniquePeptide = TRUE)

lf_procdat <- dataProcess(lf_msstats)


## TMT data with 2 groups and pairwise design ----

tmt_msstats <- MaxQtoMSstatsTMTFormat(evidence = tmt_evid,
                                      proteinGroups = tmt_pgr,
                                      annotation = tmt_annot,
                                      which.proteinid = "Leading.razor.protein",
                                      rmProt_Only.identified.by.site = TRUE)

tmt_procdat <- proteinSummarization(tmt_msstats)

## Label-free data with 2 groups, non-pairwise design -----

lf2_msstats <- MaxQtoMSstatsFormat(evidence = lf2_evid,
                                   annotation = lf2_annot,
                                   proteinGroups = lf2_pgr,
                                   proteinID = "Leading.razor.protein",
                                   useUniquePeptide = TRUE) %>%
            dplyr::filter(Condition != "unknown")

lf2_procdat <- dataProcess(lf2_msstats)

## Reference for the generation of the design object manually ----

### Develop idea of the function for 'automatically' create the design object ----

## REFERENCE

patientfct2 <- factor(sample_datafilsel$ParticipantID)

group2 <- factor(sample_datafilsel$Group, 
                 levels = c("Normal","Tumor"))

patientfct2
group2

design2

##


paired = TRUE
BiologicalReplicates = 5
n_samples = 10
Conditions = c("Normal","Tumor")




samp <- rep("samp", length(Conditions))
seqbioreps <- seq(1,BiologicalReplicates)

bioreps <- factor(paste(rep(seqbioreps, each = length(Conditions)), samp, sep = "_"))

invid <- factor(paste(rep(Conditions, times = BiologicalReplicates)),
                levels = c("Normal","Tumor"))

design <- model.matrix(~bioreps+invid)

design



### Testing how to extract experimental design information and transform into wide format. ----

## Test first with two-gruop data (unpaired) ----

# Working with lf2_procdat 

annotlf1 <-  lf2_annot %>% 
      dplyr::mutate(Condition = factor(Condition,
                                       levels = (unique(Condition)))) %>%
      dplyr::mutate(Sample_ID = paste0(Condition,"_",BioReplicate)) %>% 
      dplyr::select(Sample_ID, BioReplicate, Condition)


lf2proc1 <- lf2_procdat$RunlevelData %>% 
            dplyr::mutate(GROUP_ORIGINAL = factor(GROUP_ORIGINAL,
                                            levels = (unique(GROUP_ORIGINAL))),
                          Protein = as.character(Protein),
                          Sample_ID = paste0(GROUP_ORIGINAL,"_",SUBJECT_ORIGINAL)) %>%
            dplyr::filter(GROUP_ORIGINAL != "unknown",
                          str_detect(Protein,"^Biognosys_pep",  negate = TRUE)) %>% 
            dplyr::select(Sample_ID, Protein, LogIntensities)


rsh_lf2proc1 <- tidyr::pivot_wider(lf2proc1,
                                   names_from = Sample_ID,
                                   values_from = LogIntensities,
                                   values_fn = list(LogIntensities = median)) %>% 
            dplyr::select(ID = Protein, annotlf1$Sample_ID) %>% na.omit()


group <- factor(annotlf1$Condition)

design <- model.matrix(~group)


## Test with TMT data (paired design) ----

Conditions = c("prim", "rec")

annottmt1 <- dplyr::filter(tmt_annot,
                           Condition != "Empty",
                           Condition != "Norm") %>% 
      dplyr::mutate(Condition = factor(Condition,
                                       levels = Conditions)) %>%
      dplyr::mutate(Sample_ID = paste0(Condition,"_",BioReplicate)) %>% 
      dplyr::select(Sample_ID, BioReplicate, Condition) %>%
      dplyr::arrange(BioReplicate) %>%
      dplyr::distinct()


cond <- factor(annottmt1$Condition)
biorep <- factor(annottmt1$BioReplicate)

design <- model.matrix(~biorep+cond)


tmt1proc1 <- tmt_procdat %>%
            dplyr::mutate(Condition = factor(Condition,
                                             levels = Conditions),
                          Protein = as.character(Protein),
                          Sample_ID = paste0(Condition,"_",BioReplicate)) %>%
            dplyr::filter(str_detect(Protein, "^Biognosys_pep", negate = TRUE)) %>%
            dplyr::select(Sample_ID, Protein, Abundance)

rsh_tmtproc1 <- tidyr::pivot_wider(tmt1proc1,
                                   names_from = Sample_ID,
                                   values_from = Abundance,
                                   values_fn = list(Abundance = median)) %>%
            dplyr::select(ID = Protein, annottmt1$Sample_ID)


## Test with TMT data (non-paired design) ----

Conditions = c("prim", "rec")

annottmt1 <- dplyr::filter(tmt_annot,
                           Condition != "Empty",
                           Condition != "Norm") %>% 
      dplyr::mutate(Condition = factor(Condition,
                                       levels = Conditions)) %>%
      dplyr::mutate(Sample_ID = paste0(Condition,"_",BioReplicate)) %>% 
      dplyr::select(Sample_ID, BioReplicate, Condition) %>%
      dplyr::arrange(Condition) %>%
      dplyr::distinct()


cond <- factor(annottmt1$Condition)

design <- model.matrix(~cond)


tmt1proc1 <- tmt_procdat %>%
      dplyr::mutate(Condition = factor(Condition,
                                       levels = Conditions),
                    Protein = as.character(Protein),
                    Sample_ID = paste0(Condition,"_",BioReplicate)) %>%
      dplyr::filter(str_detect(Protein, "^Biognosys_pep", negate = TRUE)) %>%
      dplyr::select(Sample_ID, Protein, Abundance)

rsh_tmtproc1 <- tidyr::pivot_wider(tmt1proc1,
                                   names_from = Sample_ID,
                                   values_from = Abundance,
                                   values_fn = list(Abundance = median)) %>%
      dplyr::select(ID = Protein, annottmt1$Sample_ID)



