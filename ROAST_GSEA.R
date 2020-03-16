#### ROAST GSEA data analysis script v 0.1 ----
### Miguel Cosenza 18.02.2020 

### Load required packages ----

library(limma)
library(reactome.db)
library(qdapTools)

### Load dataset ----

inputdata <- read.delim("input_limma_example.txt",
                        stringsAsFactors = FALSE) %>% na.omit() %>%
            dplyr::rename(UNIPROT = Name)

load(file = here::here("CPTAC_CCRCC_Proteome_tabular_reduced.Rda"))

inputdata <- selected_patdata %>% mutate(Name = row.names(data_protsel)) %>%
            na.omit() %>% dplyr::rename(SYMBOL = Name)
            #dplyr::select(Name, 2:dim(data_protsel)[2]-1) %>% 
            #dplyr::rename(SYMBOL = Name) %>% na.omit()

### Prepare data to ROAST ----

### Prep matrix ----

genenames <- inputdata$SYMBOL

mappednames <- clusterProfiler::bitr(geneID = genenames,
                                     fromType = "SYMBOL",
                                     toType = c("ENTREZID", "UNIPROT",
                                                "GENENAME"),
                                     OrgDb = org.Hs.eg.db) %>% na.omit()

matnewname <- left_join(inputdata,
                        mappednames,
                        by = "SYMBOL") 

matnewname2 <- dplyr::select(matnewname,
                             -c(SYMBOL, UNIPROT, GENENAME)) %>%  
            dplyr::rename(Name = ENTREZID) %>% 
            dplyr::filter(is.na(Name) == FALSE) %>% 
            dplyr::distinct()

genesindata <- matnewname2$Name

matnewname3 <- dplyr::select(matnewname2,
                             -Name)

newmat <- as.matrix(matnewname3)

newmat[1:5,1:5]

dim(newmat)

### Prep index ----

library(reactome.db)

reactlistentrez <- as.list(reactomePATHID2EXTID)

head(reactlistentrez)[1:2]

index2 <- ids2indices(gene.sets = reactlistentrez,
            identifiers = genesindata,
            remove.empty = TRUE) 

thematching <- function(x){
      intersect(x, genesindata)
}

getfromgeneset <- lapply(reactlistentrez,thematching)

islengthzero <- sapply(getfromgeneset,
                       function(x){length(x)>0})

getfromgeneset2 <- getfromgeneset[islengthzero]

index <- getfromgeneset2

idsgene <- genesindata

"R-HSA-525793"

index[names(index) == "R-HSA-525793"]

load(file = here::here("exp_design_CPTACCRCC_Proteomics_reduced.Rda"))

## Define experimental design ####

condition1 <- 5 # number of samples associated to the first condition (treatment,  stage, patient, etc...)
condition2 <- 7 # number of samples associated to the second condition 

experiment <- c(rep(2,condition1),rep(1,condition2))
design <- model.matrix(~experiment)

n.rot <- 9999
roastres1 <- roast(newmat,contrast=ncol(designCPTACCRCCProteomics),
                   design=designCPTACCRCCProteomics_reduced, 
                   nrot = n.rot, index=index2) ## The comparison of interest (Condition1 vs Condition2)

roastres1 = roastres1[order(rownames(roastres1)),]




#### Try function roastReactome -----

## Load required packages
require(organism, character.only = TRUE) || stop(paste("package", organism, "is required", sep=" "))
require(limma) || stop("Package limma is required")
require(reactome.db) || stop("Package reactome.db is required")
require(clusterProfiler) || stop("Package clusterProfiler is required")
require(dplyr)

## Generate matrix to roast ----  

## Data table 
input_roast <- read.delim(file = "roast_sample_input_cptac_ccrcc_reduced.txt")

## Data design 
load(file = "sample_idesign_cptac_ccrcc_reduced.Rda")


## Define 'function' parameters ----

data <- input_roast

geneIDtype <- "SYMBOL"

organism <-"org.Hs.eg.db"

designMatrix <- sample_idesign_cptac_ccrcc_reduced

minSetSize <- 50
maxSetSize <- 150

n_rotations = 9999

## Run function ----

genenames <- data$ID

mapentrez <- clusterProfiler::bitr(geneID = genenames,
                                   fromType = geneIDtype,
                                   toType = "ENTREZID",
                                   OrgDb = org.Hs.eg.db) %>% na.omit()

names(data) <- c(geneIDtype, names(data)[2:eval(dim(data)[2])])

premat1 <- dplyr::left_join(data,
                            mapentrez,
                            by = geneIDtype) %>% na.omit()

premat2 <- dplyr::select(premat1,
                         -all_of(geneIDtype)) %>%  
      dplyr::rename(Name = ENTREZID) %>% 
      dplyr::filter(is.na(Name) == FALSE) %>% 
      dplyr::distinct()

genesindata <- premat2$Name

premat3 <- dplyr::select(premat2,
                         -Name)  

matrix1 <- as.matrix(premat3)

## Prep index for roast ----

reactlistentrez <- as.list(reactomePATHID2EXTID)

index <- limma::ids2indices(gene.sets = reactlistentrez,
                            identifiers = genesindata,
                            remove.empty = TRUE) 

lenindex <- sapply(index, length)

sublogi1 <- between(lenindex, minSetSize, maxSetSize) 
index2 <- index[sublogi1]  

lenindex2 <- sapply(index2, length) # test temp

## Run roast ----

roast_out <- roast(y = matrix1,
                   contrast= ncol(designMatrix),
                   design = designMatrix, 
                   nrot = n_rotations, 
                   index = index2)


### Develop idea of roastGO function -----

input_roast <- read.delim(file = "roast_sample_input_cptac_ccrcc_reduced.txt",
                         stringsAsFactors = FALSE)

load(file = "exp_design_CPTACCRCC_Proteomics_reduced.Rda")

## Define 'function' parameters ----

data <- input_roast
ontology = "MF"
geneIDtype <- "SYMBOL"
organism <-"org.Hs.eg.db"
designMatrix <- designCPTACCRCCProteomics_reduced

minSetSize <- 50
maxSetSize <- 150

n_rotations = 999
pval_threshold = 0.05

## Run function ----

premat2 <- dplyr::filter(data, is.na(ID) == FALSE) %>% 
           dplyr::distinct() %>% na.omit()

genesindata <- premat2$ID

premat3 <- dplyr::select(premat2,
                         -ID) 

matrix1 <- as.matrix(premat3)

## Prep index for roast ----

goterm_n_ontol <- AnnotationDbi::Ontology(GO.db::GOTERM)

gos_to_test <- goterm_n_ontol[goterm_n_ontol == ontology]  

goterm_n_id <- AnnotationDbi::mapIds(GO.db,
                                     keys = names(gos_to_test),
                                     column = "TERM",
                                     keytype = "GOID")

goterm_n_iddf <- data.frame(GOID = names(goterm_n_id),
                            GOTERM = goterm_n_id)


golist <- suppressMessages(
      AnnotationDbi::mapIds(org.Hs.eg.db, keys=names(gos_to_test), 
                                column=geneIDtype,
                                keytype="GOALL", multiVals='list'))

index <- limma::ids2indices(gene.sets = golist,
                            identifiers = genesindata,
                            remove.empty = TRUE) 

lenindex <- sapply(index, length)

sublogi1 <- between(lenindex, minSetSize, maxSetSize) 
index2 <- index[sublogi1]  

genesinterm <- qdapTools::list2df(index2,
                                  col1 = "ENTREZID",
                                  col2 = "GOID") %>%
            dplyr::mutate(ENTREZID = as.character(ENTREZID))

symb1 <- clusterProfiler::bitr(genesinterm$ENTREZID,
                               fromType = "ENTREZID",
                               toType = "SYMBOL",
                               OrgDb = org.Hs.eg.db,
                               drop = FALSE)

genesintermread <- left_join(genesinterm, symb1,
                             by = "ENTREZID") %>% 
                  left_join(., goterm_n_iddf,
                            by = "GOID")

roast_out <- roast(y = matrix1,
                   contrast= ncol(designMatrix),
                   design = designMatrix, 
                   nrot = n_rotations, 
                   index = index2)

roast_out2 <- dplyr::mutate(roast_out,
                           GOID = row.names(roast_out)) %>% 
            dplyr::left_join(.,goterm_n_iddf,
                             by = "GOID")

roast_out2 <- dplyr::filter(roast_out2,
                           FDR <= pval_threshold)


roastResult <- list(roastOutput = roast_out2,
                    GenesPerTerm = genesintermread)

return(roastResult)


### Develop roastMSigDB function ----

## Define 'function' parameters ----

data <- input_roast
geneIDtype <- "gene_symbol" # this can be "gene_symbol" or "entrez_gene"
organism <- "Homo sapiens" # this can be any resulting from calling msigdbr::msigdbr_show_species()
category <- "C2" # Any of the main categories presented here: https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp
subcategory <- NULL
specific_category <- "NABA"
designMatrix <- designCPTACCRCCProteomics_reduced

minSetSize <- 1
maxSetSize <- 200

n_rotations = 999

## Run function ----

### Prep matrix ----

premat2 <- dplyr::filter(data, is.na(ID) == FALSE) %>% 
      dplyr::distinct() %>% na.omit()

genesindata <- premat2$ID

premat3 <- dplyr::select(premat2,
                         -ID) 

matrix1 <- as.matrix(premat3)

### Prep index for roast -----

require(msigdbr)

msigtab <- msigdbr::msigdbr(species = organism,
                            category = category,
                            subcategory = subcategory)

if (isEmpty(specific_category) == FALSE){
   
   if (specific_category == "NABA"){
      msigtab <- dplyr::filter(msigtab,
                               str_detect(gs_name, "^NABA_"))
   }
}   
   list_msig <- msigtab %>% split(x = .$gene_symbol, f = .$gs_id)
   
   index <- limma::ids2indices(gene.sets = list_msig,
                               identifiers = genesindata,
                               remove.empty = TRUE)
   
   lenindex <- sapply(index, length)
   
   sublogi1 <- between(lenindex, minSetSize, maxSetSize) 
   index2 <- index[sublogi1]
   
## Run roast ----

roast_out <- roast(y = matrix1,
                   contrast= ncol(designMatrix),
                   design = designMatrix, 
                   nrot = n_rotations, 
                   index = index2)
   
## Develop roastKEGG function -----
   
## Define 'function' parameters ----

data =  input_roast
geneIDtype = "SYMBOL"
orgDB = "org.Hs.eg.db"
organism = "hsa" # here the sintax should correspond with the KEGG sintax
design = sample_idesign_cptac_ccrcc_reduced
n_rotations = 999
minSetSize = 1
maxSetSize = 1000
pvalueCutoff = 0.05
   
### Function ----

require(KEGGREST) || stop("Package KEGGREST is required")
require(orgDB, character.only = TRUE) || stop(paste("package", orgDb, "is required", sep=" "))
require(limma) || stop("Package limma is required")
suppressMessages(require(clusterProfiler)) || stop("Package clusterProfiler is required")
require(dplyr) || stop("Package dplyr is required")
require(stringr) || stop("Package stringr is required")

## Generate matrix to roast ----

if (geneIDtype != "ENTREZID"){
   
   genenames <- data$ID
   
   mapentrez <- clusterProfiler::bitr(geneID = genenames,
                                      fromType = geneIDtype,
                                      toType = "ENTREZID",
                                      OrgDb = org.Hs.eg.db) %>% na.omit()
   
   names(data) <- c(geneIDtype, names(data)[2:eval(dim(data)[2])])
   
   premat1 <- dplyr::left_join(data,
                               mapentrez,
                               by = geneIDtype)
   
   premat2 <- dplyr::select(premat1,
                            -geneIDtype) %>%  
      dplyr::rename(Name = ENTREZID) %>% 
      dplyr::filter(is.na(Name) == FALSE) %>% 
      dplyr::distinct() %>% na.omit()
   
   genesindata <- premat2$Name
   
   premat3 <- dplyr::select(premat2,
                            -Name) 
   
   matrix1 <- as.matrix(premat3)
} else {
   
   premat2 <- dplyr::filter(data, is.na(ID) == FALSE) %>% 
      dplyr::distinct() %>% na.omit()
   
   genesindata <- premat2$ID
   
   premat3 <- dplyr::select(premat2,
                            -ID) 
   
   matrix1 <- as.matrix(premat3)
   
}

## Prep index for roast ----

path2entrez <- KEGGREST::keggLink("pathway", organism)

entrezcor <- str_remove_all(names(path2entrez), "^hsa:") %>% 
             str_trim()

path2entrez <- str_remove_all(path2entrez, "^path:") %>% 
               str_trim()

names(path2entrez) <- entrezcor

list_path2entrez <- split(x = names(path2entrez),
                          f = path2entrez)

keggid_to_term <- keggList("pathway","hsa")

keggid <- str_remove_all(names(keggid_to_term), "^path:") %>%
                   str_trim()

keggterm <- str_remove_all(keggid_to_term, " \\- Homo sapiens \\(human\\)") %>%
            str_trim()

keggidtoterm_df <- data.frame(KEGGID = keggid,
                              KEGGTERM = keggterm)


index <- limma::ids2indices(gene.sets = list_path2entrez,
                            identifiers = genesindata,
                            remove.empty = TRUE)

leindex <- sapply(index, length)

sublogi1 <- between(leindex, minSetSize, maxSetSize) 
index2 <- index[sublogi1] 


genesinterm <- qdapTools::list2df(index2,
                                  col1 = "ENTREZID",
                                  col2 = "KEGGID") %>%
      dplyr::mutate(ENTREZID = as.character(ENTREZID))

suppressWarnings(
      symb1 <- clusterProfiler::bitr(genesinterm$ENTREZID,
                                     fromType = "ENTREZID",
                                     toType = "SYMBOL",
                                     OrgDb = orgDB,
                                     drop = FALSE)
)

suppressWarnings(
      genesintermread <- left_join(genesinterm, symb1,
                                   by = "ENTREZID") %>% 
            left_join(., keggidtoterm_df,
                      by = "KEGGID")
)

## Run roast ----

roast_out <- roast(y = matrix1,
                   contrast= ncol(designMatrix),
                   design = designMatrix, 
                   nrot = n_rotations, 
                   index = index2)

roast_out2 <- dplyr::mutate(roast_out,
                            KEGGID = row.names(roast_out)) %>% 
   dplyr::left_join(.,keggidtoterm_df,
                    by = "KEGGID")

roast_out2 <- dplyr::filter(roast_out2,
                            FDR <= pval_threshold)

genesintermread <- dplyr::filter(genesintermread,
                                 KEGGID %in% roast_out2$KEGGID)

roastResult <- list(roastOutput = roast_out2,
                    GenesPerTerm = genesintermread)

return(roastResult)

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
