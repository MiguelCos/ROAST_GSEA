#### ROAST GSEA data analysis script v 0.1 ----
### Miguel Cosenza 18.02.2020 

### Load required packages ----

library(limma)
library(reactome.db)
library(qdapTools)
library(tidyverse)

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
pvalueCutoff = 0.01

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
                             by = "GOID") %>% 
            dplyr::rename(CategoryID = GOID, CategoryTerm = GOTERM)

roast_out2 <- dplyr::filter(roast_out2,
                           FDR <= pvalueCutoff)

fdrnterm <- dplyr::select(roast_out2,
                          CategoryTerm, FDR, NGenes)

limma_out <- lmFit(object = matrix1,
                   design = designMatrix)

limma_out2 <- eBayes(limma_out)

limma_tab <- topTable(limma_out2, number = dim(matrix1)[1])

log2FCs <- dplyr::mutate(limma_tab,
                         ENTREZID = row.names(limma_tab)) %>% 
   dplyr::filter(ENTREZID %in% symb1$ENTREZID) %>% 
   dplyr::select(log2FC = eval(dim(.)[2]-5),ENTREZID) %>% 
   dplyr::left_join(., genesintermread, by = "ENTREZID") %>% ## Look for a way to extract the z-score values as they are used by the ROAST algorithm
   dplyr::select(log2FC, ENTREZID, SYMBOL, CategoryID = GOID, CategoryTerm = GOTERM) %>% 
   dplyr::left_join(.,fdrnterm, by = "CategoryTerm") %>%
   dplyr::filter(is.na(FDR) == FALSE)

meanFCTerm <- log2FCs %>% 
   dplyr::group_by(CategoryID, CategoryTerm) %>% 
   dplyr::summarise(meanlog2FC = mean(log2FC)) %>% dplyr::ungroup() %>% 
   dplyr::filter(CategoryID %in% roast_out2$CategoryID) %>% 
   dplyr::select(-CategoryTerm, CategoryID)


roast_out3 <- dplyr::left_join(roast_out2,
                               meanFCTerm,
                               by = "CategoryID")


roastResult <- list(roastOutput = roast_out3,
                    GenesPerTerm = genesintermread,
                    log2FCs = log2FCs)

return(roastResult)

## Develop propChangeplot ----

test_roastGO

roastOutput = test_roastGO
show_n_terms = 25

# propChangeplot: Process the data ----

roastOutput <- roastOutput$roastOutput

toproplot <- dplyr::select(roastOutput,
                           NGenes, PropDown, PropUp, Direction, CategoryTerm,
                           FDR) %>%
               dplyr::top_n(n = show_n_terms,
                            wt = NGenes) %>%
            dplyr::mutate(DiffProp = abs(PropUp - PropDown),
                          PropDown = -PropDown) %>%
            tidyr::pivot_longer(cols = c(PropDown, PropUp),
                                names_to = "PropDirection",
                                values_to = "Proportion") %>%
            dplyr::group_by(CategoryTerm, NGenes) 
            

proplot <- ggplot(data = toproplot,
                  aes(x=fct_reorder(CategoryTerm, Proportion), y=Proportion))+
         coord_flip() +
         geom_line(aes(group = CategoryTerm))+
         geom_hline(yintercept = 0, color = "red")+
         #and the oints
         geom_point(aes(color=FDR, size = NGenes))+
         geom_text(data = dplyr::filter(toproplot, PropDirection == "PropUp"), 
             aes(label=round(Proportion,2)),
             hjust = -0.85)+
         geom_text(data = dplyr::filter(toproplot, PropDirection == "PropDown"), 
                   aes(label=round(Proportion,2)),
                   hjust = +1.85)+
         scale_y_continuous(expand=c(0.2,0), limits=c(-1, 1))+
         scale_x_discrete(labels = function(x)  str_wrap(x, width = 40))+
         labs(title = "Proportion of Up- or Down-regulated Proteins by Category", 
              subtitle= "Positive values = Proportion of up-regulated proteins in category",
               x="Biological category", y = "Proportion of Proteins")+
         theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 10),
               panel.background = element_blank(),
               panel.grid.major = element_blank(),
               panel.border = element_rect(colour = "black", fill=NA, size=1.5),
               axis.title=element_text(size=10,face="bold"),
         legend.justification = c(0, 1))

### Develop ridgeplotRoast function ----

test_roastGO

roastOutput = test_roastGO
show_n_terms = 25

# ridgeplotRoast: Process the data ----

require(ggridges)
require(ggplot2)

tofil <- roastOutput$roastOutput
toridge <- roastOutput$log2FCs

toproplot <- dplyr::select(roastOutput,
                           NGenes, PropDown, PropUp, Direction, CategoryTerm,
                           FDR) %>%
   dplyr::top_n(n = show_n_terms,
                wt = NGenes) %>%
   dplyr::mutate(DiffProp = abs(PropUp - PropDown),
                 PropDown = -PropDown) %>%
   tidyr::pivot_longer(cols = c(PropDown, PropUp),
                       names_to = "PropDirection",
                       values_to = "Proportion") %>%
   dplyr::group_by(CategoryTerm, NGenes) 

datatab <- dplyr::filter(toridge,
                         CategoryTerm %in% unique(toproplot$CategoryTerm)) %>% 
            dplyr::arrange(-NGenes)

ridges <- ggplot(data = datatab, aes(x = log2FC, y = CategoryTerm, fill = FDR))+
                  scale_fill_gradient(low = "#477af8", high = "#ff3333", name = "FDR")+
                  geom_density_ridges()+
                  xlab("Log2(Fold-change)")+ 
                  ylab("Biological Category")+
                  theme(axis.text.x = element_text(angle = 0, hjust = 0.5,size = 10, color = "black"),
                        panel.background = element_blank(),
                        panel.grid.major = element_blank(),
                        panel.border = element_rect(colour = "black", fill=NA, size=1.5),
                        axis.title=element_text(size=10,face="bold"))

return(ridges)


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
exclusionList = TRUE 
   
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

keggid_to_term <- keggList("pathway", organism)

keggid <- str_remove_all(names(keggid_to_term), "^path:") %>%
   str_trim()

keggterm <- str_remove_all(keggid_to_term, "\\ \\-.*") %>%
   str_trim()

keggidtoterm_df <- data.frame(KEGGID = keggid,
                              KEGGTERM = keggterm)

if (exclusionList == TRUE){
   getClass <- function(pathway){
      qery <- KEGGREST::keggGet(pathway)
      
      class <- qery[[1]]$CLASS
      
      return(class)
   }
   
   pathclss <- sapply(keggid, getClass)
   
   pathclss <- qdapTools::list2df(pathclss,
                                  col1 = "Class",
                                  col2 = "PathwayID") %>% 
      tidyr::separate(col = Class, into = c("Class", "Subclass"), sep = "; ")
   
   
   exclusion_subclass <- c("Drug resistance: antineoplastic",
                           "Endocrine system", "Aging",
                           "Circulatory system", "Xenobiotics biodegradation and metabolism",
                           "Drug resistance: antineoplastic", "Nervous system", "Sensory system",
                           "Excretory system", "Digestive system")
   
   exclusion_class <- c("Human Diseases")
   
   aftrexclud <- dplyr::filter(pathclss,
                               !Subclass %in% exclusion_subclass,
                               !Class %in% exclusion_class)
   
   keggid <- subset(keggid, keggid %in% aftrexclud$PathwayID)
   
   keggidtoterm_df <- dplyr::filter(keggidtoterm_df,
                                    KEGGID %in% keggidfil)
   
   path2entrez <- KEGGREST::keggLink("pathway", organism)
   
   entrezcor <- str_remove_all(names(path2entrez), 
                               paste0("^",organism,":")) %>% 
      str_trim()
   
   path2entrez <- str_remove_all(path2entrez, "^path:") %>% 
      str_trim()
   
   names(path2entrez) <- entrezcor
   
   path2entrez <- subset(path2entrez, path2entrez %in% keggidtoterm_df$KEGGID)
   
   list_path2entrez <- split(x = names(path2entrez),
                             f = path2entrez)
   
} else {

path2entrez <- KEGGREST::keggLink("pathway", organism)

entrezcor <- str_remove_all(names(path2entrez), 
                            paste0("^",organism,":")) %>% 
   str_trim()

path2entrez <- str_remove_all(path2entrez, "^path:") %>% 
   str_trim()

names(path2entrez) <- entrezcor

path2entrez <- subset(path2entrez, path2entrez %in% keggidtoterm_df$KEGGID)

list_path2entrez <- split(x = names(path2entrez),
                          f = path2entrez)


index <- limma::ids2indices(gene.sets = list_path2entrez,
                            identifiers = genesindata,
                            remove.empty = TRUE)

leindex <- sapply(index, length)

sublogi1 <- between(leindex, minSetSize, maxSetSize) 
index2 <- index[sublogi1] 
}

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
