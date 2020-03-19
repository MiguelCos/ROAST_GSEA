## Temp script for testing functions ----

## Data table 
input_roast <- read.delim(file = "roast_sample_input_cptac_ccrcc_reduced.txt",
                          stringsAsFactors = FALSE)

## Data design 
load(file = "sample_idesign_cptac_ccrcc_reduced.Rda")

sample_idesign_cptac_ccrcc_reduced

## Test function roastReactome ----

test_roastReactome <- roastReactome(data = input_roast,
                                    geneIDtype = "SYMBOL", 
                                    organism = "org.Hs.eg.db", 
                                    design = sample_idesign_cptac_ccrcc_reduced,
                                    n_rotations = 999,
                                    minSetSize = 50,
                                    maxSetSize = 200)

## Test function roastGO ----
source("roastGO.R")
source("simplifyGO.R")
test_roastGO <- roastGO(data = input_roast,
                        geneIDtype = "SYMBOL",
                        ontology = "MF",
                        organism = "org.Hs.eg.db", 
                        design = sample_idesign_cptac_ccrcc_reduced,
                        n_rotations = 999,
                        minSetSize = 100,
                        maxSetSize = 500,
                        pvalueCutoff = 0.01)

test_simplifyGO <- simplifyGO(test_roastGO, cutoff = 0.7, by = "PValue")

clprofout <- clusterProfiler::enrichGO(input_roast$ID,
                                       OrgDb = org.Hs.eg.db,
                                       keyType = "SYMBOL",
                                       ont = "MF")
      
## Test function propChangePlot ----

test_propChangeplot <- propChangePlot(test_roastGO)

## Test function ridgeplotRoast ----

test_ridgeplotRoast <- ridgeplotRoast(test_roastGO)

## Test function roastGO ----

test_roastMSigDB <- roastMSigDB(data = input_roast,
                            geneIDtype = "gene_symbol", 
                            organism = "Homo sapiens",
                            category = "C6", 
                            subcategory = NULL,
                            specific_category = NULL,
                            design = sample_idesign_cptac_ccrcc_reduced,
                            n_rotations = 999,
                            minSetSize = 1,
                            maxSetSize = 1000)


## Test function roastKEGG ----

test_roastKEGG <- roastKEGG(data = input_roast,
                            geneIDtype = "SYMBOL",
                            orgDB = "org.Hs.eg.db",
                            organism = "hsa", # here the syntax should correspond with the KEGG sintax
                            design = sample_idesign_cptac_ccrcc_reduced,
                            n_rotations = 999,
                            minSetSize = 150,
                            maxSetSize = 300,
                            pvalueCutoff = 0.05)
