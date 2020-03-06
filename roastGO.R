## roastGO function ----

roastGO <- function(data,
                    geneIDtype = "SYMBOL",
                    ontology = "MF",
                    organism = "org.Hs.eg.db", 
                    design = designMatrix,
                    n_rotations = 9999,
                    minSetSize = 1,
                    maxSetSize = 1506){
      
      ## Load required packages ----
      require(organism, character.only = TRUE) || stop(paste("package", organism, "is required", sep=" "))
      require(limma) || stop("Package limma is required")
      require(GO.db) || stop("Package reactome.db is required")
      require(AnnotationDbi)
      require(dplyr)
      
      ## Generate matrix for roast ----
      
      enenames <- data$ID
      
      premat2 <- dplyr::filter(data, is.na(ID) == FALSE) %>% 
                  dplyr::distinct() %>% na.omit()
      
      genesindata <- premat2$ID
      
      premat3 <- dplyr::select(premat2, -ID) 
      
      matrix1 <- as.matrix(premat3)
      
      ## Prep index for roast ----
      
      goterm_n_ontol <- AnnotationDbi::Ontology(GO.db::GOTERM)
      gos_to_test <- goterm_n_ontol[goterm_n_ontol == ontology]  
      
      golist <- AnnotationDbi::mapIds(org.Hs.eg.db, keys=names(gos_to_test), 
                                     column=geneIDtype,
                                     keytype="GOALL", multiVals='list')
      
      index <- limma::ids2indices(gene.sets = golist,
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
      
      return(roast_out)
      
}
