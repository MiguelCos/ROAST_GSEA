## roastGO function ----

roastGO <- function(data,
                    geneIDtype = "SYMBOL",
                    ontology = "MF",
                    organism = "org.Hs.eg.db", 
                    design = designMatrix,
                    n_rotations = 9999,
                    minSetSize = 1,
                    maxSetSize = 1506,
                    pvalueCutoff = 0.05){
      
      ## Load required packages ----
      require(organism, character.only = TRUE) || stop(paste("package", organism, "is required", sep = " "))
      require(limma) || stop("Package limma is required")
      require(GO.db) || stop("Package reactome.db is required")
      require(AnnotationDbi)
      require(dplyr)
      
      ## Generate matrix for roast ----
      
      premat2 <- dplyr::filter(data, is.na(ID) == FALSE) %>% 
            dplyr::distinct() %>% na.omit()
      
      genesindata <- premat2$ID
      
      premat3 <- dplyr::select(premat2,
                               -ID) 
      
      matrix1 <- as.matrix(premat3)
      
      ## Prep index for roast ----
      
      suppressMessages(
      goterm_n_ontol <- AnnotationDbi::Ontology(GO.db::GOTERM)
      )
      
      gos_to_test <- goterm_n_ontol[goterm_n_ontol == ontology]  
      
      suppressMessages(
            goterm_n_id <- AnnotationDbi::mapIds(GO.db,
                                                 keys = names(gos_to_test),
                                                 column = "TERM",
                                                 keytype = "GOID")
      )
      
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
      
      suppressWarnings(
      symb1 <- clusterProfiler::bitr(genesinterm$ENTREZID,
                                     fromType = "ENTREZID",
                                     toType = "SYMBOL",
                                     OrgDb = org.Hs.eg.db,
                                     drop = FALSE)
      )
      
      suppressWarnings(
      genesintermread <- left_join(genesinterm, symb1,
                                   by = "ENTREZID") %>% 
         left_join(., goterm_n_iddf,
                   by = "GOID")
      )
      
      ## Run roast ----  
      
      roast_out <- roast(y = matrix1,
                         contrast= ncol(design),
                         design = design, 
                         nrot = n_rotations, 
                         index = index2)
      
      suppressWarnings(
      roast_out2 <- dplyr::mutate(roast_out,
                                  GOID = row.names(roast_out)) %>% 
            dplyr::left_join(.,goterm_n_iddf,
                             by = "GOID")
      )
      
      roast_out2 <- dplyr::filter(roast_out2,
                                  FDR <= pval_threshold)
      
      roastResult <- list(roastOutput = roast_out2,
                          GenesPerTerm = genesintermread)
      
      return(roastResult)
      
      
}
