## roastGO function ----

roastGO <- function(data,
                    geneIDtype = "SYMBOL",
                    ontology = "MF",
                    organism = "org.Hs.eg.db", 
                    design,
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
   
      suppressWarnings(
      suppressMessages(
   goterm_n_ontol <- AnnotationDbi::Ontology(GO.db::GOTERM)
      ))
   
   gos_to_test <- goterm_n_ontol[goterm_n_ontol == ontology]  
   
      suppressWarnings(
      suppressMessages(
   goterm_n_id <- AnnotationDbi::mapIds(GO.db,
                                        keys = names(gos_to_test),
                                        column = "TERM",
                                        keytype = "GOID")
      ))
      
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
   suppressMessages(
   symb1 <- clusterProfiler::bitr(genesinterm$ENTREZID,
                                  fromType = "ENTREZID",
                                  toType = "SYMBOL",
                                  OrgDb = org.Hs.eg.db,
                                  drop = FALSE)
   ))
   
   suppressWarnings(
      suppressMessages(
   genesintermread <- left_join(genesinterm, symb1,
                                by = "ENTREZID") %>% 
      left_join(., goterm_n_iddf,
                by = "GOID")
      ))
   # Run ROAST ----
   
   roast_out <- roast(y = matrix1,
                      contrast= ncol(design),
                      design = design, 
                      nrot = n_rotations, 
                      index = index2)
   
   # Process ROAST output ----
   suppressWarnings(
      suppressMessages(
   roast_out2 <- dplyr::mutate(roast_out,
                               GOID = row.names(roast_out)) %>% 
      dplyr::left_join(.,goterm_n_iddf,
                       by = "GOID") %>% 
      dplyr::rename(CategoryID = GOID, CategoryTerm = GOTERM)
   
      ))
   
   roast_out2 <- dplyr::filter(roast_out2,
                               FDR <= pvalueCutoff)
   
   fdrnterm <- dplyr::select(roast_out2,
                             CategoryTerm, FDR, NGenes)
   suppressWarnings(
      suppressMessages(
   limma_out <- lmFit(object = matrix1,
                      design = design)
      ))
   
   limma_out2 <- eBayes(limma_out)
   
   suppressWarnings(
      suppressMessages(
   limma_tab <- topTable(limma_out2, number = dim(matrix1)[1])
      ))
   
   suppressWarnings(
      suppressMessages(
   log2FCs <- dplyr::mutate(limma_tab,
                            ENTREZID = row.names(limma_tab)) %>% 
      dplyr::filter(ENTREZID %in% symb1$ENTREZID) %>% 
      dplyr::select(log2FC = eval(dim(.)[2]-5),ENTREZID) %>% 
      dplyr::left_join(., genesintermread, by = "ENTREZID") %>% ## Look for a way to extract the z-score values as they are used by the ROAST algorithm
      dplyr::select(log2FC, ENTREZID, SYMBOL, CategoryID = GOID, CategoryTerm = GOTERM) %>% 
      dplyr::left_join(.,fdrnterm, by = "CategoryTerm") %>%
      dplyr::filter(is.na(FDR) == FALSE)
      ))
   
   suppressWarnings(
      suppressMessages(
   meanFCTerm <- log2FCs %>% 
      dplyr::group_by(CategoryID, CategoryTerm) %>% 
      dplyr::summarise(meanlog2FC = mean(log2FC)) %>% dplyr::ungroup() %>% 
      dplyr::filter(CategoryID %in% roast_out2$CategoryID) %>% 
      dplyr::select(-CategoryTerm, CategoryID)
      ))
   
   suppressWarnings(
      suppressMessages(
   roast_out3 <- dplyr::left_join(roast_out2,
                                  meanFCTerm,
                                  by = "CategoryID")
      ))
   
   
   roastResult <- list(roastOutput = roast_out3,
                       GenesPerTerm = genesintermread,
                       log2FCs = log2FCs,
                       ontology = ontology,
                       organism = organism)
   
   return(roastResult)
      
}
