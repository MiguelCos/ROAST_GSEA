### ROAST - Reactome function ----

roastReactome <- function(data, geneIDtype = "SYMBOL", 
                          organism = "org.Hs.eg.db", 
                          design,
                          n_rotations = 9999,
                          minSetSize = 1,
                          maxSetSize = 1506,
                          pvalueCutoff = 0.05,
                          exclusionList = TRUE){
      
      ## Load required packages
      require(organism, character.only = TRUE) || stop(paste("package", organism, "is required", sep=" "))
      require(limma) || stop("Package limma is required")
      require(reactome.db) || stop("Package reactome.db is required")
      require(clusterProfiler) || stop("Package clusterProfiler is required")
      require(dplyr)
      require(stringr)
   
   ## Generate matrix to roast ----
   
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
   
   ## Prep index for roast ----
   
   suppressMessages(suppressWarnings(
      toexclude <- AnnotationDbi::select(reactome.db,
                                         keys = keys(reactome.db, keytype = "PATHID"),
                                         keytype =  "PATHID",
                                         columns = c("PATHID", "PATHNAME"))
   ))
   
   if(exclusionList == TRUE){
      
      
      exluded1 <-  toexclude %>% 
         dplyr::filter(str_detect(PATHNAME, "Homo sapiens: ")) %>%
         dplyr::filter(str_detect(PATHNAME, "disease"))
      
      exluded2 <- toexclude %>% 
         dplyr::filter(str_detect(PATHNAME, "Homo sapiens: ")) %>%
         dplyr::filter(str_detect(PATHNAME, "Disease"))
      
      exluded3 <- toexclude %>% 
         dplyr::filter(str_detect(PATHNAME, "Homo sapiens: ")) %>%
         dplyr::filter(str_detect(PATHNAME, "disorder"))
      
      exluded4 <- toexclude %>% 
         dplyr::filter(str_detect(PATHNAME, "Homo sapiens: ")) %>%
         dplyr::filter(str_detect(PATHNAME, "infection"))
      
      exluded5 <- toexclude %>% 
         dplyr::filter(str_detect(PATHNAME, "Homo sapiens: ")) %>%
         dplyr::filter(str_detect(PATHNAME, "Infection"))
      
      
      reactome_excluded <- bind_rows(exluded1,
                                     exluded2,
                                     exluded3,
                                     exluded4,
                                     exluded5) 
      
      pathexcl <- reactome_excluded$PATHID
      
      reactlistentrez1 <- as.list(reactomePATHID2EXTID)
      
      reactlistentrez <- reactlistentrez1[which(!names(reactlistentrez1) %in% pathexcl)]
      
      index <- limma::ids2indices(gene.sets = reactlistentrez,
                                  identifiers = genesindata,
                                  remove.empty = TRUE)
      
      lenindex <- sapply(index, length)
      
      sublogi1 <- between(lenindex, minSetSize, maxSetSize) 
      index2 <- index[sublogi1]  
      
      pathterm_n_iddf <- toexclude %>% dplyr::filter(str_detect(PATHNAME, species)) %>% 
                           dplyr::mutate(PATHNAME = str_remove(PATHNAME, ".*: "))
      
      
   } else {
   
      reactlistentrez <- as.list(reactomePATHID2EXTID)
   
      index <- limma::ids2indices(gene.sets = reactlistentrez,
                                identifiers = genesindata,
                                remove.empty = TRUE)
   
      lenindex <- sapply(index, length)
   
      sublogi1 <- between(lenindex, minSetSize, maxSetSize) 
      index2 <- index[sublogi1]  
      
      pathterm_n_iddf <- toexclude %>% dplyr::filter(str_detect(PATHNAME, species)) %>% 
                           dplyr::mutate(PATHNAME = str_remove(PATHNAME, ".*: "))
   
   }
   
   genesinterm <- qdapTools::list2df(index2,
                                     col1 = "ENTREZID",
                                     col2 = "PATHID") %>%
                  dplyr::mutate(ENTREZID = as.character(ENTREZID))
   
   
   suppressWarnings(suppressMessages(
      
      symb1 <- clusterProfiler::bitr(genesinterm$ENTREZID,
                                        fromType = "ENTREZID",
                                        toType = "SYMBOL",
                                        OrgDb = eval(as.name(organism)),
                                        drop = FALSE)
      ))
   
   suppressWarnings(suppressMessages(
         genesintermread <- left_join(genesinterm, symb1,
                                      by = "ENTREZID") %>% 
            left_join(., pathterm_n_iddf,
                      by = "PATHID")
      ))
   # Run roast ----
   
   roast_out <- fry(y = matrix1,
                      contrast= ncol(design),
                      design = design, 
                      nrot = n_rotations, 
                      index = index2)
   
   # Process ROAST output ----
   suppressWarnings(suppressMessages(
      
         roast_out2 <- dplyr::mutate(roast_out,
                                     PATHID = row.names(roast_out)) %>% 
            dplyr::left_join(.,pathterm_n_iddf,
                             by = "PATHID") %>% 
            dplyr::rename(CategoryID = PATHID, CategoryTerm = PATHNAME)
         
      ))
   
   roast_out2 <- dplyr::filter(roast_out2,
                               FDR <= pvalueCutoff)
   
   fdrnterm <- dplyr::select(roast_out2,
                             CategoryTerm, FDR, NGenes)
   
   genesintermread <- dplyr::filter(genesintermread,
                                    PATHID %in% roast_out2$CategoryID)
   
   # Run limma to get log2FC values per protein and category ----
   
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
   
   # Get log2FC information from Limma and reformat output ----
   
   suppressWarnings(
      suppressMessages(
         log2FCs <- dplyr::mutate(limma_tab,
                                  ENTREZID = row.names(limma_tab)) %>% 
            dplyr::filter(ENTREZID %in% symb1$ENTREZID) %>% 
            dplyr::select(log2FC = eval(dim(.)[2]-5),ENTREZID) %>% 
            dplyr::left_join(., genesintermread, by = "ENTREZID") %>% ## Look for a way to extract the z-score values as they are used by the ROAST algorithm
            dplyr::select(log2FC, ENTREZID, SYMBOL, CategoryID = PATHID, CategoryTerm = PATHNAME) %>% 
            dplyr::left_join(.,fdrnterm, by = "CategoryTerm") %>%
            dplyr::filter(is.na(FDR) == FALSE)
      ))
   
   suppressWarnings(suppressMessages(
      
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
                       log2FCs = log2FCs)
   
   if(exclusionList == TRUE){
      message("Pathways associated with human diseases, disorders and infections where excluded from de analysis")
      
      roastResult$ExcludedPathways <- reactome_excluded
   }
   
   return(roastResult)
}
