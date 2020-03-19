### roastKEGG function ----  


roastKEGG <- function(data,
                      geneIDtype = "SYMBOL",
                      orgDB = "org.Hs.eg.db",
                      organism = "hsa", # here the sintax should correspond with the KEGG sintax
                      design = designMatrix,
                      n_rotations = 999,
                      minSetSize = 1,
                      maxSetSize = 1000,
                      pvalueCutoff = 0.05,
                      exclusionList = TRUE) {
      
   ## Load required packages ----
      
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
                                          KEGGID %in% keggid)
         
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
         
      }
      
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
      
      ## Process output ----
      
      genesintermread <- dplyr::filter(genesintermread,
                                       KEGGID %in% roast_out2$KEGGID)
      
      
      roastResult <- list(roastOutput = roast_out2,
                          GenesPerTerm = genesintermread)
      
      if(exclusionList == TRUE){
         message("KEGG pathways associated with the next category classes were excluded from the analysis:",
                 exclusion_class,"; ",paste(exclusion_subclass, sep = " ", collapse = "; "))
      }
      
      return(roastResult)
      
}
