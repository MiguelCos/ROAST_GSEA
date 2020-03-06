### roastKEGG function ----  


roastKEGG <- function(data,
                      geneIDtype = "SYMBOL",
                      orgDB = "org.Hs.eg.db",
                      organism = "hsa", # here the sintax should correspond with the KEGG sintax
                      design = designMatrix,
                      n_rotations = 999,
                      minSetSize = 1,
                      maxSetSize = 1000,
                      pvalueCutoff = 0.05) {
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
      
      path2entrez <- KEGGREST::keggLink("pathway", organism)
      
      entrezcor <- str_remove_all(names(path2entrez), "^hsa:") %>% 
         str_trim()
      
      path2entrez <- str_remove_all(path2entrez, "^path:") %>% 
         str_trim()
      
      names(path2entrez) <- entrezcor
      
      list_path2entrez <- split(x = names(path2entrez),
                                f = path2entrez)
      
      index <- limma::ids2indices(gene.sets = list_path2entrez,
                                  identifiers = genesindata,
                                  remove.empty = TRUE)
      
      leindex <- sapply(index, length)
      
      sublogi1 <- between(leindex, minSetSize, maxSetSize) 
      index2 <- index[sublogi1] 
      
      ## Run roast ----
      
      roast_out <- roast(y = matrix1,
                         contrast= ncol(designMatrix),
                         design = designMatrix, 
                         nrot = n_rotations, 
                         index = index2)
      
      
      return(roast_outfil)
      
      
}
