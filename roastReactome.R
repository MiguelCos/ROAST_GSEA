### ROAST - Reactome function ----

roastReactome <- function(data, geneIDtype = "SYMBOL", 
                          organism = "org.Hs.eg.db", 
                          design = designMatrix,
                          n_rotations = 9999,
                          minSetSize = 1,
                          maxSetSize = 1506){
      
      ## Load required packages
      require(organism, character.only = TRUE) || stop(paste("package", organism, "is required", sep=" "))
      require(limma) || stop("Package limma is required")
      require(reactome.db) || stop("Package reactome.db is required")
      require(clusterProfiler) || stop("Package clusterProfiler is required")
      require(dplyr)
   
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
   
   reactlistentrez <- as.list(reactomePATHID2EXTID)
   
   index <- limma::ids2indices(gene.sets = reactlistentrez,
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
