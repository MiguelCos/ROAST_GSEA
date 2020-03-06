## roastMSigDB function ----

roastMSigDB <- function(data,
                        geneIDtype = "gene_symbol", # this can be "gene_symbol" or "entrez_gene"
                        organism = "Homo sapiens", # this can be any resulting from calling msigdbr::msigdbr_show_species()
                        category = "H", # Any of the main categories presented here: https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp
                        subcategory = NULL,
                        specific_category = NULL,
                        design = designMatrix,
                        n_rotations = 999,
                        minSetSize = 1,
                        maxSetSize = 200){
      
      ## Load required packages ----  
      require(limma) || stop("Package limma is required")
      require(dplyr) || stop("Package dplyr is required")
      require(msigdbr) || stop("Package msigdbr is required")
      require(stringr) || stop("Package stringr is required")
      
      ## Generate matrix for roast ----
      
      premat2 <- dplyr::filter(data, is.na(ID) == FALSE) %>% 
            dplyr::distinct() %>% na.omit()
      
      genesindata <- premat2$ID
      
      premat3 <- dplyr::select(premat2,
                               -ID) 
      
      matrix1 <- as.matrix(premat3)
      
      ## Prep index for roast ----
      
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
      
      return(roast_out)
      
}
                        