## roastMSigDB function ----
data = input_roast
geneIDtype = "SYMBOL" # this can be "gene_symbol" or "entrez_gene"
organism = "Homo sapiens" # this can be any resulting from calling msigdbr::msigdbr_show_species()
orgDB = "org.Hs.eg.db"
category = "H" # Any of the main categories presented here: https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp
subcategory = NULL # Any of the subcategies presented here: https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp
specific_category = NULL # (i.e. "NABA" for ECM proteins) A string defining a specific gene set present in the MSigDB but is not separated with its own subcategory, 
design = sample_idesign_cptac_ccrcc_reduced
n_rotations = 999
minSetSize = 15
maxSetSize = 200
pvalueCutoff = 0.05

roastMSigDB <- function(data,
                        geneIDtype = "gene_symbol", # this can be "gene_symbol" or "entrez_gene"
                        organism = "Homo sapiens", # this can be any resulting from calling msigdbr::msigdbr_show_species()
                        category = "H", # Any of the main categories presented here: https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp
                        subcategory = NULL,
                        specific_category = NULL,
                        design,
                        n_rotations = 999,
                        minSetSize = 1,
                        maxSetSize = 200,
                        pvalueCutoff = 0.05){
      
      ## Load required packages ----  
      require(orgDB, character.only = TRUE) || stop(paste("package", organism, "is required", sep=" "))
      require(limma) || stop("Package limma is required")
      require(dplyr) || stop("Package dplyr is required")
      require(msigdbr) || stop("Package msigdbr is required")
      require(stringr) || stop("Package stringr is required")

      
      ## Generate matrix for roast ----
      
      if(geneIDtype != "SYMBOL" & geneIDtype != "ENTREZID"){
      genenames <- data$ID
      
      mapentrez <- clusterProfiler::bitr(geneID = genenames,
                                         fromType = geneIDtype,
                                         toType = "SYMBOL",
                                         OrgDb = eval(as.name(orgDB))) %>% na.omit()
      
      
      names(data) <- c(geneIDtype, names(data)[2:eval(dim(data)[2])])
      
      data <- dplyr::left_join(data,
                               mapentrez,
                               by = geneIDtype)
      
      }
      
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
         
         msigtab <- dplyr::filter(msigtab,
                                  str_detect(gs_name, paste0("^",specific_category,"_")))
         
      } 
      
      pathterm_n_iddf <- dplyr::select(msigtab,
                                       gs_name, gs_id) %>% distinct()
                        
      
      list_msig <- msigtab %>% split(x = .$gene_symbol, 
                                     f = .$gs_id)
      
      index <- limma::ids2indices(gene.sets = list_msig,
                                  identifiers = genesindata,
                                  remove.empty = TRUE)
      
      lenindex <- sapply(index, length)
      
      sublogi1 <- between(lenindex, minSetSize, maxSetSize) 
      index2 <- index[sublogi1]  

      genesinterm <- qdapTools::list2df(index2,
                                        col1 = "human_gene_symbol",
                                        col2 = "gs_id") %>%
                     dplyr::mutate(human_gene_symbol = as.character(human_gene_symbol))
      
      ## Run roast ----
      
      roast_out <- fry(y = matrix1,
                         contrast= ncol(design),
                         design = design, 
                         nrot = n_rotations, 
                         index = index2)
      
      # Process ROAST output ----
      suppressWarnings(suppressMessages(
         
         roast_out2 <- dplyr::mutate(roast_out,
                                     gs_id = row.names(roast_out)) %>% 
            dplyr::left_join(.,pathterm_n_iddf,
                             by = "gs_id") %>% 
            dplyr::rename(CategoryID = gs_id, CategoryTerm = gs_name)
         
      ))
      
      roast_out2 <- dplyr::filter(roast_out2,
                                  FDR <= pvalueCutoff)
      
      fdrnterm <- dplyr::select(roast_out2,
                                CategoryTerm, FDR, NGenes)
      
      genesintermread <- dplyr::filter(genesintermread,
                                       gs_id %in% roast_out2$CategoryID) %>%
                        dplyr::rename(ID = 3)
      
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
      
     # suppressWarnings(
      #   suppressMessages(
            log2FCs <- dplyr::mutate(limma_tab,
                                     ID = row.names(limma_tab))  #%>% 
               dplyr::filter(ID %in% genesintermread$ID)# %>% 
               dplyr::select(log2FC = eval(dim(.)[2]-5),ID) %>% 
               dplyr::left_join(., genesintermread, by = "ID") # %>% ## Look for a way to extract the z-score values as they are used by the ROAST algorithm
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

      
      
      return(roast_out)
      
}
                        