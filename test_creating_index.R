## Temp3: check how the ids2indices function works ----


library(limma)
library(qusage)
library(clusterProfiler)
library(org.Hs.eg.db)


msigdb <- qusage::read.gmt("c2.all.v7.0.symbols.gmt")
gomfsym <- qusage::read.gmt("c5.mf.v7.0.symbols.gmt") 
ggallsym <- qusage::read.gmt("c5.all.v7.0.symbols.gmt")


lapply(msigdb,
       head)[1:10]

lapply(gomfsym,
       head)[1:10]


datatab <- read.delim("input_limma_example.txt")

head(datatab)

datasymbols <- clusterProfiler::bitr(geneID = datatab$Name,
                                     fromType = "UNIPROT",
                                     toType = "SYMBOL",
                                     OrgDb = org.Hs.eg.db)

head(datasymbols)

testids2ind1 <- ids2indices(gene.sets = gomfsym,
                            identifiers = datasymbols$SYMBOL,
                            remove.empty = TRUE) 

head(testids2ind1)


## Test how to set a minimum and maximum set size ----

reactlistentrez <- as.list(reactomePATHID2EXTID)

index <- limma::ids2indices(gene.sets = reactlistentrez,
                            identifiers = genesindata,
                            remove.empty = TRUE) 

elemlen1 <- sapply(index, length)


max(elemlen1)
min(elemlen1)


clogi <- between(elemlen1, 200, 1000)


head(elemlen1)
head(clogi)

subindex <- index[clogi]


minset <- NA
maxset <- NULL

if(maxset > 0 & minset > 0){
      print("good")
} else {
      print("no good")}

