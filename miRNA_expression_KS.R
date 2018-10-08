library(GEOquery)
GSE16357 <- getGEO('GSE16357', GSEMatrix = TRUE)
eset_miRNA <- GSE16357[[2]]
miRNA_expression_set <- exprs(eset_miRNA)
str(eset_miRNA)
gene.names  <- as.character(eset_miRNA@featureData@data[, "GENE_NAME"])
rownames(miRNA_expression_set) <- gene.names
colnames(miRNA_expression_set) <- c("skin-biopsy-patient1-repA", "skin-biopsy-patient2-repA", 
        "skin-biopsy-patient3-repA","AIDS_KS-biopsy-patient4-repA",	"AIDS_KS-biopsy-patient4-repB", 
        "AIDS_KS-biopsy-patient5-repA", "AIDS_KS-biopsy-patient6-repA", "AIDS_KS-biopsy-patient6-repB",
        "AIDS_KS-biopsy-patient7-repA", "AIDS_KS-biopsy-patient7-repB",	"AIDS_KS-biopsy-patient8-repA")
