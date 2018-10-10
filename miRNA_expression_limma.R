library(limma)
library(statmod)

URL <- "https://www.ebi.ac.uk/arrayexpress/files/E-GEOD-16357"
SDRF.file <- "E-GEOD-16357.sdrf.txt"
Data.file <- "E-GEOD-16357.raw.1.zip"
download.file(paste(URL,SDRF.file,sep="/"), SDRF.file)
download.file(paste(URL,Data.file,sep="/"), Data.file)
unzip(Data.file)

SDRF <- read.delim("E-GEOD-16357.sdrf.txt",check.names=FALSE,stringsAsFactors=FALSE)
SDRF <- SDRF[c(1:11), ]
SDRF[,c("Array Data File","Characteristics[disease state]")]
Disease_state <- SDRF[c(1:11), "Characteristics[disease state]"]
str(Disease_state)
levels <- c("normal","Kaposi sarcoma")
Disease_state <- factor(Disease_state, levels=levels)
str(Disease_state)

x <- read.maimages(SDRF[,"Array Data File"],
                   source="agilent", green.only=TRUE, other.columns="gIsWellAboveBG")
dim(x)
str(x)

y <- backgroundCorrect(x, method="normexp")
y <- normalizeBetweenArrays(y, method="quantile")
Control <- y$genes$ControlType==1L
IsExpr <- rowSums(y$other$gIsWellAboveBG > 0) >= 4
yfilt <- y[!Control & IsExpr, ]
dim(yfilt)

design <- model.matrix(~Disease_state)
fit <- lmFit(yfilt,design)
fit <- eBayes(fit,trend=TRUE,robust=TRUE)
summary(decideTests(fit[,-1]))
topTable(fit,n=18)
