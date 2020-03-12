# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Thu Mar 12 13:52:19 EDT 2020

################################################################
#   Differential expression analysis with limma
library(BiocManager)
install()
library(Biobase)
library(GEOquery)
library(limma)

# load series and platform data from GEO

gset <- getGEO("GSE3251", GSEMatrix =TRUE, AnnotGPL=TRUE)
head(gset)
# $GSE3251_series_matrix.txt.gz
# ExpressionSet (storageMode: lockedEnvironment)
#   assayData: 15151 features, 16 samples 
#   element names: exprs 
# protocolData: none
# phenoData
#   sampleNames: GSM73236 GSM73237 ... GSM73251 (16 total)
#   varLabels: title geo_accession ... data_row_count (35 total)
#   varMetadata: labelDescription
# featureData
#   featureNames: H3001A01_1 H3001A02_1 ... H3159G07_1 (15151 total)
#   fvarLabels: ID Gene title ... GO:Component ID (21 total)
#   fvarMetadata: Column Description labelDescription
# experimentData: use 'experimentData(object)'
#   pubMedIds: 17079202 
# Annotation: GPL2813 

length(gset) # 1
if (length(gset) > 1) idx <- grep("GPL2813", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))
head(gset)

# group names for all samples
dim(pData(gset)) # 16  35
pData(gset)[,1:6]

group <- as.factor(rep(c("PDGF_6D", "noF_6D", "PDGF_12h", "noF_12h"), each = 4))
group

colnames(pData(gset))


#gsms <- "0001131103222323"
#sml <- c()
#for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# log2 transform
ex <- exprs(gset)
head(ex)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# set up the data and proceed with analysis
#gset$description <- group

boxplot(ex)
boxplot(exprs(gset))

library(RColorBrewer)
colors <- brewer.pal(4, "Dark2")
plotMDS(ex, col=colors[group])
plotMDS(exprs(gset), col=colors[group])

#library(oligo)
#gset.rma <- rma(gset)

#sml <- paste("G", sml, sep="")    # set group names
#fl <- as.factor(sml)
#gset$description <- fl
gset$description <- as.factor(c(rep("PDGF_6D",3), rep("noF_6D",2), "noF_12h", rep("noF_6D",2), "PDGF_6D", "noF_12h", rep("PDGF_12h",3), "noF_12h", "PDGF_12h", "noF_12h"))
gset$description

design <- model.matrix(~ 0 + description, gset)
design
colnames(design) <- gsub("description", "", colnames(design))
design

cont.matrix <- makeContrasts(
  neg12VSpos12 = noF_12h - PDGF_12h, 
  neg6dVSpos6d = noF_6D - PDGF_6D,
  levels= colnames(design))
cont.matrix

fit.des <- lmFit(gset, design)
fit.contr <- contrasts.fit(fit.des, cont.matrix)
#fit2 <- eBayes(fit2, 0.01)
efit.contr <- eBayes(fit.contr, proportion = 0.01)
summary(decideTests(efit.contr)) ## SIGNIFICANCE IS AT p.val < 0.05 rather than adj.p.val < 0.05
#       neg12VSpos12 neg6dVSpos6d
#Down              0            0
#NotSig        15151        15151
#Up                0            0

efit2.contr <- eBayes(fit.contr)
summary(decideTests(efit2.contr)) ## SIGNIFICANCE IS AT p.val < 0.05 rather than adj.p.val < 0.05
#       neg12VSpos12 neg6dVSpos6d
#Down              0            0
#NotSig        15151        15151
#Up                0            0

#tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)
res.neg12VSpos12 <- topTable(efit.contr, coef = 1, n = Inf)
head(res.neg12VSpos12)

#tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","F","Gene.symbol","Gene.title"))
res.neg12VSpos12 <- subset(res.neg12VSpos12, select=c("ID","Gene.symbol","logFC", "AveExpr", "t", "P.Value","adj.P.Val","B"))
head(res.neg12VSpos12)
write.table(res.neg12VSpos12, file="res.neg12VSpos12.txt", sep="\t",
            row.names=F, col.names = T, quote = F)

res.neg6dVSpos6d <- topTable(efit.contr, coef = 2, n = Inf)
head(res.neg6dVSpos6d)

#tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","F","Gene.symbol","Gene.title"))
res.neg6dVSpos6d <- subset(res.neg6dVSpos6d, select=c("ID","Gene.symbol","logFC", "AveExpr", "t", "P.Value","adj.P.Val","B"))
head(res.neg6dVSpos6d)
write.table(res.neg6dVSpos6d, file="res.neg6dVSpos6d.txt", sep="\t",
            row.names=F, col.names = T, quote = F)


################################################################
#   Boxplot for selected GEO samples
library(Biobase)
library(GEOquery)

# load series and platform data from GEO

gset <- getGEO("GSE3251", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL2813", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# group names for all samples in a series
gsms <- "0001131103222323"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
sml <- paste("G", sml, sep="")  # set group names

# order samples by group
ex <- exprs(gset)[ , order(sml)]
sml <- sml[order(sml)]
fl <- as.factor(sml)
labels <- c("PDGF_6d","noF_6d","PDGF_12h","noF_12h")

# set parameters and draw the plot
palette(c("#dfeaf4","#dfeaf4","#dfeaf4","#f4dfdf", "#AABBCC"))
dev.new(width=4+dim(gset)[[2]]/5, height=6)
par(mar=c(2+round(max(nchar(sampleNames(gset)))/2),4,2,1))
title <- paste ("GSE3251", '/', annotation(gset), " selected samples", sep ='')
boxplot(ex, boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gset$description)
legend("topleft", labels, fill=palette(), bty="n")
