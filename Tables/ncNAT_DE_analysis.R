#PCA to determine the similiarity of the LIHC not tumor samples to the liver non-tumor and tumor LIHC samples
library("limma")
library('ggfortify')
library(readr)
library("DESeq2")
library("RColorBrewer")

prAnnotationTable <-read_delim(file='promoterCoordinates.gencode.v19_annotated.txt', delim='\t')
rownames(prAnnotationTable) <- prAnnotationTable$promoterId

load('ensembl.v75.annotation.Rdata')

#loads in the table associating the analysis with the tissue data and tumor status
sampleData.complete <- readRDS(file='sampleDescription_complete_matched.rds')

sample1="LIHC"
sample2="Liver"
isPromoter = FALSE

#loads in the expression of the different promotors/genes
abundance.TCGA <- read.delim(paste('data/absolute.TCGA.',sample1,'.gencode.v19.tsv',sep=''), row.names=1, sep=' ',quote='\"', check.names=F)
abundance.TCGA = abundance.TCGA[,names(abundance.TCGA) != 'geneId']
abundance.GTEx <- read.delim(paste('data/absolute.GTEx.',sample2,'.gencode.v19.tsv',sep=''), row.names=1, sep=' ',quote='\"', check.names=F)
abundance.GTEx = abundance.GTEx[,names(abundance.GTEx) != 'geneId']

#loads in the expression data for sidebars
abundanceGene.TCGA <- read.delim(paste('data/abundance.TCGA.',sample1,'.gencode.v19.tsv',sep=""), row.names=1, sep=' ',quote='\"', check.names=F)
abundanceGene.TCGA = abundanceGene.TCGA[,names(abundance.TCGA) != 'geneId']
abundanceGene.GTEx <- read.delim(paste('data/abundance.GTEx.',sample2,'.gencode.v19.tsv',sep=''), row.names=1, sep=' ',quote='\"', check.names=F)
abundanceGene.GTEx = abundanceGene.GTEx[,names(abundance.GTEx) != 'geneId']

#combines the LIHC and Liver data
abundance.TCGAGTEx <- cbind(abundance.TCGA,abundance.GTEx)
abundanceGene.TCGAGTEx <- cbind(abundanceGene.TCGA,abundanceGene.GTEx)

#gets the tissue and tumor status matching the LIHC and Liver data
sampleData.TCGA <- sampleData.complete[match(colnames(abundance.TCGA), sampleData.complete$analysisId),]
sampleData.GTEx <- sampleData.complete[match(colnames(abundance.GTEx), sampleData.complete$analysisId),]
sampleData.TCGAGTEx <- rbind(sampleData.TCGA,sampleData.GTEx)

ensemblAnnotationTxMatched <- (ensemblAnnotation[gsub('\\..*','',rownames(abundanceGene.TCGA)),])

# get all antisense RNA genes found in the gene biotype column in the Annotation files = "antisense"
if(isPromoter){
  abundance.TCGAGTEx.asRNA = abundance.TCGAGTEx[which(prAnnotationTable$gene_biotype=='antisense'),]
} else {
  abundance.TCGAGTEx.asRNA = abundance.TCGAGTEx[which(ensemblAnnotationTxMatched$gene_biotype=='antisense'),]
  
}
abundance.TCGAGTEx.asRNA = na.omit(abundance.TCGAGTEx.asRNA)

# get the samples which are tumours and which arn't (isTumor = TRUE/FALSE)
abundance.TCGAGTEx.tumor=abundance.TCGAGTEx[,which(sampleData.TCGAGTEx$isTumor==TRUE)]
abundance.TCGAGTEx.not_tumor=abundance.TCGAGTEx[,which(sampleData.TCGAGTEx$isTumor==FALSE)]
abundanceGene.TCGAGTEx.tumor=abundanceGene.TCGAGTEx[,which(sampleData.TCGAGTEx$isTumor==TRUE)]
abundanceGene.TCGAGTEx.not_tumor=abundanceGene.TCGAGTEx[,which(sampleData.TCGAGTEx$isTumor==FALSE)]

#quick check of HNF4A promoters
# prIndex.hnf4a <- which(prAnnotationTable$promoterId=='prmtr.11785')
# prIndex.hnf4a2 <- which(prAnnotationTable$promoterId=='prmtr.11786')
# prIndex.hnf4aas1 <- which(prAnnotationTable$promoterId=='prmtr.81103')
# prIndex.hnf4aas2  <- which(prAnnotationTable$promoterId=='prmtr.81104')
# geneIndex.hnf4aas1 <- which(ensemblAnnotationTxMatched$hgnc_symbol=='HNF4A-AS1')
# plot(as.double(abundance.TCGAGTEx.tumor[prIndex.hnf4a,]),as.double(abundance.TCGAGTEx.tumor[prIndex.hnf4a2,]), col=c("red"), xlab = "HNF4A P2", ylab ="HNF4A P1", main="LIHC red=tumor blue=not_tumor")
# points(as.double(abundance.TCGAGTEx.not_tumor[prIndex.hnf4a,]),as.double(abundance.TCGAGTEx.not_tumor[prIndex.hnf4a2,]), col=c("blue"))
# 
# plot(as.double(abundance.TCGAGTEx.tumor[prIndex.hnf4a,]),as.double(abundance.TCGAGTEx.tumor[prIndex.hnf4aas2,]), col=c("red"), xlab = "HNF4A P2", ylab ="HNF4A-AS P1")
# points(as.double(abundance.TCGAGTEx.not_tumor[prIndex.hnf4a,]),as.double(abundance.TCGAGTEx.not_tumor[prIndex.hnf4aas2,]), col=c("blue"))
# 
# plot(as.double(abundance.TCGAGTEx.tumor[prIndex.hnf4a,]),as.double(abundance.TCGAGTEx.tumor[prIndex.hnf4aas1,]), col=c("red"), xlab = "HNF4A P2", ylab ="HNF4A-AS P2")
# points(as.double(abundance.TCGAGTEx.not_tumor[prIndex.hnf4a,]),as.double(abundance.TCGAGTEx.not_tumor[prIndex.hnf4aas1,]), col=c("blue"))
# 
# plot(as.double(abundance.TCGAGTEx.tumor[prIndex.hnf4a2,]),as.double(abundance.TCGAGTEx.tumor[prIndex.hnf4aas2,]), col=c("red"), xlab = "HNF4A P1", ylab ="HNF4A-AS P1")
# points(as.double(abundance.TCGAGTEx.not_tumor[prIndex.hnf4a,]),as.double(abundance.TCGAGTEx.not_tumor[prIndex.hnf4aas2,]), col=c("blue"))
# 
# plot(as.double(abundance.TCGAGTEx.tumor[prIndex.hnf4a2,]),as.double(abundance.TCGAGTEx.tumor[prIndex.hnf4aas1,]), col=c("red"), xlab = "HNF4A P1", ylab ="HNF4A-AS P2")
# points(as.double(abundance.TCGAGTEx.not_tumor[prIndex.hnf4a,]),as.double(abundance.TCGAGTEx.not_tumor[prIndex.hnf4aas1,]), col=c("blue"))

#remove batch effects between the TCGA and liver data
# create the col data matrix listing if a sample is cancer or not
coldata=matrix(ncol = 2, nrow = length(colnames(abundance.TCGAGTEx)))
colnames(coldata)<-c('condition','batch')
rownames(coldata) <- colnames(abundance.TCGAGTEx)
coldata[unlist(sampleData.TCGAGTEx[sampleData.TCGAGTEx$isTumor==TRUE,1]),1] = 'cancer'
coldata[unlist(sampleData.TCGAGTEx[sampleData.TCGAGTEx$isTumor==FALSE,1]),1] = 'not_cancer'
coldata[unlist(sampleData.TCGAGTEx[sampleData.TCGAGTEx$study=='TCGA',1]),2] = 'TCGA'
coldata[unlist(sampleData.TCGAGTEx[sampleData.TCGAGTEx$study=='GTEx',1]),2] = 'GTEx'


#PCA plot

abundance.TCGAGTEx.nobatch=removeBatchEffect(abundance.TCGAGTEx, batch = coldata[,2], design=model.matrix(~coldata[,1]))

#testt=t(abundance.LIHCLiver[rowSums(abundance.LIHCLiver[,-1])>0,])
#testt=t(data.nobatch)
#data = t(abundance.TCGAGTEx[rowSums(abundance.TCGAGTEx[,-1])>0,])
data = t(abundance.TCGAGTEx)
data=data[, which(colSums(data) != 0)]
sampletype=rownames(data)
sampletype[is.element(sampletype,colnames(abundance.GTEx))]='GTEx'
sampletype[is.element(sampletype,colnames(abundance.TCGAGTEx.tumor))]='TCGA-tumor'
sampletype[is.element(sampletype,colnames(abundance.TCGAGTEx.not_tumor))]='TCGA non-tumor'

#res.pca <- prcomp(data, scale=TRUE)
testtdf=cbind(data,sampletype)
#autoplot(res.pca, data=testtdf,colour='sampletype')

#Differential Expression Analysis

dds <- DESeqDataSetFromMatrix(countData = na.omit(round(abundance.TCGAGTEx.asRNA)), #need to round because it is an import from Kalisto, rounding bias on read counts is minimal compared to biological and technical variation
                              colData = coldata,
                              design = ~ condition + batch)
dds <- DESeq(dds)
res <- results(dds, name="condition_not_cancer_vs_cancer")
resOrdered <- res[order(res$padj),]

write.csv(as.data.frame(resOrdered), 
          file="ncNAT_DE.csv")