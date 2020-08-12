#load in data
library(readr)

#load in the different cancer dataset acronyms
cancerDatasets <- read_csv(file='cancer_datasets.csv', col_names = TRUE)

prAnnotationTable <-read_delim(file='promoterCoordinates.gencode.v19_annotated.txt', delim='\t')
rownames(prAnnotationTable) <- prAnnotationTable$promoterId

load('ensembl.v75.annotation.Rdata')

#loads in the table associating the analysis with the tissue data and tumor status
sampleData.complete <- readRDS(file='sampleDescription_complete_matched.rds')

abundance=list()
abundance.tumor=list()
abundance.not_tumor=list()
abundanceGene=list()
abundanceGene.tumor=list()
abundanceGene.not_tumor=list()
sampleNames=c()

for (i in 1:nrow(cancerDatasets)){
  #loads in the expression of the different promotors/genes
  sample1=unlist(cancerDatasets[i,1])
  sample2=unlist(cancerDatasets[i,2])
  print(sample1)
  sampleNames=c(sampleNames,sample1)
  
  #sample1="LIHC"
  #sample2="Liver"
  if(is.na(sample2)){next}
  if(sample1=="DLBC" || sample1=="LAML" || sample1=="SKCM"){next}
  
  abundance.TCGA <- read.delim(paste('data/promoter_expression_TCGA/absolute.TCGA.',sample1,'.gencode.v19.tsv',sep=''), row.names=1, sep=' ',quote='\"', check.names=F)
  abundance.TCGA = abundance.TCGA[,names(abundance.TCGA) != 'geneId']
  abundance.GTEx <- read.delim(paste('data/promoter_expression_GTEx/absolute.GTEx.',sample2,'.gencode.v19.tsv',sep=''), row.names=1, sep=' ',quote='\"', check.names=F)
  abundance.GTEx = abundance.GTEx[,names(abundance.GTEx) != 'geneId']
  
  #loads in the expression data for sidebars
  abundanceGene.TCGA <- read.delim(paste('data/geneExpression/TCGA/kallisto-v0.44.0/biasCorrection.TRUE/abundance.TCGA.',sample1,'.gencode.v19.tsv',sep=""), row.names=1, sep=' ',quote='\"', check.names=F)
  abundanceGene.TCGA = abundanceGene.TCGA[,names(abundanceGene.TCGA) != 'geneId']
  abundanceGene.TCGA = abundanceGene.TCGA[,colnames(abundance.TCGA)]
  abundanceGene.GTEx <- read.delim(paste('data/geneExpression/GTEx/kallisto-v0.44.0/biasCorrection.TRUE/abundance.GTEx.',gsub(" ","", sample2, fixed = TRUE),'.gencode.v19.tsv',sep=''), row.names=1, sep=' ',quote='\"', check.names=F)
  abundanceGene.GTEx = abundanceGene.GTEx[,names(abundanceGene.GTEx) != 'geneId']
  abundanceGene.GTEx = abundanceGene.GTEx[,colnames(abundance.GTEx)]
  
  
  #combines the LIHC and Liver data
  abundance.TCGAGTEx <- cbind(abundance.TCGA,abundance.GTEx)
  abundanceGene.TCGAGTEx <- cbind(abundanceGene.TCGA,abundanceGene.GTEx)
  fraction=.1
  
  # get the samples which are tumours and which arn't (isTumor = TRUE/FALSE)
  #gets the tissue and tumor status matching the LIHC and Liver data
  sampleData.TCGA <- sampleData.complete[match(colnames(abundance.TCGA), sampleData.complete$analysisId),]
  sampleData.GTEx <- sampleData.complete[match(colnames(abundance.GTEx), sampleData.complete$analysisId),]
  sampleData.TCGAGTEx <- rbind(sampleData.TCGA,sampleData.GTEx)
  abundance.TCGAGTEx.tumor=abundance.TCGAGTEx[,which(sampleData.TCGAGTEx$isTumor==TRUE)]
  abundance.TCGAGTEx.not_tumor=abundance.TCGAGTEx[,which(sampleData.TCGAGTEx$isTumor==FALSE)]
  abundanceGene.TCGAGTEx.tumor=abundanceGene.TCGAGTEx[,which(sampleData.TCGAGTEx$isTumor==TRUE)]
  abundanceGene.TCGAGTEx.not_tumor=abundanceGene.TCGAGTEx[,which(sampleData.TCGAGTEx$isTumor==FALSE)]
  
  abundance[[i]]=abundance.TCGAGTEx
  abundance.tumor[[i]]=abundance.TCGAGTEx.tumor
  abundance.not_tumor[[i]]=abundance.TCGAGTEx.not_tumor
  abundanceGene[[i]]=abundanceGene.TCGAGTEx
  abundanceGene.tumor[[i]]=abundanceGene.TCGAGTEx.tumor
  abundanceGene.not_tumor[[i]]=abundanceGene.TCGAGTEx.not_tumor
}

ensemblAnnotationTxMatched <- (ensemblAnnotation[gsub('\\..*','',rownames(abundanceGene.TCGA)),])
