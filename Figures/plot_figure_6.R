#HNF4A and HNF4A AS behavior in other cancers
#load in all the cancer data and tabulate all the data. 
library(readr)

sample1="LIHC"
sample2="Liver"
isPromoter = FALSE

#load in the different cancer dataset acronyms
cancerDatasets <- read_csv(file='cancer_datasets.csv', col_names = TRUE)

abundance.TCGA <- read.delim(paste('data/absolute.TCGA.',sample1,'.gencode.v19.tsv',sep=''), row.names=1, sep=' ',quote='\"', check.names=F)
abundance.TCGA = abundance.TCGA[,names(abundance.TCGA) != 'geneId']
abundanceGene.TCGA <- read.delim(paste('data/abundance.TCGA.',sample1,'.gencode.v19.tsv',sep=""), row.names=1, sep=' ',quote='\"', check.names=F)
abundanceGene.TCGA = abundanceGene.TCGA[,names(abundance.TCGA) != 'geneId']

prAnnotationTable <-read_delim(file='promoterCoordinates.gencode.v19_annotated.txt', delim='\t')
rownames(prAnnotationTable) <- prAnnotationTable$promoterId

load('ensembl.v75.annotation.Rdata')

#loads in the table associating the analysis with the tissue data and tumor status
sampleData.complete <- readRDS(file='sampleDescription_complete_matched.rds')

ensemblAnnotationTxMatched <- (ensemblAnnotation[gsub('\\..*','',rownames(abundanceGene.TCGA)),])

geneIndex.hnf4aas1 <- which(ensemblAnnotationTxMatched$hgnc_symbol=='HNF4A-AS1')
prIndex.hnf4a <- which(prAnnotationTable$promoterId=='prmtr.11785')
prIndex.hnf4a2 <- which(prAnnotationTable$promoterId=='prmtr.11786')
prIndex.hnf4aas1 <- which(prAnnotationTable$promoterId=='prmtr.81103')
prIndex.hnf4aas2  <- which(prAnnotationTable$promoterId=='prmtr.81104')

p1_all=c()
p1_tumor=c()
p1_non_tumor=c()
p2_all=c()
p2_tumor=c()
p2_non_tumor=c()
as_all=c()
as_tumor=c()
as_non_tumor=c()
samples=c()
samples_tumor = c() 
samples_non_tumor =c()
plotname=c()
all=c()
tumor=c()
non_tumor=c()

# abundance=readRDS("abundance.rds")
# abundance.not_tumor=readRDS("abundanceNotTumor.rds")
# abundance.tumor=readRDS("abundanceTumor.rds")
# abundanceGene=readRDS("abundanceGene.rds")
# abundanceGene.not_tumor=readRDS("abundanceGeneNotTumor.rds")
# abundanceGene.tumor=readRDS("abundanceGeneTumor.rds")
# sampleNames=readRDS("samplenames.rds")

for (i in 1:nrow(cancerDatasets)){
  #loads in the expression of the different promotors/genes
  sample1=unlist(cancerDatasets[i,1])
  sample2=unlist(cancerDatasets[i,2])
  print(sample1)
  if(is.na(sample2)){next}
  if(sample1=="DLBC" || sample1=="LAML"|| sample1=="SKCM"){next}
  
  if(exists("abundance")){
    #if data is already loaded
    abundance.TCGAGTEx=abundance[[i]]
    abundance.TCGAGTEx.tumor=abundance.tumor[[i]]
    abundance.TCGAGTEx.not_tumor=abundance.not_tumor[[i]]
    abundanceGene.TCGAGTEx=abundanceGene[[i]]
    abundanceGene.TCGAGTEx.tumor=abundanceGene.tumor[[i]]
    abundanceGene.TCGAGTEx.not_tumor=abundanceGene.not_tumor[[i]]
  } else {
  
    abundance.TCGA <- read.delim(paste('data/promoter_expression_TCGA/absolute.TCGA.',sample1,'.gencode.v19.tsv',sep=''), row.names=1, sep=' ',quote='\"', check.names=F)
    #abundance.TCGA <- read.delim(paste('data/absolute.TCGA.',sample1,'.gencode.v19.tsv',sep=''), row.names=1, sep=' ',quote='\"', check.names=F)
    abundance.TCGA = abundance.TCGA[,names(abundance.TCGA) != 'geneId']
    abundance.GTEx <- read.delim(paste('data/promoter_expression_GTEx/absolute.GTEx.',sample2,'.gencode.v19.tsv',sep=''), row.names=1, sep=' ',quote='\"', check.names=F)
    #abundance.GTEx <- read.delim(paste('data/absolute.GTEx.',sample2,'.gencode.v19.tsv',sep=''), row.names=1, sep=' ',quote='\"', check.names=F)
    abundance.GTEx = abundance.GTEx[,names(abundance.GTEx) != 'geneId']
    
    #loads in the expression data for sidebars
    abundanceGene.TCGA <- read.delim(paste('data/geneExpression/TCGA/kallisto-v0.44.0/biasCorrection.TRUE/abundance.TCGA.',sample1,'.gencode.v19.tsv',sep=""), row.names=1, sep=' ',quote='\"', check.names=F)
    #abundanceGene.TCGA <- read.delim(paste('data/abundance.TCGA.',sample1,'.gencode.v19.tsv',sep=""), row.names=1, sep=' ',quote='\"', check.names=F)
    abundanceGene.TCGA = abundanceGene.TCGA[,names(abundanceGene.TCGA) != 'geneId']
    abundanceGene.TCGA = abundanceGene.TCGA[,colnames(abundance.TCGA)]
    abundanceGene.GTEx <- read.delim(paste('data/geneExpression/GTEx/kallisto-v0.44.0/biasCorrection.TRUE/abundance.GTEx.',gsub(" ","", sample2, fixed = TRUE),'.gencode.v19.tsv',sep=''), row.names=1, sep=' ',quote='\"', check.names=F)
    #abundanceGene.GTEx <- read.delim(paste('data/abundance.GTEx.',sample2,'.gencode.v19.tsv',sep=''), row.names=1, sep=' ',quote='\"', check.names=F)
    abundanceGene.GTEx = abundanceGene.GTEx[,names(abundanceGene.GTEx) != 'geneId']
    abundanceGene.GTEx = abundanceGene.GTEx[,colnames(abundance.GTEx)]
    
    #combines the LIHC and Liver data
    abundance.TCGAGTEx <- cbind(abundance.TCGA,abundance.GTEx)
    abundanceGene.TCGAGTEx <- cbind(abundanceGene.TCGA,abundanceGene.GTEx)
    
    # get the samples which are tumours and which arn't (isTumor = TRUE/FALSE)
    #gets the tissue and tumor status matching the LIHC and Liver data
    sampleData.TCGA <- sampleData.complete[match(colnames(abundance.TCGA), sampleData.complete$analysisId),]
    sampleData.GTEx <- sampleData.complete[match(colnames(abundance.GTEx), sampleData.complete$analysisId),]
    sampleData.TCGAGTEx <- rbind(sampleData.TCGA,sampleData.GTEx)
    abundance.TCGAGTEx.tumor=abundance.TCGAGTEx[,which(sampleData.TCGAGTEx$isTumor==TRUE)]
    abundance.TCGAGTEx.not_tumor=abundance.TCGAGTEx[,which(sampleData.TCGAGTEx$isTumor==FALSE)]
    abundanceGene.TCGAGTEx.tumor=abundanceGene.TCGAGTEx[,which(sampleData.TCGAGTEx$isTumor==TRUE)]
    abundanceGene.TCGAGTEx.not_tumor=abundanceGene.TCGAGTEx[,which(sampleData.TCGAGTEx$isTumor==FALSE)]
  }



  #get HN4FA promoter average expression
  #p1

  if(mean(as.double(abundance.TCGAGTEx[prIndex.hnf4a2,]))<1 && mean(as.double(abundance.TCGAGTEx[prIndex.hnf4a,]))<1 && mean(as.double(abundanceGene.TCGAGTEx[geneIndex.hnf4aas1,]))<1){next}
  samples = c(samples,rep.int(sample1,length(abundance.TCGAGTEx[prIndex.hnf4a2,])))
  samples_tumor = c(samples_tumor,rep.int(sample1,length(abundance.TCGAGTEx.tumor[prIndex.hnf4a2,])))
  samples_non_tumor = c(samples_non_tumor,rep.int(sample1,length(abundance.TCGAGTEx.not_tumor[prIndex.hnf4a2,])))
  p1_all = c(p1_all,as.double(abundance.TCGAGTEx[prIndex.hnf4a2,]))
  p1_tumor = c(p1_tumor,as.double(abundance.TCGAGTEx.tumor[prIndex.hnf4a2,]))
  p1_non_tumor = c(p1_non_tumor,as.double(abundance.TCGAGTEx.not_tumor[prIndex.hnf4a2,]))
  p2_all = c(p2_all,as.double(abundance.TCGAGTEx[prIndex.hnf4a,]))
  p2_tumor = c(p2_tumor,as.double(abundance.TCGAGTEx.tumor[prIndex.hnf4a,]))
  p2_non_tumor = c(p2_non_tumor,as.double(abundance.TCGAGTEx.not_tumor[prIndex.hnf4a,]))
  as_all = c(as_all,log(as.double(abundanceGene.TCGAGTEx[geneIndex.hnf4aas1,])+1))
  as_tumor = c(as_tumor,log(as.double(abundanceGene.TCGAGTEx.tumor[geneIndex.hnf4aas1,])+1))
  as_non_tumor = c(as_non_tumor,log(as.double(abundanceGene.TCGAGTEx.not_tumor[geneIndex.hnf4aas1,])+1))
  
  #png(paste("draft_figures/cancer_vs_normal_gene/",sample1,"_",sample2,"_HNF4A_AS_gene.png",sep=""), width=1200, height=1200)
  svg(paste("draft_figures/cancer_vs_normal_gene/",sample1,"_",sample2,"_HNF4A_AS_gene.svg",sep=""), width=12, height=12)
  par(mfrow=c(2,2))
  
  #plot cancer vs normal for each cancer type showing relationship between AS and P1 
  plot(as.double(log(abundanceGene.TCGAGTEx.tumor[geneIndex.hnf4aas1,])), as.double(abundance.TCGAGTEx.tumor[prIndex.hnf4a2,]), xlab='LOG HNF4A-AS gene expression', ylab='HNF4A promoter 1 activity', col=c("lightcoral"), pch=16, cex=1.3)
  points(as.double(log(abundanceGene.TCGAGTEx.not_tumor[geneIndex.hnf4aas1,])), as.double(abundance.TCGAGTEx.not_tumor[prIndex.hnf4a2,]), col=c("steelblue2"), pch=16, cex=1.3)
  
  legend("topleft", legend=c("Tumor","Non-tumor"),fill=c("lightcoral", "steelblue2"), bty='n')
  
  correlation=cor(as.double(cbind(log(abundanceGene.TCGAGTEx.tumor[geneIndex.hnf4aas1,]+1),log(abundanceGene.TCGAGTEx.not_tumor[geneIndex.hnf4aas1,]+1))),as.double(cbind(abundance.TCGAGTEx.tumor[prIndex.hnf4a2,],abundance.TCGAGTEx.not_tumor[prIndex.hnf4a2,])), method = "pearson")
  text(par('usr')[1]+1,y=par('usr')[3]+3,labels=as.character(signif(correlation,3)), col='black')
  correlation=cor(as.double(log(abundanceGene.TCGAGTEx.tumor[geneIndex.hnf4aas1,]+1)),as.double(abundance.TCGAGTEx.tumor[prIndex.hnf4a2,]), method = "pearson")
  text(par('usr')[1]+1,y=par('usr')[3]+2,labels=as.character(signif(correlation,3)), col='lightcoral')
  correlation=cor(as.double(log(abundanceGene.TCGAGTEx.not_tumor[geneIndex.hnf4aas1,]+1)),as.double(abundance.TCGAGTEx.not_tumor[prIndex.hnf4a2,]), method = "pearson")
  text(par('usr')[1]+1,y=par('usr')[3]+1,labels=as.character(signif(correlation,3)), col='steelblue2')
  
  
  #plot cancer vs normal for each cancer type showing relationship between AS and P2
  plot(log(as.double(abundanceGene.TCGAGTEx.tumor[geneIndex.hnf4aas1,])), as.double(abundance.TCGAGTEx.tumor[prIndex.hnf4a,]), xlab='LOG HNF4A-AS gene expression', ylab='HNF4A promoter 2 activity', col=c("lightcoral"), pch=16, cex=1.3)
  points(log(as.double(abundanceGene.TCGAGTEx.not_tumor[geneIndex.hnf4aas1,])), as.double(abundance.TCGAGTEx.not_tumor[prIndex.hnf4a,]), col=c("steelblue2"), pch=16, cex=1.3)
  
  correlation=cor(as.double(cbind(log(abundanceGene.TCGAGTEx.tumor[geneIndex.hnf4aas1,]+1),log(abundanceGene.TCGAGTEx.not_tumor[geneIndex.hnf4aas1,]+1))),as.double(cbind(abundance.TCGAGTEx.tumor[prIndex.hnf4a,],abundance.TCGAGTEx.not_tumor[prIndex.hnf4a,])), method = "pearson")
  text(par('usr')[1]+1,y=par('usr')[3]+3,labels=as.character(signif(correlation,3)), col='black')
  correlation=cor(as.double(log(abundanceGene.TCGAGTEx.tumor[geneIndex.hnf4aas1,]+1)),as.double(abundance.TCGAGTEx.tumor[prIndex.hnf4a,]), method = "pearson")
  text(par('usr')[1]+1,y=par('usr')[3]+2,labels=as.character(signif(correlation,3)), col='lightcoral')
  correlation=cor(as.double(log(abundanceGene.TCGAGTEx.not_tumor[geneIndex.hnf4aas1,]+1)),as.double(abundance.TCGAGTEx.not_tumor[prIndex.hnf4a,]), method = "pearson")
  text(par('usr')[1]+1,y=par('usr')[3]+1,labels=as.character(signif(correlation,3)), col='steelblue2')
  
  #plot cancer vs normal for each cancer type comparing p1 and p2
  plot(as.double(abundance.TCGAGTEx.tumor[prIndex.hnf4a2,]), as.double(abundance.TCGAGTEx.tumor[prIndex.hnf4a,]), xlab='HNF4A promoter 1 activity', ylab='HNF4A promoter 2 activity', col=c("lightcoral"), pch=16, cex=1.3)
  points(as.double(abundance.TCGAGTEx.not_tumor[prIndex.hnf4a2,]), as.double(abundance.TCGAGTEx.not_tumor[prIndex.hnf4a,]), col=c("steelblue2"), pch=16, cex=1.3)
  
  correlation=cor(as.double(cbind(abundance.TCGAGTEx.tumor[prIndex.hnf4a2,],abundance.TCGAGTEx.not_tumor[prIndex.hnf4a2,])),as.double(cbind(abundance.TCGAGTEx.tumor[prIndex.hnf4a,],abundance.TCGAGTEx.not_tumor[prIndex.hnf4a,])), method = "pearson")
  text(par('usr')[1]+1,y=par('usr')[3]+3,labels=as.character(signif(correlation,3)), col='black')
  correlation=cor(as.double(abundance.TCGAGTEx.tumor[prIndex.hnf4a2,]),as.double(abundance.TCGAGTEx.tumor[prIndex.hnf4a,]), method = "pearson")
  text(par('usr')[1]+1,y=par('usr')[3]+2,labels=as.character(signif(correlation,3)), col='lightcoral')
  correlation=cor(as.double(abundance.TCGAGTEx.not_tumor[prIndex.hnf4a2,]),as.double(abundance.TCGAGTEx.not_tumor[prIndex.hnf4a,]), method = "pearson")
  text(par('usr')[1]+1,y=par('usr')[3]+1,labels=as.character(signif(correlation,3)), col='steelblue2')
  
  #plot cancer vs normal for each cancer type comparing AS and the P1-P2 ratio
  p1t=as.double(abundance.TCGAGTEx.tumor[prIndex.hnf4a2,])
  p2t=as.double(abundance.TCGAGTEx.tumor[prIndex.hnf4a,])
  plot(log(as.double(abundanceGene.TCGAGTEx.tumor[geneIndex.hnf4aas1,])),p1t/(p1t+p2t), xlab='LOG HNF4A-AS gene expression', ylab='P1/(P1+P2)', col=c("lightcoral"), pch=16, cex=1.3)
  p1nt=as.double(abundance.TCGAGTEx.not_tumor[prIndex.hnf4a2,])
  p2nt=as.double(abundance.TCGAGTEx.not_tumor[prIndex.hnf4a,])
  points(log(as.double(abundanceGene.TCGAGTEx.not_tumor[geneIndex.hnf4aas1,])),p1nt/(p1nt+p2nt), col=c("steelblue2"), pch=16, cex=1.3)
  
  dev.off()
} 

library(ggplot2)
library(ggpubr)

allname=c(rep.int("p1",length(p1_all)),rep.int("p2",length(p2_all)),rep.int("as",length(as_all)))
samplescombined=c(samples,samples,samples)  
all=c(p1_all,p2_all,as_all)
df = data.frame(samplescombined, all, allname)
plot1=ggplot(df, aes(x=samplescombined, y =all, fill=allname)) + 
  geom_violin(position=position_dodge(1), scale = "width", trim=FALSE) + 
  geom_boxplot(width=0.1, position=position_dodge(1)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

tumor=c(p1_tumor,p2_tumor,as_tumor)
samplescombined=c(samples_tumor,samples_tumor,samples_tumor)  
allname=c(rep.int("p1",length(p1_tumor)),rep.int("p2",length(p2_tumor)),rep.int("as",length(as_tumor)))
df = data.frame(samplescombined, tumor, allname)
plot2=ggplot(df, aes(x=samplescombined, y =tumor, fill=allname)) + 
  geom_violin(position=position_dodge(1), scale = "width", trim=FALSE) + 
  geom_boxplot(width=0.1, position=position_dodge(1)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

non_tumor=c(p1_non_tumor,p2_non_tumor,as_non_tumor)
samplescombined=c(samples_non_tumor,samples_non_tumor,samples_non_tumor)  
allname=c(rep.int("p1",length(p1_non_tumor)),rep.int("p2",length(p2_non_tumor)),rep.int("as",length(as_non_tumor)))
df = data.frame(samplescombined, non_tumor, allname)
plot3=ggplot(df, aes(x=samplescombined, y =non_tumor, fill=allname)) + 
  geom_violin(position=position_dodge(1), scale = "width", trim=FALSE) + 
  geom_boxplot(width=0.1, position=position_dodge(1)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggarrange(plot1,plot2,plot3, labels =c("all","tumor","non_tumor"), ncol=1,nrow=3)

