# other asRNA promoter switches
library(readr)
library("DESeq2")

asThreshold = 1
prThreshold = 1
pThreshold = 0.001
fcThreshold = 1.5
corThreshold = .5

#load in data
prAnnotationTable <-read_delim(file='promoterCoordinates.gencode.v19_annotated.txt', delim='\t')
rownames(prAnnotationTable) <- prAnnotationTable$promoterId

load('ensembl.v75.annotation.Rdata')

#loads in the table associating the analysis with the tissue data and tumor status
sampleData.complete <- readRDS(file='sampleDescription_complete_matched.rds')

#load in the different cancer dataset acronyms
cancerDatasets <- read_csv(file='cancer_datasets.csv', col_names = TRUE)

sample1="LIHC"
isPromoter=FALSE
abundance.TCGA <- read.delim(paste('data/promoter_expression_TCGA/absolute.TCGA.',sample1,'.gencode.v19.tsv',sep=''), row.names=1, sep=' ',quote='\"', check.names=F)
abundance.TCGA.asRNA = abundance.TCGA[which(prAnnotationTable$gene_biotype=='antisense'),]
abundanceGene.TCGA <- read.delim(paste('data/geneExpression/TCGA/kallisto-v0.44.0/biasCorrection.TRUE/abundance.TCGA.',sample1,'.gencode.v19.tsv',sep=""), row.names=1, sep=' ',quote='\"', check.names=F)

ensemblAnnotationTxMatched <- (ensemblAnnotation[gsub('\\..*','',rownames(abundanceGene.TCGA)),])

isCandidate <- function(pr){
  #check its expressed in tissue
  if(mean(log(as.double(abundanceGene.TCGAGTEx.asRNA[pr,])+1))<asThreshold){
    return(list(p1=NA, p2=NA, p1t=NA, p2t=NA, p1nt=NA, p2nt=NA))
  }
  
  #get cognate gene
  if (isPromoter){
    asGene=prAnnotationTable[pr,9]
  } else {
    asGene=pr
  }
  
  cognateGene = as.character(sgRNAPairs[unlist(asGene),3])
  #get gene promoters of the cognate gene
  cognatePromoters=prAnnotationTable[which(prAnnotationTable$ensembl_gene_id==cognateGene),]
  
  promoter1 = NA
  promoter2 = NA
  promoter1t = NA
  promoter2t = NA
  promoter1nt = NA
  promoter2nt = NA
  
  #if it has more than one promoter
  if (nrow(cognatePromoters)>1){
    
    #check that it has only 2 active promoters
    activePromoterCount=0
    activePromoterCountt=0
    activePromoterCountnt=0
    for (j in 1:nrow(cognatePromoters)){
      mean1=mean(as.double(abundance.tumor[[i]][unlist(cognatePromoters[j,6]),]))
      mean2= mean(as.double(abundance.not_tumor[[i]][unlist(cognatePromoters[j,6]),]))
      #meanall = mean(as.double(abundance[[i]][unlist(cognatePromoters[j,6]),]))
      measuredValues = length(which((as.double(abundance[[i]][unlist(cognatePromoters[j,6]),]))>0))
      measuredValuest = length(which((as.double(abundance.tumor[[i]][unlist(cognatePromoters[j,6]),]))>0))
      measuredValuesnt = length(which((as.double(abundance.not_tumor[[i]][unlist(cognatePromoters[j,6]),]))>0))
      if (is.na(measuredValues)){measuredValues=0}
      if (is.na(measuredValuest)){measuredValuest=0}
      if (is.na(measuredValuesnt)){measuredValuesnt=0}
      if (is.na(mean1)){mean1=0}
      if (is.na(mean2)){mean2=0}
      if ((mean1>prThreshold || mean2>prThreshold) && measuredValues >= 10 ){
        activePromoterCount = activePromoterCount+1
        if(activePromoterCount==1){promoter1=unlist(cognatePromoters[j,6])}
        if(activePromoterCount==2){promoter2=unlist(cognatePromoters[j,6])}
        if(activePromoterCount>2){
          promoter1=NA
          promoter2=NA
          }
      }
      if (mean1>prThreshold && measuredValuest >= 10){
        activePromoterCountt = activePromoterCountt+1
        if(activePromoterCountt==1){promoter1t=unlist(cognatePromoters[j,6])}
        if(activePromoterCountt==2){promoter2t=unlist(cognatePromoters[j,6])}
        if(activePromoterCountt>2){
          promoter1t=NA
          promoter2t=NA
          }
      }
      if (mean2>prThreshold && measuredValuesnt >= 10){
        activePromoterCountnt = activePromoterCountnt+1
        if(activePromoterCountnt==1){promoter1nt=unlist(cognatePromoters[j,6])}
        if(activePromoterCountnt==2){promoter2nt=unlist(cognatePromoters[j,6])}
        if(activePromoterCountnt>2){
          promoter1nt=NA
          promoter2nt=NA
        }
      }
    }
  }
  return(list(p1=promoter1, p2=promoter2, p1t=promoter1t, p2t=promoter2t, p1nt=promoter1nt, p2nt=promoter2nt))
}

plotCor <- function(asRNA=NA, p1data=NA, p2data=NA, asRNAnt=NA, p1datant=NA, p2datant=NA){
  #plot cancer vs normal for each cancer type showing relationship between AS and P1 
  
  # asRNA = as.double(log(abundanceGene.tumor[[i]][pr,]))
  # p1data = as.double(abundance.tumor[[i]][promoter1,])
  # p2data = as.double(abundance.tumor[[i]][promoter2,])
  # asRNAnt = as.double(log(abundanceGene.not_tumor[[i]][pr,]))
  # p1datant =  as.double(abundance.not_tumor[[i]][promoter1,])
  # p2datant =  as.double(abundance.not_tumor[[i]][promoter2,])
  
  plot(asRNA, p1data, xlab=paste(pr,' (log) gene expression',sep = ""), ylab=paste(promoter1,' promoter activity',sep = ""), col=c("lightcoral"), pch=16, cex=1.3, 
       xlim=c(min(c(asRNA[!is.na(asRNA) & !is.infinite(asRNA)],asRNAnt[!is.na(asRNAnt) & !is.infinite(asRNAnt)])),
              max(c(asRNA[!is.na(asRNA) & !is.infinite(asRNA)],asRNAnt[!is.na(asRNAnt) & !is.infinite(asRNAnt)]))), 
       ylim=c(min(c(p1data[!is.na(p1data) & !is.infinite(p1data)],p1datant[!is.na(p1datant) & !is.infinite(p1datant)])),
              max(c(p1data[!is.na(p1data) & !is.infinite(p1data)],p1datant[!is.na(p1datant) & !is.infinite(p1datant)]))))
  points(asRNAnt, p1datant, col=c("steelblue2"), pch=16, cex=1.3)
  
  plot(asRNA, p2data, xlab=paste(pr,' (log) gene expression',sep = ""), ylab=paste(promoter2,' promoter activity',sep = ""), col=c("lightcoral"), pch=16, cex=1.3,
       xlim=c(min(c(asRNA[!is.na(asRNA) & !is.infinite(asRNA)],asRNAnt[!is.na(asRNAnt) & !is.infinite(asRNAnt)])),
              max(c(asRNA[!is.na(asRNA) & !is.infinite(asRNA)],asRNAnt[!is.na(asRNAnt) & !is.infinite(asRNAnt)]))), 
       ylim=c(min(c(p2data[!is.na(p2data) & !is.infinite(p2data)],p2datant[!is.na(p2datant) & !is.infinite(p2datant)])),
              max(c(p2data[!is.na(p2data) & !is.infinite(p2data)],p2datant[!is.na(p2datant) & !is.infinite(p2datant)]))))
  points(asRNAnt, p2datant, col=c("steelblue2"), pch=16, cex=1.3)
  
  plot(p1data, p2data, xlab=paste(promoter1,' promoter activity',sep = ""), ylab=paste(promoter2,' promoter activity',sep = ""), col=c("lightcoral"), pch=16, cex=1.3, 
       xlim=c(min(c(p1data[!is.na(p1data) & !is.infinite(p1data)],p1datant[!is.na(p1datant) & !is.infinite(p1datant)])),
              max(c(p1data[!is.na(p1data) & !is.infinite(p1data)],p1datant[!is.na(p1datant) & !is.infinite(p1datant)]))), 
       ylim=c(min(c(p2data[!is.na(p2data) & !is.infinite(p2data)],p2datant[!is.na(p2datant) & !is.infinite(p2datant)])),
              max(c(p2data[!is.na(p2data) & !is.infinite(p2data)],p2datant[!is.na(p2datant) & !is.infinite(p2datant)]))))
  
  points(p1datant, p2datant, col=c("steelblue2"), pch=16, cex=1.3)
  
  legend("topleft", legend=c("Tumor","Non-tumor"),fill=c("lightcoral", "steelblue2"), bty='n')
  
  # correlation=cor(as.double(cbind(log(abundanceGene.TCGAGTEx.tumor[geneID,]+1),log(abundanceGene.TCGAGTEx.not_tumor[geneID,]+1))),as.double(cbind(abundance.TCGAGTEx.tumor[p1ID,],abundance.TCGAGTEx.not_tumor[p1ID,])), method = "pearson")
  # text(par('usr')[1]+1,y=par('usr')[3]+3,labels=as.character(signif(correlation,3)), col='black')
  # correlation=cor(as.double(log(abundanceGene.TCGAGTEx.tumor[geneID,]+1)),as.double(abundance.TCGAGTEx.tumor[p1ID,]), method = "pearson")
  # text(par('usr')[1]+1,y=par('usr')[3]+2,labels=as.character(signif(correlation,3)), col='lightcoral')
  # correlation=cor(as.double(log(abundanceGene.TCGAGTEx.not_tumor[geneID,]+1)),as.double(abundance.TCGAGTEx.not_tumor[p1ID,]), method = "pearson")
  # text(par('usr')[1]+1,y=par('usr')[3]+1,labels=as.character(signif(correlation,3)), col='steelblue2')
}

testCandidate <- function(pr, abundanceGene, abundance, promoter1, promoter2, i){
  #get lowest 25% of asRNA values and matching P1/P2 values from all/tumor/non_tumor
  #sort by as promoter row
  sig = FALSE
  fullOutput=c(pr, unlist(sampleNames[i]), promoter1, promoter2, unlist(prAnnotationTable[promoter1,'hgnc_symbol']), -1, -1, -1, -1, -1, -1, 2, 2)
  if (is.null(abundanceGene[[i]][pr,])){
    output=c(pr, unlist(sampleNames[i]), promoter1, promoter2, unlist(prAnnotationTable[promoter1,'hgnc_symbol']), -1, -1, 2)
    return(list(results=output, lowp1=c(), highp1=c(), lowp2=c(), highp2=c(), sig=sig, fullOutput=fullOutput))
  }
  abundance.TCGAGTEx.sorted=abundanceGene[[i]][,order(abundanceGene[[i]][pr,])]
  #get first 25% of col names
  limit1=round(ncol(abundanceGene[[i]])*fraction)+1
  limit2=ncol(abundanceGene[[i]])-round(ncol(abundanceGene[[i]])*fraction)
  if(length(1:limit1)<3){limit1=3}
  if(length(limit2:ncol(abundanceGene[[i]]))<3){limit2=ncol(abundanceGene[[i]])-2}
  lowCols=colnames(abundance.TCGAGTEx.sorted[pr,1:limit1])
  highCols=colnames(abundance.TCGAGTEx.sorted[pr,limit2:ncol(abundanceGene[[i]])])
  #get those col names of p1
  lowp1=as.double(abundance[[i]][promoter1,lowCols])
  highp1=as.double(abundance[[i]][promoter1,highCols])
  #p2
  lowp2=as.double(abundance[[i]][promoter2,lowCols])
  highp2=as.double(abundance[[i]][promoter2,highCols])
  
  lowp1=lowp1[!is.na(lowp1)]
  highp1=highp1[!is.na(highp1)]
  lowp2=lowp2[!is.na(lowp2)]
  highp2=highp2[!is.na(highp2)]
  
  output=c(pr, unlist(sampleNames[i]), promoter1, promoter2, unlist(prAnnotationTable[promoter1,'hgnc_symbol']), mean(highp1)/mean(lowp1), mean(highp2)/mean(lowp2), 2)
  if(length(lowp1)>3 & length(highp1)>3 & length(lowp2)>3 & length(highp2)>3){
    #calculate mean difference
    #p1p=wilcox.test(highp1, lowp1, alternative = "greater")
    p1p=wilcox.test(highp1, lowp1)
    #p1p$p.value)
    #p2p=wilcox.test(highp2, lowp2, alternative = "greater")
    p2p=wilcox.test(highp2, lowp2)
    
    cor1=cor(unlist(log(abundanceGene[[i]][pr,]+1)), unlist(abundance[[i]][promoter1,]), method = "pearson", use="complete.obs")
    cor2=cor(unlist(log(abundanceGene[[i]][pr,]+1)), unlist(abundance[[i]][promoter2,]), method = "pearson", use="complete.obs")
    fc1=mean(highp1)/mean(lowp1)
    fc2=mean(highp2)/mean(lowp2)
    
    fullOutput=c(pr, unlist(sampleNames[i]), promoter1, promoter2, unlist(prAnnotationTable[promoter1,'hgnc_symbol']), mean(highp1), mean(lowp1), mean(highp2), mean(lowp2), p1p$p.value, p2p$p.value, cor1, cor2)
    
    if(is.na(p1p$p.value) | is.na(p2p$p.value) | is.na(fc1) | is.na(cor1)){
      output=c(pr, unlist(sampleNames[i]), promoter1, promoter2, unlist(prAnnotationTable[promoter1,'hgnc_symbol']), -1, -1, 2)
      return(list(results=output, lowp1=lowp1, highp1=highp1, lowp2=lowp2, highp2=highp2, sig=sig, fullOutput=fullOutput))
    }
    
    #candidate #tissue #p1 #p2 #difference  #pvalue
    #+/0 -/0
    if(p1p$p.value<pThreshold & p2p$p.value>pThreshold & (fc1 > fcThreshold | fc1 < (1/fcThreshold)) & (cor1 > corThreshold | cor1 < -corThreshold)){
      sig="P1"
      #print("t p1")
      output=c(pr, unlist(sampleNames[i]), promoter1, promoter2, unlist(prAnnotationTable[promoter1,'hgnc_symbol']), mean(highp1), mean(lowp1), mean(highp2), mean(lowp2), p1p$p.value, p2p$p.value, cor1, cor2)
    } #0/+ 0/- 
    else if(p2p$p.value<pThreshold & p1p$p.value>pThreshold & (fc2 > fcThreshold | fc2 < (1/fcThreshold)) & (cor2 > corThreshold | cor2 < -corThreshold)){
      sig="P2"
      #print("t p2")
      output=c(pr, unlist(sampleNames[i]), promoter2, promoter1, unlist(prAnnotationTable[promoter1,'hgnc_symbol']), mean(highp2)/mean(lowp2), mean(highp1)/mean(lowp1), p2p$p.value, p1p$p.value, cor2, cor1)
    }else {
      output=c(pr, unlist(sampleNames[i]), promoter1, promoter2, unlist(prAnnotationTable[promoter1,'hgnc_symbol']), mean(highp1),mean(lowp1), mean(highp2), mean(lowp2), p1p$p.value, p2p$p.value, cor1, cor2)
    }
    #p=boxplot(lowp1,highp1, lowp2, highp2, names=c("low AS (P1)", "high AS (P1)", "low AS (P2)", "high AS (P2)"), at = c(1,2,4,5), col = c("orange","red"), ylab = "(promoter expression)", main=paste(sample1," all", as.character(sig)))
  }
  return(list(results=output, lowp1=lowp1, highp1=highp1, lowp2=lowp2, highp2=highp2, sig=sig, fullOutput=fullOutput))
}

#candidate #tissue #p1 #p2 #difference  #pvalue 
results=matrix(ncol=8)
colnames(results)=c('candidate','tissue','p1','p2', 'gene name', 'P1 mean FC', 'P2 mean FC','pvalue')

results.tumor=matrix(ncol=8)
colnames(results.tumor)=c('candidate','tissue','p1','p2','gene name', 'P1 mean FC', 'P2 mean FC','pvalue')

results.not_tumor=matrix(ncol=8)
colnames(results.not_tumor)=c('candidate','tissue','p1','p2','gene name', 'P1 mean FC', 'P2 mean FC','pvalue')

resultsAll=matrix(ncol=13)
colnames(resultsAll)=c('candidate','tissue','p1','p2', 'gene name', 'P1 low AS', 'P1 high AS', 'P2 low AS', 'P2 high AS', 'pvalue', 'cor1', 'cor2')

resultsAll.tumor=matrix(ncol=13)
colnames(resultsAll.tumor)=c('candidate','tissue','p1','p2', 'gene name', 'P1 low AS', 'P1 high AS', 'P2 low AS', 'P2 high AS', 'pvalue', 'cor1', 'cor2')

resultsAll.not_tumor=matrix(ncol=13)
colnames(resultsAll.not_tumor)=c('candidate','tissue','p1','p2', 'gene name', 'P1 low AS', 'P1 high AS', 'P2 low AS', 'P2 high AS', 'pvalue', 'cor1', 'cor2')

heatmapMatrix = matrix(ncol = 8)
colnames(heatmapMatrix) = c("gene","sample","FC_P1_All", "FC_P1_tumor", "FC_P1_not_tumor", "FC_P2_All", "FC_P2_tumor", "FC_P2_not_tumor")
#if data is already retrieved 

hitIDs=c()
hitIDst=c()
hitIDsnt=c()

sgRNAPairs=read.csv("sgRNA_cognate_pairs.csv")
rownames(sgRNAPairs)<-sgRNAPairs[,2]

for (i in 1:length(abundance)){
#for (i in 16:16){
  if(is.null(abundance[[i]])){next}
  #if(mean(as.double(abundance[[i]][prIndex.hnf4a2,]))<1 && mean(as.double(abundance[[i]][prIndex.hnf4a,]))<1 && mean(as.double(abundanceGene[[i]][geneIndex.hnf4aas1,]))<1){next}
  print(sampleNames[i])
  # get all antisense RNA genes found in the gene biotype column in the Annotation files = "antisense"
  abundance.TCGAGTEx.asRNA = abundance[[i]][which(prAnnotationTable$gene_biotype=='antisense'),]
  abundance.TCGAGTEx.tumor.asRNA= abundance.tumor[[i]][which(prAnnotationTable$gene_biotype=='antisense'),]
  abundance.TCGAGTEx.not_tumor.asRNA = abundance.not_tumor[[i]][which(prAnnotationTable$gene_biotype=='antisense'),]
  abundanceGene.TCGAGTEx.asRNA = abundanceGene[[i]][which(ensemblAnnotationTxMatched$gene_biotype=='antisense'),]
  abundanceGene.TCGAGTEx.tumor.asRNA = abundanceGene.tumor[[i]][which(ensemblAnnotationTxMatched$gene_biotype=='antisense'),]
  abundanceGene.TCGAGTEx.not.tumor.asRNA = abundanceGene.not_tumor[[i]][which(ensemblAnnotationTxMatched$gene_biotype=='antisense'),]
  
  abundance.TCGAGTEx.asRNA = na.omit(abundance.TCGAGTEx.asRNA)
  abundanceGene.TCGAGTEx.asRNA = na.omit(abundanceGene.TCGAGTEx.asRNA)
  
  rownames(abundanceGene.TCGAGTEx.asRNA) = gsub('\\..*','',rownames(abundanceGene.TCGAGTEx.asRNA))
  
  #get all asRNA that target a gene with 2 or more promoters

  rownames(prAnnotationTable) <- prAnnotationTable$promoterId
  
  fraction=.1
  #for each asRNA
  
  
  numResults = nrow(results)
  numResults.tumor = nrow(results.tumor)
  numResults.not_tumor = nrow(results.not_tumor)
  
  hitst = 0
  hitsnt = 0
  hitsAll = 0
  hitst2 = 0
  hitsnt2 = 0
  hitsAll2 = 0
  candidatest = 0
  candidatesnt = 0
  candidatesAll = 0
  
  for (j in 1:length(rownames(abundanceGene.TCGAGTEx.asRNA))){
    pr = gsub('\\..*','',rownames(abundanceGene.TCGAGTEx.asRNA)[j])
    result=isCandidate(pr)
    if(is.na(result$p1) & is.na(result$p1t) & is.na(result$p1nt)){next}
    
    if(!is.na(result$p1) & !is.na(result$p2)){
      hitsAll=hitsAll+1
      hitIDs = c(hitIDs, pr)
      promoter1=result$p1
      promoter2=result$p2
      candidateResults = testCandidate(pr, abundanceGene, abundance,promoter1, promoter2, i)
      if(candidateResults$fullOutput[12]!=2){resultsAll = rbind(resultsAll, candidateResults$fullOutput)}
      hitsAll2=hitsAll2+1
      if(candidateResults$sig!=FALSE){
        candidatesAll=candidatesAll+1
        results = rbind(results,candidateResults$results)}
      }
    
    if(!is.na(result$p1t) & !is.na(result$p2t)){
      hitst=hitst+1
      hitIDst = c(hitIDst, pr)
      promoter1t=result$p1t
      promoter2t=result$p2t
      candidateResultst = testCandidate(pr, abundanceGene.tumor, abundance.tumor,promoter1t, promoter2t, i)
      if(candidateResultst$fullOutput[12]!=2){resultsAll.tumor = rbind(resultsAll.tumor, candidateResultst$fullOutput)}
      hitst2=hitst2+1
      if(candidateResultst$sig!=FALSE){
        candidatest=candidatest+1
        results.tumor = rbind(results.tumor, candidateResultst$results)
      }
    }
    
    if(!is.na(result$p1nt) & !is.na(result$p2nt)){
      hitsnt=hitsnt+1
      hitIDsnt = c(hitIDsnt, pr)
      promoter1nt=result$p1nt
      promoter2nt=result$p2nt
      candidateResultsnt = testCandidate(pr, abundanceGene.not_tumor, abundance.not_tumor, promoter1nt, promoter2nt, i)
      if(candidateResultsnt$fullOutput[12]!=2){resultsAll.not_tumor = rbind(resultsAll.not_tumor, candidateResultsnt$fullOutput)}
      hitsnt2=hitsnt2+1
      if(candidateResultsnt$sig!=FALSE){
        candidatesnt=candidatesnt+1
        results.not_tumor = rbind(results.not_tumor, candidateResultsnt$results)
      }
    }
    
    # if (candidateResults$sig!=FALSE | candidateResultsnt$sig!=FALSE | candidateResultst$sig!=FALSE){
    #   hitIDs = c(hitIDs, pr)
    # }
    # if (candidateResults$sig!=FALSE){
    #   hitIDs = c(hitIDs, pr)
    # }
    # if (candidateResultsnt$sig!=FALSE){
    #   hitIDsnt = c(hitIDsnt, pr)
    # }
    # if (candidateResultst$sig!=FALSE){
    #   hitIDst = c(hitIDst, pr)
    # }
    
    #resultsNum = nrow(results)
    #results.tumor.Num = nrow(results.tumor)
    #results.not_tumor.Num = nrow(results.not_tumor)
    
    #print(candidateResults$results[6])
    
    #if at least one is significant
    #if (nrow(results)-resultsNum == 1 | nrow(results.tumor)-results.tumor.Num == 1 | nrow(results.not_tumor)-results.not_tumor.Num == 1){
    #if (candidateResults$sig!=FALSE | candidateResultsnt$sig!=FALSE | candidateResultst$sig!=FALSE){
      
      #heatmapMatrix=rbind(heatmapMatrix, c(pr, sampleNames[i], candidateResults$results[6], candidateResultst$results[6], candidateResultsnt$results[6], candidateResults$results[7], candidateResultst$results[7], candidateResultsnt$results[7], candidateResults$results[8], candidateResultst$results[8], candidateResultsnt$results[8], candidateResults$results[9], candidateResultst$results[9], candidateResultsnt$results[9]))
      
      # png(paste("draft_figures/promoter_switch_candidates/",pr,"_",sampleNames[i],".png"),height = 1000, width = 1500)
      # par(mfrow=c(2,3))
      # boxplot(candidateResults$lowp1, candidateResults$highp1, candidateResults$lowp2, candidateResults$highp2,names=c("low AS (P1)", "high AS (P1)", "low AS (P2)", "high AS (P2)"), at = c(1,2,4,5), col = c("orange","red"), ylab = "(promoter expression)", main=paste(sampleNames[i]," all", as.character(candidateResults$sig)))
      # boxplot(candidateResultst$lowp1, candidateResultst$highp1, candidateResultst$lowp2, candidateResultst$highp2,names=c("low AS (P1)", "high AS (P1)", "low AS (P2)", "high AS (P2)"), at = c(1,2,4,5), col = c("orange","red"), ylab = "(promoter expression)", main=paste(sampleNames[i]," tumor", as.character(candidateResultst$sig)))
      # boxplot(candidateResultsnt$lowp1, candidateResultsnt$highp1, candidateResultsnt$lowp2, candidateResultsnt$highp2,names=c("low AS (P1)", "high AS (P1)", "low AS (P2)", "high AS (P2)"), at = c(1,2,4,5), col = c("orange","red"), ylab = "(promoter expression)", main=paste(sampleNames[i]," not tumor", as.character(candidateResultsnt$sig)))
      # plotCor(asRNA=as.double(log(abundanceGene.tumor[[i]][pr,])),p1data=as.double(abundance.tumor[[i]][promoter1,]),p2data=as.double(abundance.tumor[[i]][promoter2,]),asRNAnt=as.double(log(abundanceGene.not_tumor[[i]][pr,])),p1datant=as.double(abundance.not_tumor[[i]][promoter1,]),p2datant=as.double(abundance.not_tumor[[i]][promoter2,]))
      # dev.off()
    #}
  }
  print(hitst)
  print(hitsnt)
  print(hitsAll)
  print(candidatest)
  print(candidatesnt)
  print(candidatesAll)
  print(hitst2)
  print(hitsnt2)
  print(hitsAll2)
  # print(nrow(results)-numResults)
  # print(nrow(results.tumor)-numResults.tumor)
  # print(nrow(results.not_tumor)-numResults.not_tumor)
}

#cor(abundance[[i]][pr,],abundanceGene[[i]][promoter1,])
print(length(unique(hitIDs)))
print(length(unique(hitIDst)))
print(length(unique(hitIDsnt)))
results = results[2:nrow(results),]
results.tumor = results.tumor[2:nrow(results.tumor),]
results.not_tumor = results.not_tumor[2:nrow(results.not_tumor),]
resultsAll = resultsAll[2:nrow(resultsAll),]
resultsAll.tumor = resultsAll.tumor[2:nrow(resultsAll.tumor),]
resultsAll.not_tumor = resultsAll.not_tumor[2:nrow(resultsAll.not_tumor),]

write.csv(as.data.frame(results), file="promoter_switch_wilcoxon_all_gene.csv")
write.csv(as.data.frame(results.tumor), file="promoter_switch_wilcoxon_tumor_gene.csv")
write.csv(as.data.frame(results.not_tumor), file="promoter_switch_wilcoxon_not_tumor_gene.csv")
write.csv(as.data.frame(resultsAll), file="promoter_switch_wilcoxon_all_gene_all_hits.csv")
write.csv(as.data.frame(resultsAll.tumor), file="promoter_switch_wilcoxon_tumor_gene_all_hits.csv")
write.csv(as.data.frame(resultsAll.not_tumor), file="promoter_switch_wilcoxon_not_tumor_gene_all_hits.csv")

#violin plots of candidates

plotViolin <- function(results, abundance,abundanceGene,type){
for (j in 1:nrow(results)){
    p1_all=c()
    p2_all=c()
    as_all=c()
    samples=c()
    for (i in 1:length(abundance)){
      if(is.null(abundance[[i]])){next}
      p1_all=c(p1_all,as.double(abundance[[i]][results[j,3],]))
      p2_all=c(p2_all,as.double(abundance[[i]][results[j,4],]))
      as_all=c(as_all,log(as.double(abundanceGene[[i]][results[j,1],])+1))
      samples=c(samples,rep(sampleNames[[i]],length(abundance[[i]][results[j,3],])))
    }
    allname=c(rep.int("p1",length(p1_all)),rep.int("p2",length(p2_all)),rep.int("as",length(as_all)))
    samplescombined=c(samples,samples,samples)  
    all=c(p1_all,p2_all,as_all)
    df = data.frame(samplescombined, all, allname)
    ggplot(df, aes(x=samplescombined, y =all, fill=allname)) + 
      geom_violin(position=position_dodge(1), scale = "width", trim=FALSE) + 
      geom_boxplot(width=0.1, position=position_dodge(1)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    ggsave(paste("draft_figures/promoter_switch_candidates/",results[j,1],"_",results[j,3],'_',results[j,4],'_',type,"_violin.png"), width=16)
  }
}

plotViolin(results, abundance, abundanceGene,"all")
plotViolin(results.tumor, abundance.tumor, abundanceGene.tumor,"tumor")
plotViolin(results.not_tumor, abundance.not_tumor, abundanceGene.not_tumor,"not_tumor")

#heatmap

heatmapMatrix = heatmapMatrix[2:nrow(heatmapMatrix),]

temp=heatmapMatrix[,3:8]
mode(temp) <- 'numeric'
temp[is.na(temp)]=0
temp[is.infinite(temp)]=0
heatmap(temp)

source('R_scripts/heatmap.3.R')
#svg("Figure_1_gene_as_test.svg",width=12,height=4)
pdf("promoter_switch_heatmap.pdf",width=12,height=4)
heatmap.3(temp,
          scale = 'col',
          labRow = NA,
          #labCol = cognateGenes[,6],#cognateGenes[,6],#$unlist(geneSymbols),
          Rowv=NA,
          Colv=NA,
          RowSideColors=t(colorderv),
          symbreaks = FALSE,
          #breaks=colors,
          #col=brewer.pal(9,"Blues"),
)
dev.off()
print(nrow(results)/numCandidates)
print(nrow(results.tumor)/numCandidates)
print(nrow(results.not_tumor)/numCandidates)

results <- apply(results,2,as.character)
results.tumor <- apply(results.tumor,2,as.character)
results.not_tumor <- apply(results.not_tumor,2,as.character)


# check TGCT samples
# get the asRNA genes
sampleNames[28]
abundanceGene.asRNA = abundanceGene[[28]][which(ensemblAnnotationTxMatched$gene_biotype=='antisense'),]
abundanceGene.tumor.asRNA = abundanceGene.tumor[[28]][which(ensemblAnnotationTxMatched$gene_biotype=='antisense'),]
abundanceGene.not.tumor.asRNA = abundanceGene.not_tumor[[28]][which(ensemblAnnotationTxMatched$gene_biotype=='antisense'),]

rownames(abundanceGene.TCGAGTEx.asRNA) = gsub('\\..*','',rownames(abundanceGene.TCGAGTEx.asRNA))

# get average expression of the asRNA genes
mean(unlist(abundanceGene.asRNA))
mean(unlist(abundanceGene.tumor.asRNA))
mean(unlist(abundanceGene.not.tumor.asRNA))

for(i in 1:length(sampleNames)){
  print(sampleNames[[i]])
  #abundanceGene.asRNA = abundanceGene[[16]][which(ensemblAnnotationTxMatched$gene_biotype=='antisense'),]
  #abundanceGene.tumor.asRNA = abundanceGene.tumor[[16]][which(ensemblAnnotationTxMatched$gene_biotype=='antisense'),]
  #abundanceGene.not.tumor.asRNA = abundanceGene.not_tumor[[i]][which(ensemblAnnotationTxMatched$gene_biotype=='antisense'),]
  #mean(unlist(abundanceGene.asRNA))
  #mean(unlist(abundanceGene.tumor.asRNA))
  print(mean(unlist(abundanceGene[[i]])))
  #print(mean(unlist(abundanceGene.not.tumor.asRNA)))
}

# get average expression of the cognate genes
cognateGene = as.character(sgRNAPairs[unlist(asGene),3])
#get gene promoters of the cognate gene
cognatePromoters=prAnnotationTable[which(prAnnotationTable$ensembl_gene_id==cognateGene),]




#####compare all hits#####
#get all the gene/tissues that match across tumor/not tumor
studyNames=unique(resultsAll.tumor[,'Study'])
#tumor/nottumor plus/null/minus
pp=0
pn=0
pm=0
np=0
nn=0
nm=0
mp=0
mn=0
mm=0

for (i in 1:length(studyNames)){
  #add to count for if they are +/- null/- etc.... for correlation
  tumorTemp=resultsAll.tumor[which(resultsAll.tumor[,'Study']==studyNames[i]),]
  notTumorTemp=resultsAll.not_tumor[which(resultsAll.not_tumor[,'Study']==studyNames[i]),]
  rownames(tumorTemp)=tumorTemp[,1]
  rownames(notTumorTemp)=notTumorTemp[,1]
  for (j in 1:nrow(tumorTemp)){
    if(!(tumorTemp[j,1] %in% notTumorTemp[,1])){next}
    tp1=NA
    #tumorTemp[j,6:13]=as.numeric(tumorTemp[j,6:13])
    #notTumorTemp[tumorTemp[j,1],6:13]=as.numeric(notTumorTemp[tumorTemp[j,1],6:13])
    if (as.numeric(tumorTemp[j,10])<pThreshold & (as.numeric(tumorTemp[j,6])/as.numeric(tumorTemp[j,7])) > fcThreshold & as.numeric(tumorTemp[j,12]) > corThreshold){tp1="p"}
    else if (as.numeric(tumorTemp[j,10])<pThreshold & (as.numeric(tumorTemp[j,7])/as.numeric(tumorTemp[j,6])) > fcThreshold & as.numeric(tumorTemp[j,12]) < -corThreshold){tp1="m"}
    else{tp1="n"}
    tp2=NA
    if (as.numeric(tumorTemp[j,11])<pThreshold & (as.numeric(tumorTemp[j,8])/as.numeric(tumorTemp[j,9])) > fcThreshold & as.numeric(tumorTemp[j,13]) > corThreshold){tp2="p"}
    else if (as.numeric(tumorTemp[j,11])<pThreshold & (as.numeric(tumorTemp[j,8])/as.numeric(tumorTemp[j,9])) > fcThreshold & as.numeric(tumorTemp[j,13]) < -corThreshold){tp2="m"}
    else{tp2="n"}
    ntp1=NA
    if (as.numeric(notTumorTemp[tumorTemp[j,1],10])<pThreshold & (as.numeric(notTumorTemp[tumorTemp[j,1],6])/as.numeric(notTumorTemp[tumorTemp[j,1],7])) > fcThreshold & as.numeric(notTumorTemp[tumorTemp[j,1],12]) > corThreshold){ntp1="p"}
    else if (as.numeric(notTumorTemp[tumorTemp[j,1],10])<pThreshold & (as.numeric(notTumorTemp[tumorTemp[j,1],7])/as.numeric(notTumorTemp[tumorTemp[j,1],6])) > fcThreshold & as.numeric(notTumorTemp[tumorTemp[j,1],12]) < -corThreshold){ntp1="m"}
    else{ntp1="n"}
    ntp2=NA
    if (as.numeric(notTumorTemp[tumorTemp[j,1],11])<pThreshold & (as.numeric(notTumorTemp[tumorTemp[j,1],8])/as.numeric(notTumorTemp[tumorTemp[j,1],9])) > fcThreshold & as.numeric(notTumorTemp[tumorTemp[j,1],13]) > corThreshold){ntp2="p"}
    else if (as.numeric(notTumorTemp[tumorTemp[j,1],11])<pThreshold & (as.numeric(notTumorTemp[tumorTemp[j,1],8])/as.numeric(notTumorTemp[tumorTemp[j,1],9])) > fcThreshold & as.numeric(notTumorTemp[tumorTemp[j,1],13]) < -corThreshold){ntp2="m"}
    else{ntp2="n"}
    
    if(tp1=="p" & ntp1=="p" & tp2=="n" & ntp2=="n"){pp=pp+1}
    if(tp1=="p" & ntp1=="n" & tp2=="n" & ntp2=="n"){pn=pn+1}
    if(tp1=="p" & ntp1=="m" & tp2=="n" & ntp2=="n"){pm=pm+1}
    if(tp1=="n" & ntp1=="p" & tp2=="n" & ntp2=="n"){np=np+1}
    if(tp1=="n" & ntp1=="n" & tp2=="n" & ntp2=="n"){nn=nn+1}
    if(tp1=="n" & ntp1=="m" & tp2=="n" & ntp2=="n"){nm=nm+1}
    if(tp1=="m" & ntp1=="p" & tp2=="n" & ntp2=="n"){mp=mp+1}
    if(tp1=="m" & ntp1=="n" & tp2=="n" & ntp2=="n"){mn=mn+1}
    if(tp1=="m" & ntp1=="m" & tp2=="n" & ntp2=="n"){mm=mm+1}
    
    if(tp2=="p" & ntp2=="p" & tp1=="n" & ntp1=="n"){pp=pp+1}
    if(tp2=="p" & ntp2=="n" & tp1=="n" & ntp1=="n"){pn=pn+1}
    if(tp2=="p" & ntp2=="m" & tp1=="n" & ntp1=="n"){pm=pm+1}
    if(tp2=="n" & ntp2=="p" & tp1=="n" & ntp1=="n"){np=np+1}
    if(tp2=="n" & ntp2=="n" & tp1=="n" & ntp1=="n"){nn=nn+1}
    if(tp2=="n" & ntp2=="m" & tp1=="n" & ntp1=="n"){nm=nm+1}
    if(tp2=="m" & ntp2=="p" & tp1=="n" & ntp1=="n"){mp=mp+1}
    if(tp2=="m" & ntp2=="n" & tp1=="n" & ntp1=="n"){mn=mn+1}
    if(tp2=="m" & ntp2=="m" & tp1=="n" & ntp1=="n"){mm=mm+1}
    #if(p1p$p.value<pThreshold & p2p$p.value>pThreshold & (fc1 > fcThreshold | fc1 < (1/fcThreshold)) & (cor1 > corThreshold | cor1 < -corThreshold)){
      
  }
    #mean(highp1), mean(lowp1), mean(highp2), mean(lowp2), p1p$p.value, p2p$p.value, cor1, cor2)
}
print(pp)
print(pn)
print(pm)
print(np)
print(nn/2)
print(nm)
print(mp)
print(mn)
print(mm)

#this time look for pattern of association of ncNAT with both promoters

resultsAll.tumor=read.csv("promoter_switch_wilcoxon_tumor_gene_all_hits.csv")
resultsAll.not_tumor=read.csv("promoter_switch_wilcoxon_not_tumor_gene_all_hits.csv")
studyNames=unique(c(as.character(resultsAll.tumor[,"Study"]),as.character(resultsAll.not_tumor[,"Study"])))

results=matrix(ncol=length(studyNames)+1, nrow=12)
results[,1]=c("t_pp","t_pn","t_pm","t_nn","t_nm","t_mm","nt_pp","nt_pn","nt_pm","nt_nn","nt_nm","nt_mm")
for (i in 1:length(studyNames)){
  pp=0
  pn=0
  pm=0
  nn=0
  nm=0
  mm=0
  #add to count for if they are +/- null/- etc.... for correlation
  if(studyNames[i] %in% unique(resultsAll.tumor[,"Study"])){
    tumorTemp=resultsAll.tumor[which(resultsAll.tumor[,'Study']==studyNames[i]),]
    tumorTemp=tumorTemp[2:ncol(tumorTemp)]
    rownames(tumorTemp)=tumorTemp[,1]
    for (j in 1:nrow(tumorTemp)){
      p1=NA
      pval1=as.numeric(tumorTemp[j,10])
      fc1=as.numeric(tumorTemp[j,6])/as.numeric(tumorTemp[j,7])
      fc2=as.numeric(tumorTemp[j,7])/as.numeric(tumorTemp[j,6])
      cor1=as.numeric(tumorTemp[j,12])
      p2=NA
      pval2=as.numeric(tumorTemp[j,11])
      fc12=as.numeric(tumorTemp[j,8])/as.numeric(tumorTemp[j,9])
      fc22=as.numeric(tumorTemp[j,9])/as.numeric(tumorTemp[j,8])
      cor2=as.numeric(tumorTemp[j,13])
      
      #if (pval1<pThreshold & pval2>pThreshold & fc1 > fcThreshold & cor1 > corThreshold){p2="p"}
      #else if (pval1<pThreshold & pval2>pThreshold & fc2 > fcThreshold & cor1 < -corThreshold){p2="m"}
      if (pval1<pThreshold & fc1 > fcThreshold & cor1 > corThreshold){p2="p"}
      else if (pval1<pThreshold & fc2 > fcThreshold & cor1 < -corThreshold){p2="m"}
      else{p2="n"}
      
      if (pval2<pThreshold & fc12 > fcThreshold & cor2 > corThreshold){p1="p"}
      else if (pval2<pThreshold & fc22 > fcThreshold & cor2 < -corThreshold){p1="m"}
      else{p1="n"}
      
      if(p1 == 'p' & p2 == 'p'){pp=pp+1}
      if(p1 == 'n' & p2 == 'n'){nn=nn+1}
      if(p1 == 'm' & p2 == 'm'){mm=mm+1}
      if(p1 == 'p' & p2 == 'n'){pn=pn+1}
      if(p1 == 'n' & p2 == 'p'){pn=pn+1}
      if(p1 == 'p' & p2 == 'm'){pm=pm+1}
      if(p1 == 'm' & p2 == 'p'){pm=pm+1}
      if(p1 == 'n' & p2 == 'm'){nm=nm+1}
      if(p1 == 'm' & p2 == 'n'){nm=nm+1}
    }
  }

  ntpp=0
  ntpn=0
  ntpm=0
  ntnn=0
  ntnm=0
  ntmm=0
  if(studyNames[i] %in% unique(resultsAll.not_tumor[,"Study"])){
    nottumorTemp=resultsAll.not_tumor[which(resultsAll.not_tumor[,'Study']==studyNames[i]),]
    nottumorTemp=nottumorTemp[2:ncol(nottumorTemp)]
    rownames(nottumorTemp)=nottumorTemp[,1]
    for (j in 1:nrow(nottumorTemp)){
      p1=NA
      pval1=as.numeric(nottumorTemp[j,10])
      fc1=as.numeric(nottumorTemp[j,6])/as.numeric(nottumorTemp[j,7])
      fc2=as.numeric(nottumorTemp[j,7])/as.numeric(nottumorTemp[j,6])
      cor1=as.numeric(nottumorTemp[j,12])
      
      p2=NA
      pval2=as.numeric(nottumorTemp[j,11])
      fc12=as.numeric(nottumorTemp[j,8])/as.numeric(nottumorTemp[j,9])
      fc22=as.numeric(nottumorTemp[j,9])/as.numeric(nottumorTemp[j,8])
      cor2=as.numeric(nottumorTemp[j,13])

      #if (pval1<pThreshold & pval2>pThreshold & fc1 > fcThreshold & cor1 > corThreshold){p2="p"}
      #else if (pval1<pThreshold & pval2>pThreshold & fc2 > fcThreshold & cor1 < -corThreshold){p2="m"}
      if (pval1<pThreshold & fc1 > fcThreshold & cor1 > corThreshold){p2="p"}
      else if (pval1<pThreshold & fc2 > fcThreshold & cor1 < -corThreshold){p2="m"}
      else{p2="n"}
      
      #if (pval2<pThreshold & pval1>pThreshold &  fc12 > fcThreshold & cor2 > corThreshold){p1="p"}
      #else if (pval2<pThreshold & pval1>pThreshold & fc22 > fcThreshold & cor2 < -corThreshold){p1="m"}
      if (pval2<pThreshold &  fc12 > fcThreshold & cor2 > corThreshold){p1="p"}
      else if (pval2<pThreshold & fc22 > fcThreshold & cor2 < -corThreshold){p1="m"}
      else{p1="n"}
      
      if(p1 == 'p' & p2 == 'p'){ntpp=ntpp+1}
      if(p1 == 'n' & p2 == 'n'){ntnn=ntnn+1}
      if(p1 == 'm' & p2 == 'm'){ntmm=ntmm+1}
      if(p1 == 'p' & p2 == 'n'){ntpn=ntpn+1}
      if(p1 == 'n' & p2 == 'p'){ntpn=ntpn+1}
      if(p1 == 'p' & p2 == 'm'){ntpm=ntpm+1}
      if(p1 == 'm' & p2 == 'p'){ntpm=ntpm+1}
      if(p1 == 'n' & p2 == 'm'){ntnm=ntnm+1}
      if(p1 == 'm' & p2 == 'n'){ntnm=ntnm+1}
    }
  }
  
  results[,i+1]=c(pp,pn,pm,nn,nm,mm,ntpp,ntpn,ntpm,ntnn,ntnm,ntmm)
}
colnames(results)=c("behavior",studyNames)
write.csv(as.data.frame(results), file="two_promoter_behaviors_individual.csv")