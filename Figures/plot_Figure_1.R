#load in data
resOrdered = read.csv("ncNAT_DE.csv",row.names = c(1))

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

#remove batch effects between the TCGA and liver data
abundance.TCGAGTEx.nobatch=removeBatchEffect(abundance.TCGAGTEx, batch = coldata[,2], design=model.matrix(~coldata[,1]))

#plot as heatmap
candidates=resOrdered[na.omit(resOrdered)$padj<0.01,]
if(nrow(candidates)>100){
  candidates=candidates[1:100,]
}
candidates <- candidates[order(candidates$log2FoldChange,decreasing = TRUE),]
candidates = na.omit(candidates)
candidates=candidates[abs(candidates$log2FoldChange)>=1,]

#extract the rows from abudance tables 
data=as.matrix(abundance.TCGAGTEx.nobatch[rownames(candidates),])
colorder = c(colnames(abundance.TCGAGTEx.not_tumor),colnames(abundance.TCGAGTEx.tumor))
datax=data[,colorder]
heatmap(datax, Colv=NA, Rowv=NA, col=brewer.pal(9,"Reds"), scale = 'row', keep.dendro = FALSE, xlab = "Tumor - Non-tumor", ylab = "promoter expression", labRow=NA, labCol=NA) #blues
datat = t(datax)

getSidebars = function(){
  cognateGenes= matrix(nrow=nrow(candidates),ncol =8)
  colnames(cognateGenes)<-c("cognate","expression", 'correlation', 'cancerCensus', 'cancerDriver','cognate_symbol','expression2', 'correlation2')
  
  sgRNAPairs=read.csv("sgRNA_cognate_pairs.csv")
  rownames(sgRNAPairs)<-sgRNAPairs[,2]
  #match the promoter ids to gene id and then the partner gene
  if (isPromoter){
    rownames(prAnnotationTable) <- prAnnotationTable$promoterId
    promoterGenes=prAnnotationTable[rownames(candidates),9]
  } else {
    promoterGenes=rownames(candidates)
    promoterGenes=gsub('\\..*','',promoterGenes)
  }
  
  cognateGenes[,1] = as.character(sgRNAPairs[unlist(promoterGenes),3])
  
  #get the liver expression of the partner genes
  #promoterGenes_RE=gsub('\\..*','',promoterGenes)
  cognateGenes[,1] = as.character(sgRNAPairs[unlist(promoterGenes),3])
  cognateGenes[,6] = ensemblAnnotationTxMatched[cognateGenes[,1],2]
  
  #average the expression across non-cancer liver cells
  rownames(abundanceGene.TCGAGTEx.not_tumor) <- gsub('\\..*','',rownames(abundanceGene.TCGAGTEx.not_tumor))
  rownames(cognateGenes)<-rownames(candidates)
  for (i in 1:nrow(cognateGenes)){
    expression = as.numeric(mean(as.matrix(abundanceGene.TCGAGTEx.not_tumor[cognateGenes[i,1],])))
    if (is.na(expression)){
      cognateGenes[i,2]= 'snow2'
      cognateGenes[i,7]=0
      next
    }
    cognateGenes[i,7]=expression
    if (expression > 20){
      cognateGenes[i,2]= 'slategrey'
    }
    else {
      cognateGenes[i,2]= 'snow2'
    }
    #cognateGenes[i,2]= expression
  }
  
  #calculate correlation between asRNA and cognate gene
  correlations=c()
  correlationsP=c()
  for (i in 1:nrow(cognateGenes)){
    asRNAexpression=as.double(abundanceGene.TCGAGTEx[unlist(promoterGenes)[i],])
    cognateGeneExpression=as.double(abundanceGene.TCGAGTEx[cognateGenes[i,1],])
    correlation=as.double(cor(asRNAexpression,cognateGeneExpression,method = "pearson"))
    
    if (!is.na(correlation)){
      correlations=c(correlations, correlation)
      correlationsP=c(correlationsP, cor.test(asRNAexpression,cognateGeneExpression,method = "pearson")$p.value)
    }
    
    colour="white"
    cognateGenes[i,8]=correlation
    if (is.na(correlation)){
      colour="white"
      cognateGenes[i,8]=0
    }
    else if (correlation > .75){colour="royalblue4"}
    else if (correlation > .50){colour="royalblue1"}
    else if (correlation > .25){colour="skyblue1"}
    else if (correlation > -.25){colour="white"}
    else if (correlation > -.50){colour="salmon1"}
    else if (correlation > -.75){colour="tomato2"}
    else if (correlation > -1.1){colour="tomato4"}
    cognateGenes[i,3]= colour
  }
  
  #is it a cancer census gene
  canceroutput=read.csv("sgRNA_cancer_census.csv")
  rownames(canceroutput) = canceroutput[,1]
  for (i in 1:length(promoterGenes)){
    cognateGenes[i,4] = "snow2"
    if (any(rownames(canceroutput)==unlist(promoterGenes)[i]))
    {
      if (!is.na(canceroutput[unlist(promoterGenes)[i],3]))
      {
        if (grepl('L',canceroutput[unlist(promoterGenes)[i],4])){
          cognateGenes[i,4] = "lightsalmon1"
        } else {
          cognateGenes[i,4] = "slategrey"
        }
        
      }
      #cognateGenes[i,4]= canceroutput[unlist(promoterGenes[i,1]),3]
    }
  }
  
  #is it a cancer driver gene
  drivergenes=read.csv("TableS1_compendium_mutational_drivers_cdsA_only.csv")
  rownames(drivergenes)<-drivergenes[,2]
  for (i in 1:nrow(cognateGenes)){
    symbol=ensemblAnnotationTxMatched[cognateGenes[i,1],2]
    if (is.na(symbol)){cognateGenes[i,5] = "snow2"}
    else {
      if (any(drivergenes==symbol)){
        if (grepl('[Ll]iver',drivergenes[symbol,7])){
          cognateGenes[i,5] = "lightsalmon1"
        }
        else{
          cognateGenes[i,5] = "slategrey"
        }
      }
      else {
        cognateGenes[i,5] = "snow2"
      }
    }
  }
  return(cognateGenes)
}

cognateGenes=getSidebars()
sidebars = cbind(cognateGenes[,2],cognateGenes[,3],cognateGenes[,4],cognateGenes[,5])
colnames(sidebars)<-c('liver expressed cognate gene','correlation','census','driver')
colorderv=c(rep('lightgreen', times=ncol(abundance.TCGAGTEx.not_tumor)),rep('lightsalmon1',times=ncol(abundance.TCGAGTEx.tumor)))
colors = seq(-1,2,length.out=10)
if (isPromoter){
  geneSymbols=prAnnotationTable[rownames(candidates),10]
} else {
  geneSymbols=ensemblAnnotation[gsub('\\..*','',rownames(candidates)),2]
  geneSymbols=unlist(geneSymbols)
  ids=unlist(ensemblAnnotation[gsub('\\..*','',rownames(candidates)),1])
  for (i in 1:length(geneSymbols)){
    if(is.na(geneSymbols[i]) || geneSymbols[i]==""){
      geneSymbols[i]=ids[i]
    }
  }
}

write.csv(cognateGenes,file="figure_1_sidebars.csv")

source('R_scripts/heatmap.3.R')
#svg("Figure_1_gene_as_test.svg",width=12,height=4)
pdf("Figure_1_gene_as_test.pdf",width=12,height=4)
heatmap.3(datat,
          scale = 'col',
          labRow = NA,
          labCol = cognateGenes[,6],#cognateGenes[,6],#$unlist(geneSymbols),
          Rowv=NA,
          Colv=NA,
          RowSideColors=t(colorderv),
          ColSideColors=sidebars,
          symbreaks = FALSE,
          breaks=colors,
          col=brewer.pal(9,"Blues"),
)
dev.off()

expression=cbind(as.numeric(cognateGenes[,7]),as.numeric(cognateGenes[,7]))
svg("sidebar1_expression.svg",width=12,height=4)
heatmap.3(expression,
          scale = 'col',
          labRow = NA,
          Rowv=NA,
          Colv=NA,
          symbreaks = FALSE,
          breaks=colors,
          col=brewer.pal(9,"Blues"),
)
dev.off()

colors = seq(-1.1,1.1,length.out=10)
correlation=cbind(as.numeric(cognateGenes[,8]),as.numeric(cognateGenes[,8]))
svg("sidebar1_correlation.svg",width=12,height=4)
heatmap.3(correlation,
          scale = 'none',
          labRow = NA,
          Rowv=NA,
          Colv=NA,
          symbreaks = FALSE,
          breaks=colors,
          col=brewer.pal(9,"RdBu"),
)
dev.off()



