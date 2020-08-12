#run promoter_Switch_wilcoxon_gene3.R

#histogram of number of promoters per gene
p1=hist(table(prAnnotationTable$ensembl_gene_id), breaks=30)

#histogram of number of promoters per ncNAT regulated gene
p2=hist(table(prAnnotationTable[which(prAnnotationTable$ensembl_gene_id %in% sgRNAPairs$cognate_gene),]$ensembl_gene_id), breaks=30)

plot(p1, col=rgb(0,0,1,1/4), freq=FALSE)
lines(p2, col=rgb(1,0,0,1/4), freq=FALSE)

#histogram of length of ncNAT compared to those that regulate genes of two or more promoters
p1=hist(ensemblAnnotationTxMatched[unique(c(results.tumor[,1],results.not_tumor[,1])),7]-ensemblAnnotationTxMatched[unique(c(results.tumor[,1],results.not_tumor[,1])),6])
p2=hist(ensemblAnnotationTxMatched[which(ensemblAnnotationTxMatched$gene_biotype=="antisense"),7]-ensemblAnnotationTxMatched[which(ensemblAnnotationTxMatched$gene_biotype=="antisense"),6])
plot(p1, col=rgb(0,0,1,1/4), freq=FALSE)
lines(p2, col=rgb(1,0,0,1/4), freq=FALSE)

#histogram of exon number of ncNAT compared to those that regulate genes of two or more promoters

#histogram of expression of ncNAT compared to those that regulate genes of two or more promoters

#distance of promoter switching asRNA to promoter

#barplot of tissue count of candidates
results.tumor=read.csv("promoter_switch_wilcoxon_tumor_gene.csv")
results.tumor=results.tumor[2:ncol(results.tumor)]
results.not_tumor=read.csv("promoter_switch_wilcoxon_not_tumor_gene.csv")
results.not_tumor=results.not_tumor[2:ncol(results.not_tumor)]
resultsAll.tumor=read.csv("promoter_switch_wilcoxon_tumor_gene_all_hits.csv")
resultsAll.tumor= resultsAll.tumor[2:ncol(resultsAll.tumor)]
resultsAll.not_tumor=read.csv("promoter_switch_wilcoxon_not_tumor_gene_all_hits.csv")
resultsAll.not_tumor=resultsAll.not_tumor[2:ncol(resultsAll.not_tumor)]

tumorTable=table(resultsAll.tumor[,2])
not_tumorTable=table(resultsAll.not_tumor[,2])
xnames=names(table(resultsAll.not_tumor[,2]))
for (i in 1:length(xnames)){
  if(!(xnames[i] %in% names(tumorTable))){
    tumorTable[xnames[i]]=0
  }
}
xnames=names(table(resultsAll.tumor[,2]))
for (i in 1:length(xnames)){
  if(!(xnames[i] %in% names(not_tumorTable))){
    not_tumorTable[xnames[i]]=0
  }
}
tumorTable=tumorTable[sort(names(tumorTable))]
not_tumorTable=not_tumorTable[sort(names(not_tumorTable))]

svg("hits_per_sample.svg",width=6,height=4)
#barplot(rbind(table(resultsAll.tumor[,2]),table(resultsAll.not_tumor[,2])), beside=TRUE, col=c("steelblue2","lightcoral"), legend=c("tumor","not tumor"), las=2, ylab="Number of tested genes", xlab="Samples")
barplot(rbind(tumorTable,not_tumorTable), beside=TRUE, col=c("steelblue2","lightcoral"), legend=c("tumor","not tumor"), las=2, ylab="Number of tested genes", xlab="Tumor Type", ylim=c(0,300))
dev.off()

tumorTable=table(results.tumor[,2])
not_tumorTable=table(results.not_tumor[,2])
xnames=names(table(results.not_tumor[,2]))
for (i in 1:length(xnames)){
  if(!(xnames[i] %in% names(tumorTable))){
    tumorTable[xnames[i]]=0
  }
}
xnames=names(table(results.tumor[,2]))
for (i in 1:length(xnames)){
  if(!(xnames[i] %in% names(not_tumorTable))){
    not_tumorTable[xnames[i]]=0
  }
}
tumorTable=tumorTable[sort(names(tumorTable))]
not_tumorTable=not_tumorTable[sort(names(not_tumorTable))]

svg("hits_candidates_per_sample.svg",width=6,height=4)
barplot(rbind(tumorTable,not_tumorTable), beside=TRUE, col=c("steelblue2","lightcoral"), legend=c("tumor","not tumor"), las=2, ylab="Number of candidates", xlab="Tumor Type", ylim=c(0,30))
dev.off()

tumorTable=table(results.tumor[,2])
not_tumorTable=table(results.not_tumor[,2])
xnames=names(table(resultsAll.tumor[,2]))
for (i in 1:length(xnames)){
  if(!(xnames[i] %in% names(tumorTable))){
    tumorTable[xnames[i]]=0
  }
}

xnames=names(table(resultsAll.not_tumor[,2]))
for (i in 1:length(xnames)){
  if(!(xnames[i] %in% names(not_tumorTable))){
    not_tumorTable[xnames[i]]=0
  }
}

tumorTable=tumorTable[sort(names(tumorTable))]
not_tumorTable=not_tumorTable[sort(names(not_tumorTable))]

tumorTable=tumorTable/table(resultsAll.tumor[,2])
not_tumorTable=not_tumorTable/table(resultsAll.not_tumor[,2])

xnames=c(names(tumorTable),names(not_tumorTable))
for (i in 1:length(xnames)){
  if(!(xnames[i] %in% names(not_tumorTable))){
    not_tumorTable[xnames[i]]=0
  }
  if(!(xnames[i] %in% names(tumorTable))){
    tumorTable[xnames[i]]=0
  }
}
tumorTable=tumorTable[sort(names(tumorTable))]
not_tumorTable=not_tumorTable[sort(names(not_tumorTable))]
svg("hit_to_candidate_proportion_per_sample.svg",width=6,height=4)  
barplot(rbind(tumorTable,not_tumorTable), beside=TRUE, col=c("steelblue2","lightcoral"), legend=c("tumor","not tumor"), las=2, ylab="Proportion of hits/candidates", xlab="Tumor Type", ylim=c(0,0.12))
dev.off()

#number of hits per sample and candidate conversion
hitCount=read.csv('promoter_switch_candidate_count.csv')
barplot(rbind(unlist(hitCount[1,]),unlist(hitCount[2,])), beside=TRUE, col=c("steelblue2","lightcoral"), legend=c("tumor","not tumor"), xlab = "sample name", ylab = "Number of tested asRNA targets", las=2)
barplot(rbind(unlist(hitCount[7,]),unlist(hitCount[8,])), beside=TRUE, col=c("steelblue2","lightcoral"), legend=c("tumor","not tumor"), xlab = "sample name", ylab = "Significant candidates", las=2)
#percentage based
barplot(rbind(unlist(hitCount[3,])/unlist(hitCount[1,]),unlist(hitCount[4,])/unlist(hitCount[2,])), beside=TRUE, col=c("steelblue2","lightcoral"), legend=c("tumor","not tumor"), las=2, xlab = "Sample Name", ylab ="")
barplot(rbind(unlist(hitCount[7,])/unlist(hitCount[5,]),unlist(hitCount[8,])/unlist(hitCount[6,])), beside=TRUE, col=c("steelblue2","lightcoral"), legend=c("tumor","not tumor"), las=2, xlab = "Sample Name", ylab ="Significant Proportion of Promoter Switch Candidates")


df <- data.frame(A = sample(1:5, 30, T), 
                 B = sample(c('Y', 'N'), 30, T), 
                 C = rep(LETTERS[1:3], 10))
ggplot(df) + geom_bar(aes(B, fill = C), position  = 'stack', width = 0.9) + 
  facet_wrap(~A, nrow = 1) + theme(panel.spacing = unit(0, "lines"))
#hitCount[2,]-hitCount[4]
#hitCount[3,]-hitCount[5]


#histogram of number of tissues for each candidate
#tumor
p2=hist(table(results.tumor[,1]), breaks=20)
#not_tumor
p1=hist(table(results.not_tumor[,1]), breaks=20)
plot( p1, col=rgb(0,0,1,1/4)) 
plot( p2, col=rgb(1,0,0,1/4), add=T)  

#how many candidates are driver genes
drivers=read.csv('TableS1_compendium_mutational_drivers_cdsA_only.csv')
candidatest=sgRNAPairs[results.tumor[,1],"cognate_gene"]
cognatematches=candidatest[which(candidatest %in% gsub('\\..*','',drivers[,3]))]
asmatches=sgRNAPairs[which(sgRNAPairs[,"cognate_gene"] %in% cognatematches),"asRNA"]

candidatesnt=sgRNAPairs[results.not_tumor[,1],"cognate_gene"]
cognatematches=candidatesnt[which(candidatesnt %in% gsub('\\..*','',drivers[,3]))]
asmatches=sgRNAPairs[which(sgRNAPairs[,"cognate_gene"] %in% cognatematches),"asRNA"]
