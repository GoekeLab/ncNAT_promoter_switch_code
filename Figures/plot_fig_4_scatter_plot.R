gene ='ENSG00000261105'
p1='prmtr.30557'
p2='prmtr.30554'
geneID <- which(ensemblAnnotationTxMatched$ensembl_gene_id==gene)
p2ID <- which(prAnnotationTable$promoterId==p1)
p1ID <- which(prAnnotationTable$promoterId==p2)
sample = "KIRC"

for (i in 1:length(sampleNames)){
  if(sampleNames[i]==sample){break}
}

abundance.TCGAGTEx=abundance[[i]]
abundance.TCGAGTEx.tumor=abundance.tumor[[i]]
abundance.TCGAGTEx.not_tumor=abundance.not_tumor[[i]]
abundanceGene.TCGAGTEx=abundanceGene[[i]]
abundanceGene.TCGAGTEx.tumor=abundanceGene.tumor[[i]]
abundanceGene.TCGAGTEx.not_tumor=abundanceGene.not_tumor[[i]]

#png(paste("draft_figures/cancer_vs_normal_gene/",sample1,"_",sample2,"_HNF4A_AS_gene.png",sep=""), width=1200, height=1200)
svg(paste("draft_figures/figure_x/",gene,"_",sample,".svg",sep=""), width=12, height=12)
par(mfrow=c(2,2))

#plot cancer vs normal for each cancer type showing relationship between AS and P1 
plot(as.double(log(abundanceGene.TCGAGTEx.tumor[geneID,])), as.double(abundance.TCGAGTEx.tumor[p1ID,]), xlab=paste(gene,' (log) gene expression',sep = ""), ylab=paste(p1,' promoter activity',sep = ""), col=c("lightcoral"), main=sample, pch=16, cex=1.3)
points(as.double(log(abundanceGene.TCGAGTEx.not_tumor[geneID,])), as.double(abundance.TCGAGTEx.not_tumor[p1ID,]), col=c("steelblue2"), pch=16, cex=1.3)

x=as.data.frame(t(rbind(as.double(log(abundanceGene.TCGAGTEx.tumor[geneID,])), as.double(abundance.TCGAGTEx.tumor[p1ID,]))))
colnames(x) = c('asRNA','gene')
x <- x[!is.infinite(rowSums(x)),]
reg=lm(formula = gene~asRNA, data=x)
abline(reg, col=c("lightcoral"), lwd=3)
text(par('usr')[1]+2,y=par('usr')[3]+2,labels=as.character(signif(summary(reg)$r.squared,3)), col='lightcoral')

x=as.data.frame(t(rbind(as.double(log(abundanceGene.TCGAGTEx.not_tumor[geneID,])), as.double(abundance.TCGAGTEx.not_tumor[p1ID,]))))
colnames(x) = c('asRNA','gene')
x <- x[!is.infinite(rowSums(x)),]
reg=lm(formula = gene~asRNA, data=x)
abline(reg, col=c("steelblue2"), lwd=3)
text(par('usr')[1]+2,y=par('usr')[3]+3,labels=as.character(signif(summary(reg)$r.squared,3)), col='steelblue2')

legend("topleft", legend=c("Tumor","Non-tumor"),fill=c("lightcoral", "steelblue2"), bty='n')

correlation=cor(as.double(cbind(log(abundanceGene.TCGAGTEx.tumor[geneID,]+1),log(abundanceGene.TCGAGTEx.not_tumor[geneID,]+1))),as.double(cbind(abundance.TCGAGTEx.tumor[p1ID,],abundance.TCGAGTEx.not_tumor[p1ID,])), method = "pearson")
text(par('usr')[1]+1,y=par('usr')[3]+3,labels=as.character(signif(correlation,3)), col='black')
correlation=cor(as.double(log(abundanceGene.TCGAGTEx.tumor[geneID,]+1)),as.double(abundance.TCGAGTEx.tumor[p1ID,]), method = "pearson")
text(par('usr')[1]+1,y=par('usr')[3]+2,labels=as.character(signif(correlation,3)), col='lightcoral')
correlation=cor(as.double(log(abundanceGene.TCGAGTEx.not_tumor[geneID,]+1)),as.double(abundance.TCGAGTEx.not_tumor[p1ID,]), method = "pearson")
text(par('usr')[1]+1,y=par('usr')[3]+1,labels=as.character(signif(correlation,3)), col='steelblue2')


#plot cancer vs normal for each cancer type showing relationship between AS and P2
plot(log(as.double(abundanceGene.TCGAGTEx.tumor[geneID,])), as.double(abundance.TCGAGTEx.tumor[p2ID,]), xlab=paste(gene,' (log) gene expression',sep = ""), ylab=paste(p2,' promoter activity',sep = ""), col=c("lightcoral"), pch=16, cex=1.3)
points(log(as.double(abundanceGene.TCGAGTEx.not_tumor[geneID,])), as.double(abundance.TCGAGTEx.not_tumor[p2ID,]), col=c("steelblue2"), pch=16, cex=1.3)

x=as.data.frame(t(rbind(as.double(log(abundanceGene.TCGAGTEx.tumor[geneID,])), as.double(abundance.TCGAGTEx.tumor[p2ID,]))))
colnames(x) = c('asRNA','gene')
x <- x[!is.infinite(rowSums(x)),]
reg=lm(formula = gene~asRNA, data=x)
abline(reg, col=c("lightcoral"), lwd=3)
text(par('usr')[1]+2,y=par('usr')[3]+2,labels=as.character(signif(summary(reg)$r.squared,3)), col='lightcoral')

x=as.data.frame(t(rbind(as.double(log(abundanceGene.TCGAGTEx.not_tumor[geneID,])), as.double(abundance.TCGAGTEx.not_tumor[p2ID,]))))
colnames(x) = c('asRNA','gene')
x <- x[!is.infinite(rowSums(x)),]
reg=lm(formula = gene~asRNA, data=x)
abline(reg, col=c("steelblue2"), lwd=3)
text(par('usr')[1]+2,y=par('usr')[3]+3,labels=as.character(signif(summary(reg)$r.squared,3)), col='steelblue2')

correlation=cor(as.double(cbind(log(abundanceGene.TCGAGTEx.tumor[geneID,]+1),log(abundanceGene.TCGAGTEx.not_tumor[geneID,]+1))),as.double(cbind(abundance.TCGAGTEx.tumor[p2ID,],abundance.TCGAGTEx.not_tumor[p2ID,])), method = "pearson")
text(par('usr')[1]+1,y=par('usr')[3]+3,labels=as.character(signif(correlation,3)), col='black')
correlation=cor(as.double(log(abundanceGene.TCGAGTEx.tumor[geneID,]+1)),as.double(abundance.TCGAGTEx.tumor[p2ID,]), method = "pearson")
text(par('usr')[1]+1,y=par('usr')[3]+2,labels=as.character(signif(correlation,3)), col='lightcoral')
correlation=cor(as.double(log(abundanceGene.TCGAGTEx.not_tumor[geneID,]+1)),as.double(abundance.TCGAGTEx.not_tumor[p2ID,]), method = "pearson")
text(par('usr')[1]+1,y=par('usr')[3]+1,labels=as.character(signif(correlation,3)), col='steelblue2')

#plot cancer vs normal for each cancer type comparing p1 and p2
plot(as.double(abundance.TCGAGTEx.tumor[p1ID,]), as.double(abundance.TCGAGTEx.tumor[p2ID,]), xlab=paste(p1,' promoter activity',sep = ""), ylab=paste(p2,' promoter activity',sep = ""), col=c("lightcoral"), pch=16, cex=1.3)
points(as.double(abundance.TCGAGTEx.not_tumor[p1ID,]), as.double(abundance.TCGAGTEx.not_tumor[p2ID,]), col=c("steelblue2"), pch=16, cex=1.3)

correlation=cor(as.double(cbind(abundance.TCGAGTEx.tumor[p1ID,],abundance.TCGAGTEx.not_tumor[p1ID,])),as.double(cbind(abundance.TCGAGTEx.tumor[p2ID,],abundance.TCGAGTEx.not_tumor[p2ID,])), method = "pearson")
text(par('usr')[1]+1,y=par('usr')[3]+3,labels=as.character(signif(correlation,3)), col='black')
correlation=cor(as.double(abundance.TCGAGTEx.tumor[p1ID,]),as.double(abundance.TCGAGTEx.tumor[p2ID,]), method = "pearson")
text(par('usr')[1]+1,y=par('usr')[3]+2,labels=as.character(signif(correlation,3)), col='lightcoral')
correlation=cor(as.double(abundance.TCGAGTEx.not_tumor[p1ID,]),as.double(abundance.TCGAGTEx.not_tumor[p2ID,]), method = "pearson")
text(par('usr')[1]+1,y=par('usr')[3]+1,labels=as.character(signif(correlation,3)), col='steelblue2')

#plot cancer vs normal for each cancer type comparing AS and the P1-P2 ratio
p1t=as.double(abundance.TCGAGTEx.tumor[p1ID,])
p2t=as.double(abundance.TCGAGTEx.tumor[p2ID,])
plot(log(as.double(abundanceGene.TCGAGTEx.tumor[geneID,])),p1t/(p1t+p2t), xlab=paste(gene,' (log) gene expression',sep = ""), ylab=paste(p1,' expression proportion',sep = ""), col=c("lightcoral"), pch=16, cex=1.3)
p1nt=as.double(abundance.TCGAGTEx.not_tumor[p1ID,])
p2nt=as.double(abundance.TCGAGTEx.not_tumor[p2ID,])
points(log(as.double(abundanceGene.TCGAGTEx.not_tumor[geneID,])),p1nt/(p1nt+p2nt), col=c("steelblue2"), pch=16, cex=1.3)

dev.off()