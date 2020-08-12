#plot promoter switch candidates

#need to get p value and correlation for each candidate

#for each candidate get the required information
library(dplyr)
library(ggplot2)
library(ggpubr)

getPlotData <- function(candidates, abundance, abundanceGene){
  x=c()
  y=c()
  size=c()
  color=c()
  
  for (i in 1:nrow(candidates)){
    pr = candidates[i,1]
    promoter1 = candidates[i,3]
    promoter2 = candidates[i,4]
    rowi=0
    for (j in 1:length(abundance)){
      #for (i in 16:16){
      if(is.null(abundance[[j]])){next}
      rowi=rowi+1
      candidateResultst = testCandidate(pr, abundanceGene, abundance, promoter1, promoter2, j)
      
      if (candidateResultst$sig != FALSE){
        x=c(x, as.character(i))
        y=c(y, as.character(rowi))
        size = c(size,as.numeric(candidateResultst$results[6]))
        if(candidateResultst$results[6]=="NaN"){
          color=c(color,'white')
        } else if (as.numeric(candidateResultst$results[6])>1){
          color=c(color,'blue')
        } else {
          color=c(color,'red')
        }
      }
    }
  }
  size=abs(log(size,2))
  size[size>2]=2
  size[size<-2]=2
  return(list(x=x, y=y,size=size,color=color))
}

getPlotDataSimple <- function(candidates){
  x=c()
  x2=c()
  y=c()
  size=c()
  color=c()
  ylab=c()
  for (i in 1:nrow(candidates)){
    pr = candidates[i,1]
    promoter1 = candidates[i,3]
    promoter2 = candidates[i,4]
    sample=candidates[i,2]
    rowi=which(unique(candidates[,'tissue'])==sample)
    # rowi=0
    # for (j in 1:length(unique(candidates[,'tissue']))){
    #   #for (i in 16:16){
    #   if (sampleNames[[j]]==sample){break}
    # }
    # rowi=j
    #x=c(x, as.character(i))
    x=c(x, as.character(candidates[i,5]))
    x2=c(x2, as.character(candidates[i,1]))
    #if(rowi<10){y=c(y, paste('0',as.character(rowi),sep=""))}
    #else{y=c(y, as.character(rowi))}
    y=c(y, as.character(unlist(sample)))
    
    #if(!(sampleNames[[j]] %in% ylab)){ylab=c(ylab,sampleNames[[j]])}
    size = c(size,as.numeric(candidates[i,6]))
    if (as.numeric(candidates[i,6])>1){
      color=c(color,'blue')
    } else {
      color=c(color,'red')
    }
  }
  size=abs(log(size,2))
  size[size>2]=2
  size[size<-2]=2
  return(list(x=x, x2=x2, y=y,size=size,color=color, ylab=ylab))
}

plotCandidates <- function(x, y, size, color, condition){
  #create size legend
  legend_bubbles <- data.frame(
    label = c("3", "20", "40m"),
    size  = c(2,1,0.5)
  ) %>%
    mutate(radius = sqrt(size / pi))   
  
  x=factor(x, levels=unique(x))
  p <- qplot(x, y, geom="point", size=size, color=color, 
             main=paste("Promoter Switching Candidates in ", condition, " samples")) + 
        theme(legend.position = "none")
  p <- p + theme( axis.line = element_line(colour = "darkblue", 
                                    size = 1, linetype = "solid"))
  #p <- p + scale_x_discrete(breaks=c(1:nrow(results.tumor)),
  #                     labels=results.tumor[,5])
  #p <- p + scale_y_discrete(breaks=c(1:length(ylab)),labels=ylab)
  
  p <- p + geom_point(data = legend_bubbles,
                      #  The "radius/50" was trial and error. Better way?
                      aes(x = 5, y = 4 + radius, size = size),
                      shape = 21, color = "black", fill = c('red','blue','red')) +
    annotate("text", x = 4, y = 8, label = "Log 2 Fold Change", fontface = "bold") +
    annotate("text", x = 4, y = 4, label = as.character(signif(min(size),3)), fontface = "bold") + 
    annotate("text", x = 4, y = 5, label = paste('-',as.character(signif(((max(size)-min(size))/2)+min(size),3))), fontface = "bold") +
    annotate("text", x = 4, y = 6, label = paste('>',as.character(signif(max(size)),3)), fontface = "bold") +
    labs(x = "Regulated Genes") +
    labs(y = "Sample Name") +
    theme(axis.text.x=element_text(angle=90, hjust=1)) +
    scale_y_discrete(drop=F)
  
  #add comulative plots totalling the rows and columns
  table(x)
  table(y)
  
  # Horizontal bar plot
  bpx<-ggplot(data=as.data.frame(sort(table(x))), aes(x=x, y=Freq, fill="Red")) + geom_bar(stat="identity") + 
    scale_y_reverse() + theme(legend.position = "none") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    theme(axis.text.x=element_text(angle=90, hjust=1))
  
  bpy<-ggplot(data=as.data.frame(table(y)), aes(x=y, y=Freq, fill="Red")) + geom_bar(stat="identity") + coord_flip() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                  panel.background = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position = "none")
  
  ggarrange(p, bpy, bpx, 
            ncol = 2, nrow = 2, align='hv', heights = c(2, 0.7), widths= c(2, 0.7))
  
  #print(p)
}

levels=sort(unique(c(levels(results.not_tumor[,2]),levels(results.tumor[,2]))))
#output = getPlotData(results.tumor, abundance.tumor, abundanceGene.tumor)
#output = getPlotData(results.tumor[which(!duplicated(results.tumor[,1])),], abundance.tumor, abundanceGene.tumor)
output = getPlotDataSimple(results.tumor)
x = output$x
x2 = output$x2
y = as.character(output$y)
size = output$size
color = output$color

#sort x by the number of hits
sortedNames = names(sort(table(x)))
newx=c()
newx2=c()
newy=c()
newsize=c()
newcolor=c()
for (i in 1:length(sortedNames)){
  newx=c(newx,x[which(x==sortedNames[i])])
  newx2=c(newx2,x2[which(x==sortedNames[i])])
  newy=c(newy,y[which(x==sortedNames[i])])
  newsize=c(newsize,size[which(x==sortedNames[i])])
  newcolor=c(newcolor,color[which(x==sortedNames[i])])
}
newy=factor(newy, levels=levels)

#replace ids with gene symbols for ncNATs if they are present
newx3=newx2
newx3[which(ensemblAnnotation[newx2,'hgnc_symbol']!="")]=ensemblAnnotation[newx2[which(ensemblAnnotation[newx2,'hgnc_symbol']!="")],'hgnc_symbol']

#replace the 2 ncNATs with a label "2 ncNATs*
newx4=newx3
newx4[33]="2 ncNATs"
newx4[34]="2 ncNATs"
#idReplace = matrix(ncol=1, colnames=unique(newx))
#idReplace[,1] = unique(newx2)

svg("candidates_tumor_with_ncNAT_2_ncNAT.svg", width=7, height=7)
plotCandidates(newx4, newy, newsize, newcolor, "Non-Tumor")
dev.off()

#output = getPlotData(results.not_tumor, abundance.not_tumor, abundanceGene.not_tumor)
output = getPlotDataSimple(results.not_tumor)
x = output$x
x2 = output$x2
y = as.character(output$y)
size = output$size
color = output$color

#sort x by the number of hits
sortedNames = names(sort(table(x)))
newxnt=c()
newx2nt=c()
newynt=c()
newsize=c()
newcolor=c()
for (i in 1:length(sortedNames)){
  newxnt=c(newxnt,x[which(x==sortedNames[i])])
  newx2nt=c(newx2nt,x2[which(x==sortedNames[i])])
  newynt=c(newynt,y[which(x==sortedNames[i])])
  newsize=c(newsize,size[which(x==sortedNames[i])])
  newcolor=c(newcolor,color[which(x==sortedNames[i])])
}
newynt=factor(newynt, levels=levels)

#replace ids with gene symbols for ncNATs if they are present
newx3nt=newx2nt
newx3nt[which(ensemblAnnotation[newx2nt,'hgnc_symbol']!="")]=ensemblAnnotation[newx2nt[which(ensemblAnnotation[newx2nt,'hgnc_symbol']!="")],'hgnc_symbol']

#replace the 2 ncNATs with a label "2 ncNATs*
newx4nt=newx3nt
newx4nt[85:92]="2 ncNATs"
newx4nt[85:92]="2 ncNATs"

svg("candidates_not_tumor_with_ncNAT_2_ncNAT.svg", width=10, height=7)
plotCandidates(newx4nt, newynt, newsize, newcolor, "Non-Tumor")
dev.off()


#plot the overlap ncNATs
overlapGenes=intersect(unique(newx), unique(newxnt))
overlapColor=c()
overlapx=c()
overlapy=c()
for (i in 1:length(overlapGenes)){
  tempx=newy[which(newx==overlapGenes[i])]
  tempy=newynt[which(newxnt==overlapGenes[i])]
  for (j in 1:length(tempx)){
    for(x in 1:length(tempy)){
      overlapColor=c(overlapColor,overlapGenes[i])
      overlapx=c(overlapx,as.character(tempx[j]))
      overlapy=c(overlapy,as.character(tempy[x]))
    }
  }
}
overlapx=factor(overlapx, levels=sort(union(unique(newy), unique(newynt))))
overlapy=factor(overlapy, levels=sort(union(unique(newy), unique(newynt))))

df = as.data.frame(cbind(as.character(overlapx),as.character(overlapy), overlapColor))
b <- ggplot(df, aes(x = overlapx, y = overlapy, color=overlapColor))
# Basic scatter plot
b <- b + geom_point(size=7, shape=18) +
  scale_x_discrete(drop=F) +
  scale_y_discrete(drop=F) +
  xlab("In Tumor") +
  ylab("In Non Tumor") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

svg("candidates_overlap.svg", width=7, height=5)
b
dev.off()

# 
# p <- qplot(overlapx, overlapy, geom="point", color=overlapColor, 
#            main=paste("ncNATs with candidates in both tumor and non tumour tissues")) + 
#   theme(legend.position = "none") + 
#   theme( axis.line = element_line(colour = "darkblue", size = 1, linetype = "solid")) +
#   scale_x_discrete(drop=F) +
#   scale_y_discrete(drop=F) +
#   xlab("In Tumor") +
#   ylab("In Non Tumor") +
#   theme(legend.position="bottom")

#p <- p + scale_x_discrete(breaks=c(1:nrow(results.tumor)),
#                     labels=results.tumor[,5])
#p <- p + scale_y_discrete(breaks=c(1:length(ylab)),labels=ylab)

# p <- p + geom_point(data = legend_bubbles,
#                     #  The "radius/50" was trial and error. Better way?
#                     aes(x = 5, y = 4 + radius, size = size),
#                     shape = 21, color = "black", fill = c('red','blue','red')) +
#   annotate("text", x = 4, y = 8, label = "Log 2 Fold Change", fontface = "bold") +
#   annotate("text", x = 4, y = 4, label = as.character(signif(min(size),3)), fontface = "bold") + 
#   annotate("text", x = 4, y = 5, label = paste('-',as.character(signif(((max(size)-min(size))/2)+min(size),3))), fontface = "bold") +
#   annotate("text", x = 4, y = 6, label = paste('>',as.character(signif(max(size)),3)), fontface = "bold") +
#   labs(x = "Regulated Genes") +
#   labs(y = "Sample Name") +
#   theme(axis.text.x=element_text(angle=90, hjust=1))