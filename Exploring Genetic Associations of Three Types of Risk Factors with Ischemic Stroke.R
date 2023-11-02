#MR analysis
library(TwoSampleMR)
Exp_dat <- extract_instruments(outcomes = 'ieu-b-110')
Out_dat <- extract_outcome_data(snps =Exp_dat $SNP, outcomes = 'ebi-a-GCST005843')
Out_dat<-subset(Out_dat,pval.outcome>0.04)
dat <- harmonise_data(
  exposure_dat = Exp_dat,
  outcome_dat = Out_dat,
  action = 2
)
res <- mr(dat)
res
OR <- generate_odds_ratios(res)
OR
library(MRPRESSO)
mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure",
          OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = dat, NbDistribution = 1200,
          SignifThreshold = 0.05)
het <- mr_heterogeneity(dat)
het
pleio <- mr_pleiotropy_test(dat)
pleio

#Colocalization analysis
library(TwoSampleMR)
library(gwasglue)
library(gwasvcf)
library(dplyr)
library(VariantAnnotation)
library(gassocplot2)
library(ieugwasr)
library(coloc)
top <- ieugwasr::tophits('ieu-b-110') %>% arrange(p)
top

chrpos <- paste0(top$chr[1], ":", top$position[1] - 90000, "-", top$position[1] + 90000)
chrpos
out <- ieugwasr_to_coloc(id1 = 'ieu-b-110', id2 = 'ebi-a-GCST005843', chrompos = chrpos)

res <- coloc::coloc.abf(out[[1]], out[[2]])

temp <- coloc_to_gassocplot(out)
markers<-data.frame(temp$markers)
z<-data.frame(temp$z)
corr<-data.frame(temp$corr)
gassocplot2::stack_assoc_plot(markers, z, corr, traits=temp$traits)

#GEO analysis
#data extraction
library(limma)              
expFile="geneMatrix.txt"   
conFile="s1.txt"             
treatFile="s2.txt"          
setwd("")     
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
s1=read.table(conFile, header=F, sep="\t", check.names=F)
sampleName1=as.vector(s1[,1])
conData=data[,sampleName1]
s2=read.table(treatFile, header=F, sep="\t", check.names=F)
sampleName2=as.vector(s2[,1])
treatData=data[,sampleName2]
rt=cbind(conData, treatData)
qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
if(LogC){
  rt[rt<0]=0
  rt=log2(rt+1)}
data=normalizeBetweenArrays(rt)
conNum=ncol(conData)
treatNum=ncol(treatData)
Type=c(rep("Control",conNum),rep("Treat",treatNum))
outData=rbind(id=paste0(colnames(data),"_",Type),data)
write.table(outData, file="normalize.txt", sep="\t", quote=F, col.names=F)

#Expression data of associated genes
library(limma)              
expFile="normalize.txt"    
geneFile="gene.txt"        
setwd("")     
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
gene=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(as.vector(gene[,1]), rownames(data))
geneExp=data[sameGene,]
out=rbind(ID=colnames(geneExp), geneExp)
write.table(out, file="crgGeneExp.txt", sep="\t", quote=F, col.names=F)

#Differential expression analysis of associated genes
library(limma)
library(pheatmap)
library(reshape2)
library(ggpubr)
expFile="crgGeneExp.txt"      
setwd("")   
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
exp=data

Type=gsub("(.*)\\_(.*)", "\\2", colnames(data))

sigVec=c()
sigGeneVec=c()
for(i in row.names(data)){
  test=wilcox.test(data[i,] ~ Type)
  pvalue=test$p.value
  Sig=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
  if(pvalue<0.05){
    sigVec=c(sigVec, paste0(i, Sig))
    sigGeneVec=c(sigGeneVec, i)}
}
data=data[sigGeneVec,]
outTab=rbind(ID=colnames(data), data)
write.table(outTab, file="diffGeneExp.txt", sep="\t", quote=F, col.names=F)
row.names(data)=sigVec

names(Type)=colnames(data)
Type=as.data.frame(Type)
pdf(file="heatmap.pdf", width=7, height=4.5)
pheatmap(data,
         annotation=Type,
         color = colorRampPalette(c(rep("blue",2), "white", rep("red",2)))(100),
         cluster_cols =F,
         cluster_rows =T,
         scale="row",
         show_colnames=F,
         show_rownames=T,
         fontsize=7,
         fontsize_row=7,
         fontsize_col=7)
dev.off()

exp=as.data.frame(t(exp))
exp=cbind(exp, Type=Type)
data=melt(exp, id.vars=c("Type"))
colnames(data)=c("Type", "Gene", "Expression")

p=ggboxplot(data, x="Gene", y="Expression", color = "Type", 
            xlab="",
            ylab="Gene expression",
            legend.title="Type",
            palette = c("blue", "red"),
            add="point",
            width=0.8)
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=Type),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif")

pdf(file="boxplot.pdf", width=7, height=5)
print(p1)
dev.off()

#chromosomal location analysis of DEAGs
library("RCircos")      
setwd("")      

cytoBandIdeogram=read.table("refer.txt", header=T, sep="\t", check.names=F)
chr.exclude <- NULL
cyto.info <- cytoBandIdeogram
tracks.inside <- 4
tracks.outside <- 0
RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside)

rcircos.params <- RCircos.Get.Plot.Parameters()
rcircos.params$text.size=0.7
rcircos.params$point.size=5
RCircos.Reset.Plot.Parameters(rcircos.params)

pdf(file="RCircos.pdf", width=7, height=7)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()


RCircos.Gene.Label.Data=read.table("Rcircos.geneLabel.txt", header=T, sep="\t", check.names=F)
name.col <- 4
side <- "in"
track.num <- 1
RCircos.Gene.Connector.Plot(RCircos.Gene.Label.Data, track.num, side)
track.num <- 2
RCircos.Gene.Name.Plot(RCircos.Gene.Label.Data, name.col, track.num, side)
dev.off()

#Co-expression analysis of DEAGs
library(corrplot)
library(circlize)

inputFile="diffGeneExp.txt"   
setwd("")   

data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)

group=gsub("(.*)\\_(.*)", "\\2", colnames(data))
data=data[,group=="Treat",drop=F]
rt=t(data)

cor1=cor(rt)

col = c(rgb(1,0,0,seq(1,0,length=32)),rgb(0,1,0,seq(0,1,length=32)))
cor1[cor1==1]=0
c1 = ifelse(c(cor1)>=0,rgb(1,0,0,abs(cor1)),rgb(0,1,0,abs(cor1)))
col1 = matrix(c1,nc=ncol(rt))

pdf(file="circos.pdf", width=7, height=7)
par(mar=c(2,2,2,4))
circos.par(gap.degree=c(3,rep(2, nrow(cor1)-1)), start.degree = 180)
chordDiagram(cor1, grid.col=rainbow(ncol(rt)), col=col1, transparency = 0.5, symmetric = T)
par(xpd=T)
colorlegend(col, vertical = T,labels=c(1,0,-1),xlim=c(1.1,1.3),ylim=c(-0.4,0.4))
dev.off()
circos.clear()

pdf(file="corrplot.pdf", width=7, height=7)
corrplot(cor1,
         method = "pie",
         order = "hclust",
         type = "upper",
         col=colorRampPalette(c("green", "white", "red"))(50)
)
dev.off()

#Analysis of immune cells in IS samples
inputFile="normalize.txt"      
setwd("")      
source("geoCRG11.CIBERSORT.R")   

outTab=CIBERSORT("ref.txt", inputFile, perm=1000, QN=T)

outTab=outTab[outTab[,"P-value"]<0.1,]
outTab=as.matrix(outTab[,1:(ncol(outTab)-3)])
outTab=rbind(id=colnames(outTab),outTab)
write.table(outTab, file="CIBERSORT-Results.txt", sep="\t", quote=F, col.names=F)

library(reshape2)
library(ggpubr)

inputFile="CIBERSORT-Results.txt"  
setwd("")     

rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)

con=grepl("_Control", rownames(rt), ignore.case=T)
data=t(rbind(conData,treatData))


conData=rt[con,]
treatData=rt[treat,]
conNum=nrow(conData)
treatNum=nrow(treatData)

pdf(file="barplot.pdf", width=14.5, height=8.5)
col=rainbow(nrow(data), s=0.7, v=0.7)
par(las=1,mar=c(8,5,4,16),mgp=c(3,0.1,0),cex.axis=1.5)
a1=barplot(data,col=col,xaxt="n",yaxt="n",ylab="Relative Percent",cex.lab=1.8)
a2=axis(2,tick=F,labels=F)
axis(2,a2,paste0(a2*100,"%"))
par(srt=0,xpd=T)
rect(xleft = a1[1]-0.5, ybottom = -0.01, xright = a1[conNum]+0.5, ytop = -0.06,col="green")
text(a1[conNum]/2,-0.035,"Control",cex=1.8)
rect(xleft = a1[conNum]+0.5, ybottom = -0.01, xright =a1[length(a1)]+0.5, ytop = -0.06,col="red")
text((a1[length(a1)]+a1[conNum])/2,-0.035,"Treat",cex=1.8)
ytick2 = cumsum(data[,ncol(data)])
ytick1 = c(0,ytick2[-length(ytick2)])
legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(data),col=col,pch=15,bty="n",cex=1.2)
dev.off()

Type=gsub("(.*)\\_(.*)", "\\2", rownames(rt))
data=cbind(as.data.frame(t(data)), Type)
data=melt(data, id.vars=c("Type"))
colnames(data)=c("Type", "Immune", "Expression")
group=levels(factor(data$Type))
bioCol=c("#0066FF","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(group)]
boxplot=ggboxplot(data, x="Immune", y="Expression", color="Type",
                  xlab="",
                  ylab="Fraction",
                  legend.title="Type",
                  add="point",
                  width=0.8,
                  palette=bioCol)+
  rotate_x_text(50)+
  stat_compare_means(aes(group=Type),symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "")), label="p.signif")
pdf(file="immune.diff.pdf", width=8, height=6)
print(boxplot)
dev.off()

library(limma)
library(reshape2)
library(tidyverse)
library(ggplot2)

expFile="diffGeneExp.txt"          
immFile="CIBERSORT-Results.txt"     
setwd(")     

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

group=gsub("(.*)\\_(.*)", "\\2", colnames(data))
data=data[,group=="Treat",drop=F]
data=t(data)
immune=read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(data), row.names(immune))
data=data[sameSample,,drop=F]
immune=immune[sameSample,,drop=F]

outTab=data.frame()
for(cell in colnames(immune)){
  if(sd(immune[,cell])==0){next}
  for(gene in colnames(data)){
    x=as.numeric(immune[,cell])
    y=as.numeric(data[,gene])
    corT=cor.test(x,y,method="spearman")
    cor=corT$estimate
    pvalue=corT$p.value
    text=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
    outTab=rbind(outTab,cbind(Gene=gene, Immune=cell, cor, text, pvalue))
  }
}

outTab$cor=as.numeric(outTab$cor)
pdf(file="cor.pdf", width=7, height=5)
ggplot(outTab, aes(Immune, Gene)) + 
  geom_tile(aes(fill = cor), colour = "grey", size = 1)+
  scale_fill_gradient2(low = "#5C5DAF", mid = "white", high = "#EA2E2D") + 
  geom_text(aes(label=text),col ="black",size = 3) +
  theme_minimal() +    #ȥ?����?
  theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "bold"),   #x??????
        axis.text.y = element_text(size = 8, face = "bold")) +       #y??????
  labs(fill =paste0("***  p<0.001","\n", "**  p<0.01","\n", " *  p<0.05","\n", "\n","Correlation")) +   #????ͼ??
  scale_x_discrete(position = "bottom")      #X????????ʾλ??
dev.off()

#Clutering analysis of DEAGs
library(ConsensusClusterPlus)     
expFile="diffGeneExp.txt"          
workDir=""     
setwd(workDir)     

data=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
data=as.matrix(data)

group=sapply(strsplit(colnames(data),"\\_"), "[", 2)
data=data[,group=="Treat"]

maxK=9    
results=ConsensusClusterPlus(data,
                             maxK=maxK,
                             reps=50,
                             pItem=0.8,
                             pFeature=1,
                             title=workDir,
                             clusterAlg="km",
                             distance="euclidean",
                             seed=123456,
                             plot="png")

calcICL(results, title="consensusScore", plot="png")

clusterNum=2     
cluster=results[[clusterNum]][["consensusClass"]]
cluster=as.data.frame(cluster)
colnames(cluster)=c("Cluster")
cluster$Cluster=paste0("C", cluster$Cluster)
outTab=cbind(t(data), cluster)
outTab=rbind(ID=colnames(outTab), outTab)
write.table(outTab, file="cluster.txt", sep="\t", quote=F, col.names=F)

#Differential expression analysis between DEAG clusters
library(limma)
library(pheatmap)
library(reshape2)
library(ggpubr)

clusterFile="cluster.txt"     
setwd("")    

rt=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
rt=rt[order(rt$Cluster),]

data=t(rt[,1:(ncol(rt)-1),drop=F])
Type=rt[,ncol(rt),drop=F]

bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
ann_colors=list()
crgCluCol=bioCol[1:length(levels(factor(Type$Cluster)))]
names(crgCluCol)=levels(factor(Type$Cluster))
ann_colors[["Cluster"]]=crgCluCol

pdf("heatmap.pdf", width=7, height=4.5)
pheatmap(data,
         annotation=Type,
         annotation_colors = ann_colors,
         color = colorRampPalette(c(rep("blue",2), "white", rep("red",2)))(100),
         cluster_cols =F,
         cluster_rows =T,
         scale="row",
         show_colnames=F,
         show_rownames=T,
         fontsize=7,
         fontsize_row=7,
         fontsize_col=7)
dev.off()

data=melt(rt, id.vars=c("Cluster"))
colnames(data)=c("Cluster", "Gene", "Expression")

p=ggboxplot(data, x="Gene", y="Expression", color = "Cluster",
            xlab="",
            ylab="Gene expression",
            legend.title="Cluster",
            palette = crgCluCol,
            width=0.8,
            add="point")
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=Cluster),
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif")

pdf(file="boxplot.pdf", width=7, height=5)
print(p1)
dev.off()

#PCA
library(limma)
library(ggplot2)

clusterFile="cluster.txt"      
setwd("")     

rt=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
data=rt[,1:(ncol(rt)-1),drop=F]
Cluster=as.vector(rt[,ncol(rt)])

#PCA
data.pca=prcomp(data)
pcaPredict=predict(data.pca)
PCA=data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], Cluster=Cluster)
PCA.mean=aggregate(PCA[,1:2], list(Cluster=PCA$Cluster), mean)


bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
crgCluCol=bioCol[1:length(levels(factor(Cluster)))]



veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) {
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}
df_ell <- data.frame()
for(g in levels(factor(PCA$Cluster))){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(PCA[PCA$Cluster==g,],
                                                   veganCovEllipse(cov.wt(cbind(PC1,PC2),
                                                                          wt=rep(1/length(PC1),length(PC1)))$cov,
                                                                   center=c(mean(PC1),mean(PC2))))), Cluster=g))
}


pdf(file="PCA.pdf", width=6.5, height=5)
ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color = Cluster)) +
  scale_colour_manual(name="Cluster", values =crgCluCol)+
  theme_bw()+
  theme(plot.margin=unit(rep(1.5,4),'lines'))+
  geom_path(data=df_ell, aes(x=PC1, y=PC2, colour=Cluster), size=1, linetype=2)+
  annotate("text",x=PCA.mean$PC1, y=PCA.mean$PC2, label=PCA.mean$Cluster, cex=7)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

#GSVA
library(reshape2)
library(ggpubr)
library(limma)
library(GSEABase)
library(GSVA)

expFile="normalize.txt"           
clusterFile="cluster.txt"           
gmtFile="c2.cp.kegg.symbols.gmt"   
setwd("")     

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

group=gsub("(.*)\\_(.*)", "\\2", colnames(data))
data=data[,group=="Treat",drop=F]
geneSets=getGmt(gmtFile, geneIdType=SymbolIdentifier())

ssgseaScore=gsva(data, geneSets, method='gsva')
#??GSVA?Ĵ??ֽ??н???
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
ssgseaScore=normalize(ssgseaScore)

cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
nameC1=row.names(cluster[cluster$Cluster=="C1",,drop=F])
nameC2=row.names(cluster[cluster$Cluster=="C2",,drop=F])
dataC1=ssgseaScore[,nameC1,drop=F]
dataC2=ssgseaScore[,nameC2,drop=F]
conNum=ncol(dataC1)
treatNum=ncol(dataC2)
data=cbind(dataC1, dataC2)
Type=c(rep("C1",conNum), rep("C2",treatNum))


outTab=data.frame()
for(i in row.names(data)){
  test=t.test(data[i,] ~ Type)
  pvalue=test$p.value
  t=test$statistic
  if(pvalue<0.05){
    Sig=ifelse(pvalue>0.05, "Not", ifelse(t>0,"Up","Down"))
    outTab=rbind(outTab, cbind(Pathway=i, t, pvalue, Sig))
  }
}

termNum=10      
outTab=outTab[order(outTab$t),]
outTab=outTab[c(1:termNum,(nrow(outTab)-termNum):nrow(outTab)),]
pdf(file="barplot.pdf", width=9, height=6)
outTab$t=as.numeric(outTab$t)
outTab$Sig=factor(outTab$Sig, levels=c("Down", "Up"))
gg1=ggbarplot(outTab, x="Pathway", y="t", fill = "Sig", color = "white",
              palette=c("blue3", "red3"), sort.val = "asc", sort.by.groups = T,
              rotate=TRUE, legend="right", title="",
              xlab="Term", ylab="t value of GSVA score, C2 vs C1",  legend.title="Group", x.text.angle=60)
print(gg1)
dev.off()

#Machine learning models
library(caret)
library(DALEX)
library(ggplot2)
library(randomForest)
library(kernlab)
library(xgboost)
library(pROC)

set.seed(123)   
inputFile="normalize.txt"      
geneFile="interGenes.txt"      
setwd("")    
data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)

geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
data=data[as.vector(geneRT[,1]),]
row.names(data)=gsub("-", "_", row.names(data))

data=t(data)
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))
data=as.data.frame(data)
data$Type=group

inTrain<-createDataPartition(y=data$Type, p=0.7, list=F)
train<-data[inTrain,]
test<-data[-inTrain,]

control=trainControl(method="repeatedcv", number=5, savePredictions=TRUE)
mod_rf = train(Type ~ ., data = train, method='rf', trControl = control)

vm=train(Type ~., data = train, method = "svmRadial", prob.model=TRUE, trControl=control)

mod_xgb=train(Type ~., data = train, method = "xgbDART", trControl=control)

mod_glm=train(Type ~., data = train, method = "glm", family="binomial", trControl=control)

p_fun=function(object, newdata){
  predict(object, newdata=newdata, type="prob")[,2]
}
yTest=ifelse(test$Type=="Control", 0, 1)

explainer_rf=explain(mod_rf, label = "RF",
                     data = test, y = yTest,
                     predict_function = p_fun,
                     verbose = FALSE)
mp_rf=model_performance(explainer_rf)

explainer_svm=explain(mod_svm, label = "SVM",
                      data = test, y = yTest,
                      predict_function = p_fun,
                      verbose = FALSE)
mp_svm=model_performance(explainer_svm)

explainer_xgb=explain(mod_xgb, label = "XGB",
                      data = test, y = yTest,
                      predict_function = p_fun,
                      verbose = FALSE)
mp_xgb=model_performance(explainer_xgb)

explainer_glm=explain(mod_glm, label = "GLM",
                      data = test, y = yTest,
                      predict_function = p_fun,
                      verbose = FALSE)
mp_glm=model_performance(explainer_glm)

pdf(file="residual.pdf", width=6, height=6)
p1 <- plot(mp_rf, mp_svm, mp_xgb, mp_glm)
print(p1)
dev.off()

pdf(file="boxplot.pdf", width=6, height=6)
p2 <- plot(mp_rf, mp_svm, mp_xgb, mp_glm, geom = "boxplot")
print(p2)
dev.off()

pred1=predict(mod_rf, newdata=test, type="prob")
pred2=predict(mod_svm, newdata=test, type="prob")
pred3=predict(mod_xgb, newdata=test, type="prob")
pred4=predict(mod_glm, newdata=test, type="prob")
roc1=roc(yTest, as.numeric(pred1[,2]))
roc2=roc(yTest, as.numeric(pred2[,2]))
roc3=roc(yTest, as.numeric(pred3[,2]))
roc4=roc(yTest, as.numeric(pred4[,2]))
pdf(file="ROC.pdf", width=5, height=5)
plot(roc1, print.auc=F, legacy.axes=T, main="", col="red")
plot(roc2, print.auc=F, legacy.axes=T, main="", col="blue", add=T)
plot(roc3, print.auc=F, legacy.axes=T, main="", col="green", add=T)
plot(roc4, print.auc=F, legacy.axes=T, main="", col="yellow", add=T)
legend('bottomright',
       c(paste0('RF: ',sprintf("%.03f",roc1$auc)),
         paste0('SVM: ',sprintf("%.03f",roc2$auc)),
         paste0('XGB: ',sprintf("%.03f",roc3$auc)),
         paste0('GLM: ',sprintf("%.03f",roc4$auc))),
       col=c("red","blue","green","yellow"), lwd=2, bty = 'n')
dev.off()

importance_rf<-variable_importance(
  explainer_rf,
  loss_function = loss_root_mean_square
)
importance_svm<-variable_importance(
  explainer_svm,
  loss_function = loss_root_mean_square
)
importance_glm<-variable_importance(
  explainer_glm,
  loss_function = loss_root_mean_square
)
importance_xgb<-variable_importance(
  explainer_xgb,
  loss_function = loss_root_mean_square
)

pdf(file="importance.pdf", width=7, height=10)
plot(importance_rf[c(1,(ncol(data)-8):(ncol(data)+1)),],
     importance_svm[c(1,(ncol(data)-8):(ncol(data)+1)),],
     importance_xgb[c(1,(ncol(data)-8):(ncol(data)+1)),],
     importance_glm[c(1,(ncol(data)-8):(ncol(data)+1)),])
dev.off()

geneNum=5    
write.table(importance_rf[(ncol(data)-geneNum+2):(ncol(data)+1),], file="importanceGene.RF.txt", sep="\t", quote=F, row.names=F)
write.table(importance_svm[(ncol(data)-geneNum+2):(ncol(data)+1),], file="importanceGene.SVM.txt", sep="\t", quote=F, row.names=F)
write.table(importance_xgb[(ncol(data)-geneNum+2):(ncol(data)+1),], file="importanceGene.XGB.txt", sep="\t", quote=F, row.names=F)
write.table(importance_glm[(ncol(data)-geneNum+2):(ncol(data)+1),], file="importanceGene.GLM.txt", sep="\t", quote=F, row.names=F)

#Nomogram
library(rms)
library(rmda)

inputFile="normalize.txt"             
geneFile="importanceGene..txt"    
setwd("")    
data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
row.names(data)=gsub("-", "_", row.names(data))

geneRT=read.table(geneFile, header=T, sep="\t", check.names=F)
data=data[as.vector(geneRT[,1]),]
data=t(data)
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))
rt=cbind(as.data.frame(data), Type=group)
paste(colnames(data), collapse="+")

ddist=datadist(rt)
options(datadist="ddist")

lrmModel=lrm(Type~ ENC1+AUH+UBE2E4P+R3HDM1+ITPKB, data=rt, x=T, y=T)
nomo=nomogram(lrmModel, fun=plogis,
              fun.at=c(0.0001,0.1,0.3,0.5,0.7,0.9,0.99),
              lp=F, funlabel="Risk of Disease")
pdf("Nomo.pdf", width=8, height=6)
plot(nomo)
dev.off()

cali=calibrate(lrmModel, method="boot", B=1000)
pdf("Calibration.pdf", width=5.5, height=5.5)
plot(cali,
     xlab="Predicted probability",
     ylab="Actual probability", sub=F)
dev.off()

rt$Type=ifelse(rt$Type=="Control", 0, 1)
dc=decision_curve(Type ~ ENC1+AUH+UBE2E4P+R3HDM1+ITPKB, data=rt, 
                  family = binomial(link ='logit'),
                  thresholds= seq(0,1,by = 0.01),
                  confidence.intervals = 0.95)
pdf(file="DCA.pdf", width=5.5, height=5.5)
plot_decision_curve(dc,
                    curve.names="Model",
                    xlab="Threshold probability",
                    cost.benefit.axis=T,
                    col="red",
                    confidence.intervals=FALSE,
                    standardize=FALSE)
dev.off()

#Differential expression analysis between DEAG clusters
library(limma) 
library(VennDiagram)

logFCfilter=1                 
adj.P.Val.Filter=0.05         
expFile="normalize.txt"      
cluFile="Cluster.txt"     
setwd("")     

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

cluster=read.table(cluFile, header=T, sep="\t", check.names=F, row.names=1)

sameSample=intersect(colnames(data), row.names(cluster))
data=data[,sameSample]
cluster=cluster[sameSample, "Cluster"]

geneList=list()
Type=as.vector(cluster)
design=model.matrix(~0+factor(Type))
colnames(design)=levels(factor(Type))
comp=combn(levels(factor(Type)), 2)
allDiffGenes=c()

for(i in 1:ncol(comp)){
  fit=lmFit(data, design)
  contrast=paste0(comp[2,i], "-", comp[1,i])
  #print(contrast)
  cont.matrix=makeContrasts(contrast, levels=design)
  fit2=contrasts.fit(fit, cont.matrix)
  fit2=eBayes(fit2)
  
  allDiff=topTable(fit2,adjust='fdr',number=200000)
  allDiffOut=rbind(id=colnames(allDiff),allDiff)
  write.table(allDiffOut, file=paste0(contrast, ".all.txt"), sep="\t", quote=F, col.names=F)
  
  diffSig=allDiff[with(allDiff, (abs(logFC)>logFCfilter & adj.P.Val < adj.P.Val.Filter )), ]
  diffSigOut=rbind(id=colnames(diffSig),diffSig)
  write.table(diffSigOut, file=paste0(contrast, ".diff.txt"), sep="\t", quote=F, col.names=F)
  geneList[[contrast]]=row.names(diffSig)
}

venn.plot=venn.diagram(geneList,filename=NULL,fill=rainbow(length(geneList)) )
pdf(file="venn.pdf", width=5, height=5)
grid.draw(venn.plot)
dev.off()

interGenes=Reduce(intersect,geneList)
write.table(file="interGene.txt",interGenes,sep="\t",quote=F,col.names=F,row.names=F)

interGeneExp=data[interGenes,]
interGeneExp=rbind(id=colnames(interGeneExp), interGeneExp)
write.table(interGeneExp, file="interGeneExp.txt", sep="\t", quote=F, col.names=F)

#DEG cluster
library(limma)
library(ConsensusClusterPlus)

expFile="interGeneExp.txt"      
workDir=""     
setwd(workDir)    

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

maxK=9      
results=ConsensusClusterPlus(data,
                             maxK=maxK,
                             reps=50,
                             pItem=0.8,
                             pFeature=1,
                             title=workDir,
                             clusterAlg="km",
                             distance="euclidean",
                             seed=123456,
                             plot="png")

clusterNum=2        
cluster=results[[clusterNum]][["consensusClass"]]
cluster=as.data.frame(cluster)
colnames(cluster)=c("geneCluster")
letter=c("A","B","C","D","E","F","G")
uniqClu=levels(factor(cluster$geneCluster))
cluster$geneCluster=letter[match(cluster$geneCluster, uniqClu)]
clusterOut=rbind(ID=colnames(cluster), cluster)
write.table(clusterOut, file="geneCluster.txt", sep="\t", quote=F, col.names=F)

#Validation of ML model
library(caret)
library(DALEX)
library(ggplot2)
library(randomForest)
library(kernlab)
library(xgboost)
library(pROC)

set.seed(123)      
inputFile="test.normalize.txt"         
geneFile="importanceGene..txt"     
setwd("")  

data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
row.names(data)=gsub("-", "_", row.names(data))

geneRT=read.table(geneFile, header=T, sep="\t", check.names=F)
data=data[as.vector(geneRT[,1]),]

data=t(data)
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))
data=as.data.frame(data)
data$Type=group

inTrain<-createDataPartition(y=data$Type, p=0.7, list=F)
train<-data[inTrain,]
test<-data[-inTrain,]

control=trainControl(method="repeatedcv", number=5, savePredictions=TRUE)
if(geneFile=="importanceGene.RF.txt"){

  model=train(Type ~ ., data = train, method='rf', trControl = control)
}else if(geneFile=="importanceGene.SVM.txt"){

  model=train(Type ~., data = train, method = "svmRadial", prob.model=TRUE, trControl=control)
}else if(geneFile=="importanceGene.XGB.txt"){

  model=train(Type ~., data = train, method = "xgbDART", trControl=control)
}else if(geneFile=="importanceGene.GLM.txt"){
 
  model=train(Type ~., data = train, method = "glm", family="binomial", trControl=control)
}

yTest=ifelse(test$Type=="Control", 0, 1)
pred1=predict(model, newdata=test, type="prob")
roc1=roc(yTest, as.numeric(pred1[,2]))
ci1=ci.auc(roc1, method="bootstrap")
ciVec=as.numeric(ci1)
pdf(file="ROC.pdf", width=5, height=5)
plot(roc1, print.auc=T, legacy.axes=T, main="", col="red")
text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
dev.off()

