#Rscript RLE_normalization_and_RUVr_function.R read_count_table.txt Exposure_condition.txt

library(RUVSeq)
library(edgeR)
library(readr)
library(ggplot2)
library(DESeq2)
library(RColorBrewer)

args=commandArgs(T)
# path to the raw read count table
# format:
#                sample1  sample2 ... sampleN
#chr1,start,end   count    count       count


input<-paste0("/home/table/",args[1]) 
peak<-read.table(input,header=T)
fCount<-peak

# path to the file with condition group information
# format:
# Ctrl    liver_Female_adult        4       3
# Ctrl    liver_Male_adult        4       3
input<-paste0("/home/table/",args[2]) 
table<-read.table(input)

n<-ncol(fCount)

lab<-c("BA","BI","AL","DO","WK","MU","ZB")
for(lab_name in lab){
lab_table<-table[grep(lab_name,table$V1),]
lab_l<-length(lab_table$V1)

condition<-c()
for(i in 1:lab_l){
 con<-grep(lab_table[i,1],colnames(fCount))
 names<-colnames(fCount)[con]
 con2<-grep(lab_table[i,2],names)
 condition<-c(condition,names[con2])
}

countdata<-fCount[,condition]

cons<-as.character(paste(lab_table$V1,lab_table$V2,sep="_"))
x <- as.factor(c(rep(cons,lab_table$V3)))
set <- newSeqExpressionSet(as.matrix(countdata),phenoData=data.frame(x,row.names=colnames(countdata)))

design <- model.matrix(~x, data=pData(set))
y <- DGEList(counts=counts(set), group=x)
y=calcNormFactors(y, method="RLE")
y=estimateGLMCommonDisp(y, design)
y=estimateGLMTagwiseDisp(y, design)
set=newSeqExpressionSet(as.matrix(y$counts),phenoData=data.frame(x,row.names=colnames(countdata)))

#raw
colors <- brewer.pal(8, "Set2")
name1<-paste0(lab_name,"_RUVr_rawdata_RLE")
suppressMessages(pdf(paste0(lab_name,"_RUVr_rawdata_RLE.pdf")))
suppressMessages(plotRLE(set, outline=FALSE, ylim=c(-2, 2), col=colors[x],main=name1))
invisible(dev.off())
name2<-paste0(lab_name,"_RUVr_rawdata_PCA")
suppressMessages(pdf(paste0(lab_name,"_RUVr_rawdata_PCA.pdf")))
suppressMessages(plotPCA(set, col=colors[x], main=name2,cex=0.4))
invisible(dev.off())

out<-paste0(lab_name,"_raw.bed")
write.table(y$counts,file=out,sep="\t",quote=F)

fit=glmFit(y, design)
res=residuals(fit, type="deviance")
seqUQ=betweenLaneNormalization(set, which="upper")

#k<-3
run_RUVr<-function(k){
set2 <- RUVr(seqUQ, rownames(set), k=k, res)
ndddd=normCounts(set2)

name1<-paste0(lab_name,"_RUVr k=",k)
suppressMessages(pdf(paste0(lab_name,"_RUVr_k", k, "_RLE.pdf")))
suppressMessages(plotRLE(set2, outline=FALSE, ylim=c(-2, 2), col=colors[x], main=name1))
invisible(dev.off())
name2<-paste0(lab_name,"_RUVr k=",k)
suppressMessages(pdf(paste0(lab_name,"_RUVr_k", k, "_PCA.pdf")))
suppressMessages(plotPCA(set2, col=colors[x], main=name2,cex=0.4))
invisible(dev.off())

out<-paste0(lab_name,"_RUVr_k",k,".bed")
write.table(ndddd,file=out,sep="\t",quote=F)
}
for (runr in seq(1,4,1)){
  run_RUVr(runr)
}

}
