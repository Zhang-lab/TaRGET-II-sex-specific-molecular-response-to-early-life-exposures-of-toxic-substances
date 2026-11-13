#Rscript RLE_normalization_and_RUVr_function.R read_count_table.txt

library(RUVSeq)
library(edgeR)
library(readr)
library(ggplot2)
library(DESeq2)
library(RColorBrewer)

args=commandArgs(T)
input<-args[1] # path to the raw read count table
peak<-read.table(input,header=T)
fCount<-peak

n<-colnames(fCount);lab<-c();assay<-c();exp<-c();tis<-c();sex<-c();age<-c()
for(j in 1:length(n)){
   lab[j]<-strsplit(n[j],"_")[[1]][1]
   assay[j]<-strsplit(n[j],"_")[[1]][2]
   exp[j]<-strsplit(n[j],"_")[[1]][3]
   tis[j]<-strsplit(n[j],"_")[[1]][5]
   sex[j]<-strsplit(n[j],"_")[[1]][6]
   age[j]<-strsplit(n[j],"_")[[1]][7]
}
group<-paste(lab,assay,exp,sep="_");group<-group[!(duplicated(group))];lab<-lab[!duplicated(lab)]

for(lab_name in lab){
    lab_table<-group[grep(lab_name,group)]
    condition<-c();length_num<-c();cons<-c()
    for(i in 1:length(lab_table)){
        con<-grep(lab_table[i],colnames(fCount))
        names<-colnames(fCount)[con]
        con1<-grep("_F_",names)
        con2<-grep("_M_",names)
        condition<-c(condition,names[con1],names[con2])
        length_num<-c(length_num,length(con1),length(con2))
        cons<-c(cons,paste0(lab_table[i],"_F"),paste0(lab_table[i],"_M"))
    }
    countdata<-fCount[,condition]
    x <- as.factor(c(rep(cons,length_num)))
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

    run_RUVr<-function(k){
       set2 <- RUVr(seqUQ, rownames(set), k=k, res)
       ndddd=normCounts(set2)
       out<-paste0(lab_name,"_RUVr_k",k,".bed")
       write.table(ndddd,file=out,sep="\t",quote=F)

       name1<-paste0(lab_name,"_RUVr k=",k)
       suppressMessages(pdf(paste0(lab_name,"_RUVr_k", k, "_RLE.pdf")))
       suppressMessages(plotRLE(set2, outline=FALSE, ylim=c(-2, 2), col=colors[x], main=name1))
       invisible(dev.off())
       name2<-paste0(lab_name,"_RUVr k=",k)
       suppressMessages(pdf(paste0(lab_name,"_RUVr_k", k, "_PCA.pdf")))
       suppressMessages(plotPCA(set2, col=colors[x], main=name2,cex=0.4))
       invisible(dev.off())
    }

    for (runr in seq(1,4,1)){
         run_RUVr(runr)
    }

}


