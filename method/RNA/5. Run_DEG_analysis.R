#Rscript Run_DEG_analysis.R

library(RUVSeq)
library(edgeR)
library(readr)
library(ggplot2)
library(DESeq2)
library(RColorBrewer)

fCount<-read.table("Output_correct.txt",header=T)
fCount<-round(fCount*25)

n<-colnames(fCount);lab<-c();assay<-c();exp<-c();tis<-c();sex<-c();age<-c()
for(j in 1:length(n)){
   lab[j]<-strsplit(n[j],"_")[[1]][1]
   assay[j]<-strsplit(n[j],"_")[[1]][2]
   exp[j]<-strsplit(n[j],"_")[[1]][3]
   tis[j]<-strsplit(n[j],"_")[[1]][5]
   sex[j]<-strsplit(n[j],"_")[[1]][6]
   age[j]<-strsplit(n[j],"_")[[1]][7]
}
group<-paste(lab,assay,exp,sep="_");group<-group[!(duplicated(group))]
group<-group[!(group %in% group[grep("Ctrl",group)])]
sex<-sex[!(duplicated(sex))];sex<-paste0("_",sex,"_")
Ctrl<-n[grep("_Ctrl_",n)]

d<-c()

for(lab_name in group){
    lab_table<-n[grep(lab_name,n)]
    for(i in sex){
        con1<-lab_table[grep(i,lab_table)]
        con2<-Ctrl[grep(i,Ctrl)]
        countdata<-fCount[,c(con1,con2)]
        l<-length(con1)
        keep=rowSums(cpm(countdata)>=1)>=l
        countdata=countdata[keep,]
        c1<-paste0(lab_name,i)
        c2<-paste0("Ctrl",i)
        x<-factor(c(rep(c1,length(con1)),rep(c2,length(con2))),levels=c(c2,c1))
        colData<-data.frame(row.names=colnames(countdata),x)
        ddscount <- DESeqDataSetFromMatrix(countData =countdata,colData = colData,design = ~x)
        dds <- DESeq(ddscount)
        results_df <- as.data.frame(results(dds))        
        results_df$cond1 <- c1
        results_df$cond2 <- c2
        results_df$gene <- rownames(results_df)
        results_df$lab<-strsplit(lab_name,"_")[[1]][1]
        results_df$type<-"liver_adult_pseudo_corrected"
        results_df$sig<-ifelse(abs(results_df$log2FoldChange) > log2(1.5) & results_df$padj < 0.001, T, F)
        results_df$sig[is.na(results_df$sig)] <- F
        results_df$direction <- ifelse(results_df$log2FoldChange >0, "up", "down")
        results_df$group <- paste("liver_adult_pseudo_corrected",results_df$lab,results_df$cond1,sep="_")
        results_df$tissue<-"liver"
        
        out<-paste0(c1,c2,"liver_adult_PseudoCorrect.txt")
        write.table(results_df,file=out,sep="\t",quote=F,row.names=F)
        results_df_f<-results_df[abs(results_df$log2FoldChange) > log2(1.5) & results_df$padj < 0.001,]
        out<-paste0(c1,c2,"liver_adult_PseudoCorrect_filter.txt")
        write.table(results_df_f,file=out,sep="\t",quote=F,row.names=F)

        d<-rbind(d,results_df) 
    }
}

d$cond2<-ifelse(d$cond2=="Ctrl_F_", "Ctrl_Female_aged", "Ctrl_Male_aged")
write.table(d,"Treatment_DEG.txt",sep="\t",quote=F,row.names=F)


