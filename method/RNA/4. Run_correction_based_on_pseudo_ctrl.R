#Rscript Run_correction_based_on_pseudo_ctrl.R

library(edgeR)

args=commandArgs(T)

data1<-read.table("LAB_RUVr_k3_pseudo.txt",header=T)
data2<-read.table("LAB_RUVr_k3.bed",header=T)
data<-cpm(data2+1)
ave<-cpm(data1+1)

n<-colnames(data)
names<-n
lab<-c();assay<-c();exp<-c();tis<-c();sex<-c();age<-c()
for(j in 1:length(n)){
   lab[j]<-strsplit(n[j],"_")[[1]][1]
   assay[j]<-strsplit(n[j],"_")[[1]][2]
   exp[j]<-strsplit(n[j],"_")[[1]][3]
   tis[j]<-strsplit(n[j],"_")[[1]][5]
   sex[j]<-strsplit(n[j],"_")[[1]][6]
   age[j]<-strsplit(n[j],"_")[[1]][7]
}

group<-paste(lab,assay,sep="_");group<-group[!(duplicated(group))]
table<-as.data.frame(rep(group,rep(2,length(group))));colnames(table)<-"V1"
table$V2<-"Ctrl"
table$V3<-c("_F_","_M_")

n<-length(table$V1)
d<-c()
x<-c()

i<-1
a<-names[grep(as.character(table[i,1]),names)]
b<-a[grep(as.character(table[i,3]),a)]
LAB<-data[,b]
r<-b[grep("Ctrl",b)]
LAB_r<-LAB[,r]

a<-ave[,paste0("X", as.character(table[i, 3]))]
b<-rowSums(LAB_r)/ncol(LAB_r)
c<-as.data.frame(b/a)
colnames(c)<-"cor"
LAB_c<-(LAB)/c$cor
d<-LAB_c

for(i in 2:n){

  a<-names[grep(as.character(table[i,1]),names)]
  b<-a[grep(as.character(table[i,3]),a)]
  LAB<-data[,b]
  r<-b[grep("Ctrl",b)]
  LAB_r<-LAB[,r]
  
  a<-ave[,paste0("X", as.character(table[i, 3]))]
  b<-rowSums(LAB_r)/ncol(LAB_r)
  c<-as.data.frame(b/a)
  colnames(c)<-"cor"
  LAB_c<-(LAB)/c$cor
  d<-cbind(d,LAB_c)

}

d<-round(d,4)
write.table(d,"Output_correct.txt",sep="\t",quote=F)

