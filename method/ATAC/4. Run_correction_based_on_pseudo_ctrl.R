#Rscript Run_correction_based_on_pseudo_ctrl.R

library(edgeR)

args=commandArgs(T)

#table to match the sample information
# AL_ATAC Ctrl    Li_F_adt
# AL_ATAC Ctrl    Li_M_adt
# BA_ATAC Ctrl    Li_F_adt
# BA_ATAC Ctrl    Li_M_adt
table<-read.table("match_LAB.txt")

data1<-read.table("LAB_RUVr_k3_pseudo.txt",header=T)
data2<-read.table("LAB_RUVr_k3.bed",header=T)
data<-cpm(data2+1)
ave<-cpm(data1+1)

names<-colnames(data)
n<-length(table$V1)
d<-c()
x<-c()

i<-1
a<-names[grep(as.character(table[i,1]),names)]
b<-a[grep(as.character(table[i,3]),a)]
LAB<-data[,b]
r<-b[grep("Ctrl",b)]
LAB_r<-LAB[,r]

a<-ave[,as.character(table[i,3])]
b<-rowSums(LAB_r)/ncol(LAB_r)
c<-as.data.frame(b/a)
colnames(c)<-"cor"
LAB_c<-(LAB)/c$cor
#BL_c<-(Blood+1)/c$cor
d<-LAB_c
#x<-BL_c

for(i in 2:n){

  a<-names[grep(as.character(table[i,1]),names)]
  b<-a[grep(as.character(table[i,3]),a)]
  LAB<-data[,b]
  r<-b[grep("Ctrl",b)]
  LAB_r<-LAB[,r]
  
  a<-ave[,as.character(table[i,3])]
  b<-rowSums(LAB_r)/ncol(LAB_r)
  c<-as.data.frame(b/a)
  colnames(c)<-"cor"
  LAB_c<-(LAB)/c$cor
  d<-cbind(d,LAB_c)

}

d<-round(d,4)
write.table(d,"Output_correct.txt",sep="\t",quote=F)

