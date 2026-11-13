#Rscript build_pseudo_control.R
#remove the 10% low and 10% high value
data<-read.table("LAB_RUVr_k3.bed",header=T)
ave1<-c()
ave2<-c()
n<-colnames(data)
n<-n[grep("Ctrl",n)]
n1<-n[grep("Li_F",n)]
n2<-n[grep("Li_M",n)]
d1<-data[,n1]
j<-length(d1[,1])
k<-ncol(d1)
k1<-round(k*0.1)
k2<-round(k*0.9)
s1<-c()
for(i in 1:j){
  l<-as.numeric(d1[i,])
  m<-l[order(l,decreasing = TRUE)][k1:k2]
  s1[i]<-mean(m)
}

d2<-data[,n2]
j<-length(d2[,1])
k<-ncol(d2)
k1<-round(k*0.1)
k2<-round(k*0.9)
s2<-c()
for(i in 1:j){
  l<-as.numeric(d2[i,])
  m<-l[order(l,decreasing = TRUE)][k1:k2]
  s2[i]<-mean(m)
}

s<-as.data.frame(cbind(s1,s2))
colnames(s)<-c("Li_F_adt","Li_M_adt")
rownames(s)<-rownames(data)

write.table(s,"LAB_RUVr_k3_pseudo.txt",sep="\t",quote = F)


