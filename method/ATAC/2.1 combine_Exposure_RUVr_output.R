# Rscript combine_Exposure_RUVr_output.R
# Enter into the folder with the RUVr output of each exposure
file<-list.files(pattern="k3")
d<-c()
a1<-read.table(file[1],header=T)
d<-a1
for(i in 2:length(file)){
   a<-read.table(file[i],header=T)
   d<-cbind(d,a)    

}

write.table(d,"LAB_RUVr_k3.bed",sep="\t",quote=F)

