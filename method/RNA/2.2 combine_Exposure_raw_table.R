# Rscript combine_Exposure_raw_table.R 
# Enter into the folder with the raw output of each exposure

file<-list.files(pattern="raw")
d<-c()
a1<-read.table(file[1],header=T)
d<-a1
for(i in 2:length(file)){
   a<-read.table(file[i],header=T)
   d<-cbind(d,a)    

}

write.table(d,"LAB_raw.bed",sep="\t",quote=F)

