library(ggplot2)
library(pheatmap)
library(scales)
library(tidyverse)
library(ggrepel)
library(plyr)

adt_DEG_list<- list.files("../treatment_DEG_list_noChrXY/adult_DEG/")
## liver only
adt_DEG_list<- adt_DEG_list[grep("liver", adt_DEG_list)]


adt_DEG_list_up<- adt_DEG_list[grep("_up$", adt_DEG_list)]
adt_DEG_list_dn<- adt_DEG_list[grep("_down$", adt_DEG_list)]

######
get_piechart<- function(my_data, my_tag){

my_data[my_data$Freq>=2,"cond"]<- "ge2"
cond_freq<- data.frame(table(as.character(my_data$cond)))
cond_freq$label<- paste0(cond_freq$Var1, " n=", cond_freq$Freq)
cond_freq$Var1<- factor(cond_freq$Var1, levels=c("As", "BPA10ug", "BPA10mg", "DEHP", "Pb", "PM2.5_Biswal", "PM2.5_Mutlu", "TBT", "TCDD", "ge2"))
cond_freq<- cond_freq[order(cond_freq$Var1),]
cond_freq<- cond_freq %>%
  mutate(csum = rev(cumsum(rev(Freq))),
         pos = Freq/2 + lead(csum, 1),
         pos = if_else(is.na(pos), Freq/2, pos))
cond_freq[nrow(cond_freq),"pos"]<- 1
s<- ggplot(cond_freq, aes(x="", y=Freq, fill=Var1)) + geom_col(width=1, color="white") + coord_polar(theta="y")
s<- s+ theme_void() + scale_fill_brewer(palette = "Set3")
s<- s+ ggtitle(paste0("n=", nrow(my_data)));
s<- s+ theme(legend.position="none",   plot.title = element_text(hjust = 0.5))
pdf(paste0("pie_multiExpo_", my_tag, "_expoDEG.pdf"));
print(s);
dev.off();

my_data4<- my_data[my_data$Freq>=4,]
my_data4<- my_data4[order(-my_data4$Freq),]

out_data<- data.frame(matrix(0, nrow=nrow(my_data4), ncol=18))
colnames(out_data)<- paste0(rep(c("As", "BPA10ug", "BPA10mg", "DEHP", "Pb", "PM2.5_Biswal", "PM2.5_Mutlu", "TBT", "TCDD"), each=2), rep(c("_F", "_M"), 9))
rownames(out_data)<- my_data4$Var1;

for(i in 1:length(adt_DEG_list)){
    data<- read.table(paste0("../treatment_DEG_list_noChrXY/adult_DEG/", adt_DEG_list[i]), header=T, sep="\t");
    expo<- strsplit(adt_DEG_list[i], "_")[[1]][1]
    if(grepl("Biswal", adt_DEG_list[i])){
       expo<- paste0(expo, "_Biswal");
      }
    if(grepl("Mutlu", adt_DEG_list[i])){
       expo<- paste0(expo, "_Mutlu");
      }
    if(grepl("Female", adt_DEG_list[i])){
       expo<- paste0(expo, "_F");
      }else{
       expo<- paste0(expo, "_M");
      }
    pos<- grep(expo, colnames(out_data));
print(pos);

    for(j in 1:nrow(data)){
        if(data[j, "gene"] %in% rownames(out_data)){
           out_data[grep(data[j, "gene"], rownames(out_data)),pos]<- data[j, "log2FoldChange"];
          }
       }
   }

data_for_plot<- out_data;

paletteLength <- 31
myColor <- colorRampPalette(c("navy", "white", "firebrick3"))(paletteLength)
myBreaks <- c(seq(min(data_for_plot), 0, length.out=ceiling(paletteLength/2) + 1),
          seq(max(data_for_plot)/paletteLength, max(data_for_plot), length.out=floor(paletteLength/2)))

exposure<- gsub("_F$|_M$", "", colnames(data_for_plot))
sex<- ifelse(grepl("_F", colnames(data_for_plot)), "F", "M")

ann<- cbind.data.frame(sex, exposure);
rownames(ann)<-  colnames(data_for_plot)
ann$exposure<- factor(ann$exposure, levels=c("As", "BPA10ug", "BPA10mg", "DEHP", "Pb", "PM2.5_Biswal", "PM2.5_Mutlu", "TBT", "TCDD"));

frequency<- data.frame(my_data4[,"Freq"])
rownames(frequency)<-  rownames(data_for_plot)
colnames(frequency)<- ("Frequency")

freq_data<- data.frame(table(frequency$Frequency))
for_gap_row<- cumsum(rev(freq_data$Freq))
for_gap_row<- for_gap_row[1:(length(for_gap_row)-1)]

pdf(paste0("heatmap_freq_", my_tag, "_DEG_expo_LFC.pdf"), width=5.9, height=6);
pheatmap(data_for_plot, cluster_cols=F, cluster_rows=F, color=myColor, gaps_col=2*(1:((ncol(data_for_plot)-2)/2)),  annotation_col=ann, annotation_row=frequency, main=paste0("LFC of frequent DEGs n=", nrow(data_for_plot)), gaps_row= for_gap_row, breaks=myBreaks)
dev.off();
}
#####
get_count<- function(my_list){

gene_vec<- c();
my_list_F_dn<- my_list[intersect(grep("Female", my_list), grep("_down$", my_list))]

F_specific_list<- list();
M_specific_list<- list();
F_specific_count<- 1;
M_specific_count<- 1;

expo_gene_count<- data.frame(matrix(nrow=9, ncol=2))

for(i in 1:length(my_list_F_dn)){
    data_F_dn<- read.table(paste0("../treatment_DEG_list_noChrXY/adult_DEG/", my_list_F_dn[i]), header=T, sep="\t");
    gene_F_dn<- data_F_dn$gene;
    data_F_up<- read.table(paste0("../treatment_DEG_list_noChrXY/adult_DEG/", gsub("_down", "_up", my_list_F_dn[i])), header=T, sep="\t");
    gene_F_up<- data_F_up$gene;

    data_M_dn<- read.table(paste0("../treatment_DEG_list_noChrXY/adult_DEG/", gsub("Female", "Male", my_list_F_dn[i])), header=T, sep="\t");
    gene_M_dn<- data_M_dn$gene;
    data_M_up<- read.table(paste0("../treatment_DEG_list_noChrXY/adult_DEG/", gsub("_down", "_up", gsub("Female", "Male", my_list_F_dn[i]))), header=T, sep="\t");
    gene_M_up<- data_M_up$gene;

    shared_genes<- intersect(c(gene_F_up, gene_F_dn), c(gene_M_up, gene_M_dn))
    gene_vec<- c(gene_vec, shared_genes);

    expo<- strsplit(my_list_F_dn[i], "_")[[1]][1]
    lab<- strsplit(my_list_F_dn[i], "_")[[1]][7]
    if(lab %in%  c("Biswal", "Mutlu")){
       expo<- paste0(expo, "_", lab);
      }
    gene_F<- c(gene_F_up, gene_F_dn);
    gene_M<- c(gene_M_up, gene_M_dn);
    expo_gene_count[i,]<- c(expo, length(unique(c(gene_F, gene_M))));

    F_specific<- gene_F[!(gene_F %in% gene_M)]
    M_specific<- gene_M[!(gene_M %in% gene_F)]
    if(length(F_specific)>0){
       for(j in 1:length(F_specific)){
           F_specific_list[[F_specific_count]]<- c(F_specific[j], expo);
           F_specific_count<- F_specific_count +1;
          }
      }
    if(length(M_specific)>0){
       for(j in 1:length(M_specific)){
           M_specific_list[[M_specific_count]]<- c(M_specific[j], expo);
           M_specific_count<- M_specific_count +1;
          }
      }
   }
data_F_specific<-  data.frame(do.call(rbind, F_specific_list))
data_M_specific<-  data.frame(do.call(rbind, M_specific_list))

info_F_specific<- data.frame(ddply(data_F_specific, .(X1), summarise, info=list(X2)))
info_M_specific<- data.frame(ddply(data_M_specific, .(X1), summarise, info=list(X2)))

freq_F_specific<- data.frame(table(data_F_specific[,1]))
freq_M_specific<- data.frame(table(data_M_specific[,1]))

freq_F_specific$cond<- freq_F_specific$Freq;
freq_M_specific$cond<- freq_M_specific$Freq;

freq_F_specific$cond<- ifelse(freq_F_specific$Freq==1, info_F_specific[match(freq_F_specific$Var1, info_F_specific$X1),"info"],   freq_F_specific$cond)
freq_M_specific$cond<- ifelse(freq_M_specific$Freq==1, info_M_specific[match(freq_M_specific$Var1, info_M_specific$X1),"info"],   freq_M_specific$cond)

get_piechart(freq_F_specific, "Fspecific");
get_piechart(freq_M_specific, "Mspecific");

gene_freq<- data.frame(table(gene_vec))
gene_vec<- unique(gene_vec);

expo_vec<- sapply(strsplit(my_list, "_"), head, 1)
expo_vec<- ifelse(grepl("Mutlu", my_list), paste0(expo_vec,"_Mutlu"), expo_vec)
expo_vec<- ifelse(grepl("Biswal", my_list), paste0(expo_vec,"_Biswal"), expo_vec)
expo_vec<- ifelse(grepl("Female", my_list), paste0(expo_vec,"_F"), paste0(expo_vec, "_M"))
expo_vec<- unique(expo_vec)

out_data<- data.frame(matrix(0, nrow=length(gene_vec), ncol=2*length(my_list_F_dn)));
colnames(out_data)<- expo_vec;
rownames(out_data)<- gene_vec;
data_cond<- data.frame(matrix("", nrow=length(gene_vec), 3));
rownames(data_cond)<- gene_vec;
data_cond[,1]<- gene_vec;

for(i in 1:length(my_list)){
    data<- read.table(paste0("../treatment_DEG_list_noChrXY/adult_DEG/", my_list[i]), header=T, sep="\t");
    expo<- strsplit(my_list[i], "_")[[1]][1]
    if(grepl("Biswal", my_list[i])){
       expo<- paste0(expo, "_Biswal");
      }
    if(grepl("Mutlu", my_list[i])){
       expo<- paste0(expo, "_Mutlu");
      }
    if(grepl("Female", my_list[i])){
       expo<- paste0(expo, "_F");
      }else{
       expo<- paste0(expo, "_M");
      }
    pos<- grep(expo, colnames(out_data));
print(pos);

    for(j in 1:nrow(data)){
        if(data[j, "gene"] %in% rownames(out_data)){
           out_data[grep(data[j, "gene"], rownames(out_data)),pos]<- data[j, "log2FoldChange"];
           data_cond[data[j,"gene"],2]<- paste0(data_cond[data[j,"gene"],2], ",", gsub("_F$|_M$", "", expo));
          }
       }
   }
out_data$sum<- gene_freq[match(gene_vec, gene_freq$gene_vec),"Freq"]

out_data$expo<- data_cond[match(rownames(out_data), row.names(data_cond)),"X2"]
out_data$cond<- out_data$expo
for(i in 1:nrow(out_data)){
    if(out_data[i, "sum"]==1){
       out_data[i, "cond"]<- strsplit(out_data[i, "expo"], ",")[[1]][duplicated(strsplit(out_data[i, "expo"], ",")[[1]])]
      }else{
       out_data[i, "cond"]<- "ge2"
      }
   }
cond_freq<- data.frame(table(out_data$cond))
cond_freq$label<- paste0(cond_freq$Var1, " n=", cond_freq$Freq)
cond_freq$Var1<- factor(cond_freq$Var1, levels=c("As", "BPA10ug", "BPA10mg", "DEHP", "Pb", "PM2.5_Biswal", "PM2.5_Mutlu", "TBT", "TCDD", "ge2"))
cond_freq<- cond_freq[order(cond_freq$Var1),]
cond_freq<- cond_freq %>% 
  mutate(csum = rev(cumsum(rev(Freq))), 
         pos = Freq/2 + lead(csum, 1),
         pos = if_else(is.na(pos), Freq/2, pos))
cond_freq[nrow(cond_freq),"pos"]<- 1
s<- ggplot(cond_freq, aes(x="", y=Freq, fill=Var1)) + geom_col(width=1, color="white") + coord_polar(theta="y")
s<- s+ theme_void() + scale_fill_brewer(palette = "Set3")
s<- s+ ggtitle(paste0("n=", nrow(out_data)));
s<- s+ theme(plot.title = element_text(hjust = 0.5))
pdf("pie_multiExpo_shared_expoDEG.pdf", width=8);
print(s);
dev.off();

data_cond1<- data_cond[rownames(data_cond) %in% gene_freq[gene_freq$Freq==1,"gene_vec"],]
data_cond2<- data_cond1;
for(i in 1:nrow(data_cond1)){
    data_cond2[i,2]<- strsplit(data_cond1[i,2], ",")[[1]][duplicated(strsplit(data_cond1[i,2], ",")[[1]])]
    LFC_F<- out_data[rownames(data_cond2)[i], paste0(data_cond2[i,2], "_F")]
    LFC_M<- out_data[rownames(data_cond2)[i], paste0(data_cond2[i,2], "_M")]
    if(LFC_F * LFC_M <0){ 
       data_cond2[i,3]<- "inconsistent";
      }else{
       if(LFC_F >0){
          data_cond2[i,3]<- "shared_up";
         }else{
          data_cond2[i,3]<- "shared_dn";
         }
      }
   }
write.table(data_cond2, file="shared_DEG_single_expo.txt", quote=F, sep="\t", row.names=F, col.names=F);
colnames(data_cond2)<- c("gene", "expo", "type");

data_cond3<- data.frame(matrix(nrow=27, ncol=3));
colnames(data_cond3)<- c("expo", "type", "count");
data_cond3[,"expo"]<- rep(c("As","BPA10mg","BPA10ug","DEHP","Pb","PM2.5_Biswal","PM2.5_Mutlu","TBT","TCDD"), each=3)
data_cond3[,"type"]<- rep(c("shared_up", "shared_dn", "inconsistent"), 9)
for(i in 1:nrow(data_cond3)){
    
    data_cond3[i, "count"]<- nrow(data_cond2[data_cond2$expo==data_cond3[i,"expo"] & data_cond2$type== data_cond3[i, "type"],])
   }

data_cond3<- ddply(data_cond3, .(expo), transform, percent=round(count/sum(count)*100,2))
data_cond3$expo<- factor(data_cond3$expo, levels=rev(c("As","BPA10ug","BPA10mg","DEHP","Pb","PM2.5_Biswal","PM2.5_Mutlu","TBT","TCDD")))
expo_count<-aggregate(data_cond3$count, by=list(count=data_cond3$expo), FUN=sum)
data_cond3$type<- factor(data_cond3$type, levels=c( "shared_up", "shared_dn", "inconsistent"));
expo_count$total<- expo_gene_count[match(expo_count$count, expo_gene_count$X1),"X2"]


q<- q+ coord_flip()
q<- q+ scale_y_continuous(expand = c(0, 0))

q<- q+geom_text(aes(count, 0.75, label=paste0(x, "/", total), fill=NULL), data=expo_count)
q<- q+geom_text(aes(expo, 0.2, label=ifelse(type=="inconsistent", paste0(percent, "%"), ""), fill=NULL), data=data_cond3)


q<- q+ scale_fill_manual(values=c("shared_up"="#ADDFAD", "shared_dn"="#92a881", "inconsistent"="#DCDCDC"));
q<- q+ ylab("ratio")

pdf("bar_shared_expoDEG_singleExpo_direction.pdf", height=2.2, width=4);
print(q);
dev.off();
return(out_data);
}
######

adt_count<- get_count(adt_DEG_list)
adt_count<- adt_count[order(-adt_count$sum),]

count_freq<-  data.frame(table(adt_count$sum))
count_freq$Var1<- factor(count_freq$Var1, levels= rev(count_freq$Var1))
p<- ggplot(count_freq, aes(x=Var1, y=Freq) )+geom_bar(stat="identity", fill=alpha(rainbow(nrow(count_freq) ),0.2) )
p<- p+ theme_classic()
p<- p + scale_y_continuous(expand = c(0, 0))
p<- p+ ylab("DEG count") + xlab("Number of exposoure+sex")
p<- p+ ggtitle("Frequency of shared DEGs among exposures");
p<- p+ coord_flip()
p<- p+ annotate("text", size=5, y=100, label=count_freq$Freq, x=count_freq$Var1);

pdf("freq_shared_DEG_among_expo.pdf", height=3);
print(p);
dev.off();

freq_cutoff<- 2
data_for_plot<- adt_count[adt_count$sum >= freq_cutoff, -ncol(adt_count)]
data_for_plot<- data_for_plot[, paste0(rep(c("As", "BPA10ug", "BPA10mg", "DEHP", "Pb", "PM2.5_Biswal", "PM2.5_Mutlu", "TBT", "TCDD"), each=2), rep(c("_F", "_M"), 9)) ]

paletteLength <- 31
myColor <- colorRampPalette(c("navy", "white", "firebrick3"))(paletteLength)
myBreaks <- c(seq(min(data_for_plot), 0, length.out=ceiling(paletteLength/2) + 1),
          seq(max(data_for_plot)/paletteLength, max(data_for_plot), length.out=floor(paletteLength/2)))

exposure<- gsub("_F$|_M$", "", colnames(data_for_plot))
sex<- ifelse(grepl("_F", colnames(data_for_plot)), "F", "M")

ann<- cbind.data.frame(sex, exposure);
rownames(ann)<-  colnames(data_for_plot)
ann$exposure<- factor(ann$exposure, levels=c("As", "BPA10ug", "BPA10mg", "DEHP", "Pb", "PM2.5_Biswal", "PM2.5_Mutlu", "TBT", "TCDD"));

frequency<- data.frame(adt_count[adt_count$sum >= freq_cutoff,"sum"])
rownames(frequency)<-  rownames(data_for_plot)
colnames(frequency)<- ("Frequency")

freq_data<- data.frame(table(frequency$Frequency))
for_gap_row<- cumsum(rev(freq_data$Freq))
for_gap_row<- for_gap_row[1:(length(for_gap_row)-1)]

pdf("heatmap_freq_shared_DEG_expo_LFC.pdf", width=5.9, height=8);
pheatmap(data_for_plot, cluster_cols=F, cluster_rows=F, color=myColor, gaps_col=2*(1:((ncol(data_for_plot)-2)/2)),  annotation_col=ann, annotation_row=frequency, main=paste0("LFC of frequent DEGs n=", nrow(data_for_plot)), gaps_row= for_gap_row, breaks=myBreaks)
dev.off();











