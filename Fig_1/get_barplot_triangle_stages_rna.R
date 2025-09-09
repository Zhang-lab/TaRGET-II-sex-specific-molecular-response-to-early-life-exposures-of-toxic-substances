.libPaths("/home/shuhua/R/x86_64-pc-linux-gnu-library/4.0")
library(ggplot2)
library(HH)
library(tidyverse)
library(scater)
library(tidyr)
library(ggbreak)

expo_vec<- c("As", "BPA_L", "BPA_H", "DEHP", "Pb", "PM2.5_Biswal", "PM2.5_Mutlu", "TBT", "TCDD")
stage_vec<- c("wl", "adt", "aged");

data_count<- data.frame(matrix(0, ncol=7, nrow=27));
colnames(data_count)<- c("uniq_F_up","uniq_F_dn","uniq_M_up","uniq_M_dn","shared_up","shared_dn","shared_inconsistent")
rownames(data_count)<- paste(rep(expo_vec, each=3), rep(stage_vec, 9), sep="_")
#####
fill_data<- function(my_file, my_stage){

my_data<- read.table(my_file, header=T, sep="\t", row.names=1)

for(i in 1:nrow(my_data)){
    data_count[paste0(rownames(my_data)[i], "_", my_stage), ]<<- my_data[i,]
   }

}
#####

fill_data("count_liver_weanling.txt", "wl")
fill_data("count_liver_adult.txt", "adt")
fill_data("count_liver_aged.txt", "aged")

data_count$M_count<- rowSums(data_count[,c("uniq_M_up", "uniq_M_dn", "shared_up", "shared_dn", "shared_inconsistent")])
data_count$F_count<- rowSums(data_count[,c("uniq_F_up", "uniq_F_dn", "shared_up", "shared_dn", "shared_inconsistent")])
data_count<- data_count[!is.na(data_count$uniq_F_up),]
data_count$LFC<- log2(ifelse(data_count$M_count!=0, data_count$M_count, 1)/ifelse(data_count$F_count!=0, data_count$F_count, 1) )

data_count$perc<- round(100*(data_count$shared_up + data_count$shared_dn + data_count$shared_inconsistent)/(data_count$M_count + data_count$F_count -data_count$shared_up -data_count$shared_dn - data_count$shared_inconsistent),2)
data_count$perc<- ifelse(is.na(data_count$perc), 0, data_count$perc)
data_count$condition<- factor(rownames(data_count), levels=rev(rownames(data_count)))
data_count$stage<- sapply(strsplit(as.character(data_count$condition), "_"), tail, 1)

q<-  ggplot(data_count, aes(x=condition, y=perc, fill=stage)) + geom_bar(stat="identity", alpha=0.5)
q<- q+ coord_flip()
q<- q+ theme(legend.position= "none", axis.text=element_text(size=11), axis.title=element_text(size=15), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.ticks.y=element_blank())
q<- q+ scale_fill_manual(values=c("wl"="#B9D5C2", "adt"="#66B58C", "aged"= "#0E9A60"))
q<- q+ geom_vline(xintercept=0.5+ seq(3,24,3), color="darkgray")
q<- q+  scale_y_continuous(expand = c(0, 0))
q<- q+ geom_text(data=data_count, aes(label=ifelse(perc!=0, perc, "NA")), y=7)
pdf("shared_sex_expoDEG_bar.pdf", width=2.5);
print(q);
dev.off();

write.table(data.frame("name"=rownames(data_count), data_count), file="shared_uniq_expoDEG_bySex_stages", quote=F, sep="\t", row.names=F, col.names=T)

data_count$uniq_F<- data_count$uniq_F_up + data_count$uniq_F_dn;
data_count$uniq_M<- data_count$uniq_M_up + data_count$uniq_M_dn;
data_count$shared_F<- data_count$shared_up + data_count$shared_dn + data_count$shared_inconsistent
data_count$shared_M<- data_count$shared_up + data_count$shared_dn + data_count$shared_inconsistent

data_for_text<- data.frame(data_count[,c("F_count", "M_count")] %>% tibble::rownames_to_column(var = "cond") %>% pivot_longer(-cond, names_to = "Name", values_to = "count") )
data_for_plot<- data.frame(data_count[,c("uniq_F", "uniq_M", "shared_F", "shared_M")] %>% tibble::rownames_to_column(var = "cond") %>% pivot_longer(-cond, names_to = "Name", values_to = "count") )

data_for_text$cond<- factor(data_for_text$cond, rev(rownames(data_count)));
data_for_plot$cond<- factor(data_for_plot$cond, rev(rownames(data_count)));
data_for_plot$col<- paste0(sapply(strsplit(as.character(data_for_plot$cond), "_"), tail,1), "_", data_for_plot$Name)
data_for_plot$col<- factor(data_for_plot$col, levels=c( "wl_uniq_F", "wl_uniq_M",  "wl_shared_F",  "wl_shared_M", "adt_uniq_F", "adt_uniq_M",  "adt_shared_F", "adt_shared_M", "aged_uniq_F", "aged_uniq_M", "aged_shared_F", "aged_shared_M"))

p<- ggplot(data_for_plot, aes(x=cond, y=ifelse(grepl("_M", Name), count, -count)))+ geom_bar(stat="identity",  aes(fill=col, color=col))
p<- p+ scale_fill_manual(values=c("wl_uniq_F"= alpha("#E8ACB7",0.4), "adt_uniq_F"="#E8ACB7", "aged_uniq_F"= "#E16A86", "wl_uniq_M"= alpha("#B2BAE5",0.4), "adt_uniq_M"="#B2BAE5", "aged_uniq_M"= "#768BE6", "wl_shared_F"="#B9D5C2", "wl_shared_M"="#B9D5C2", "adt_shared_F"="#66B58C", "adt_shared_M"="#66B58C", "aged_shared_F"="#0E9A60", "aged_shared_M"="#0E9A60") )
p<- p+ scale_color_manual(values=c("wl_shared_F"="#A9A9A9", "wl_shared_M"="#A9A9A9", "adt_shared_F"="#A9A9A9", "adt_shared_M"="#A9A9A9", "aged_shared_F"="#A9A9A9", "aged_shared_M"="#A9A9A9", "wl_uniq_F"="#A9A9A900", "wl_uniq_M"="#A9A9A900", "adt_uniq_F"="#A9A9A900", "adt_uniq_M"="#A9A9A900", "aged_uniq_F"="#A9A9A900", "aged_uniq_M"="#A9A9A900"));


scale_ratio<- max(data_for_plot$count)/max(abs(data_count$LFC))
p<- p+ geom_hline(yintercept=c((-1)*scale_ratio,scale_ratio),  color="gray", size=0.8, linetype="dotted");

data_count1<- data_count[abs(data_count$LFC)>=1,]
p<- p+ geom_point(data=data_count1, aes(x=rownames(data_count1), y=LFC*scale_ratio), pch=2, size=5, fill="darkgray")
second_y_limit<- max(abs(data_count1$LFC));
primary_y_limit<- 2800
print(second_y_limit);
interval<- 1;
if(ceiling(second_y_limit)>5){
   interval<- 2;
  }
if(ceiling(second_y_limit)<=1){
  interval<- 0.5;
  }

p<- p+ scale_y_continuous(name="DEG count", sec.axis=sec_axis(~./scale_ratio, name="log2(M_DEG/F_DEG)"), limits=c((-1)*primary_y_limit, primary_y_limit) )

p<- p+ geom_vline(xintercept=0.5+seq(3,24,3), color="darkgray") 
p<- p+ geom_hline(yintercept=0, color="#DCDCDC", linetype="dashed");
p<- p+ coord_flip() + theme_linedraw();
p<- p+ theme(legend.position= "none", axis.text=element_text(size=11), axis.title=element_text(size=15), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
p<- p+ labs(x="Exposure in stage")

p<- p+ geom_text(data=data_for_text[data_for_text$Name=="F_count",], aes(label=ifelse(count!=0, count, "NA"), x=cond), y=-2700)
p<- p+ geom_text(data=data_for_text[data_for_text$Name=="M_count",], aes(label=ifelse(count!=0, count, "NA"), x=cond), y=2700)


pdf("barplot_wTriangle_liver_rna.pdf");
p;
dev.off();




