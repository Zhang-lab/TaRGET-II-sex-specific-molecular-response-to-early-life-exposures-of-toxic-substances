library(tidyr)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(gridExtra)
library(RColorBrewer)

'%!in%' <- function(x,y)!('%in%'(x,y))

my_tissue<- "../liver_adt";
data_tissue<- read.table(paste0(my_tissue, "_homer/motif_gene.txt"), header=F)
colnames(data_tissue)<- c("Pvalue", "logPvalue", "%_of_Targets", "%_of_background", "Best_match_motif","Match_score", "Rank", "Gene","Group")

#>= 5% target
data_tissue1<- data_tissue[as.numeric(gsub("%", "", data_tissue[,"%_of_Targets"]))>=5,]
data_tissue2<- cbind.data.frame(data_tissue1, data.frame(t(data.frame(strsplit(data_tissue1$Group, "_"))))[,c(3,5)])
colnames(data_tissue2)[(ncol(data_tissue2)-1):ncol(data_tissue2)]<- c("lab", "expo");
data_tissue2$tag<- ifelse(data_tissue2$expo=="PM2.5", paste0(data_tissue2$expo, "_", data_tissue2$lab), data_tissue2$expo);
data_tissue2$TAG<- ifelse(grepl("_F_", data_tissue2$Group), paste0(data_tissue2$tag, "_F"), paste0(data_tissue2$tag, "_M"))

data_for_plot<- data.frame(matrix(nrow=9, ncol=7));
rownames(data_for_plot)<- c("As", "BPA10ug", "BPA10mg", "DEHP", "Pb", "PM2.5_BI", "PM2.5_MU", "TBT", "TCDD")
colnames(data_for_plot)<-  c("F_up", "F_dn", "M_up", "M_dn", "shared_up", "shared_dn", "inconsistent");

for(i in 1:nrow(data_for_plot)){
   data_expo<- data_tissue2[data_tissue2$tag== rownames(data_for_plot)[i],]
   data_F_up<- data_expo[grepl("_F$", data_expo$TAG) & grepl("MORE", data_expo$Group),]
   data_F_dn<- data_expo[grepl("_F$", data_expo$TAG) & grepl("LESS", data_expo$Group),]
   data_M_up<- data_expo[grepl("_M$", data_expo$TAG) & grepl("MORE", data_expo$Group),]
   data_M_dn<- data_expo[grepl("_M$", data_expo$TAG) & grepl("LESS", data_expo$Group),]

   data_for_plot[i, "shared_up"]<- length(which(data_F_up$Best_match_motif %in% data_M_up$Best_match_motif))
   data_for_plot[i, "shared_dn"]<- length(which(data_F_dn$Best_match_motif %in% data_M_dn$Best_match_motif))
   data_for_plot[i, "inconsistent"]<- length(which(data_F_up$Best_match_motif %in% data_M_dn$Best_match_motif)) + length(which(data_F_dn$Best_match_motif %in% data_M_up$Best_match_motif))

   data_for_plot[i, "F_up"]<- length(which(data_F_up$Best_match_motif %!in% c(data_M_up$Best_match_motif, data_M_dn$Best_match_motif)) )
   data_for_plot[i, "F_dn"]<- length(which(data_F_dn$Best_match_motif %!in% c(data_M_up$Best_match_motif, data_M_dn$Best_match_motif)) )
   data_for_plot[i, "M_up"]<- length(which(data_M_up$Best_match_motif %!in% c(data_F_up$Best_match_motif, data_F_dn$Best_match_motif)) )
   data_for_plot[i, "M_dn"]<- length(which(data_M_dn$Best_match_motif %!in% c(data_F_up$Best_match_motif, data_F_dn$Best_match_motif)) )
   }

data_for_plot1<- data.frame(data_for_plot %>% tibble::rownames_to_column(var="expo") %>% pivot_longer(-expo, names_to="type", values_to="count"))
data_for_plot1$count<- as.numeric(data_for_plot1$count)
data_for_plot1$expo<- factor(data_for_plot1$expo, levels=rev(c("As", "BPA10ug", "BPA10mg", "DEHP", "Pb", "PM2.5_BI", "PM2.5_MU", "TBT", "TCDD")))
data_for_plot1$type<- factor(data_for_plot1$type, levels=c("F_up", "F_dn", "M_up", "M_dn", "shared_up", "shared_dn", "inconsistent"));
count_max<- ceiling(max(data_for_plot1$count)/10)*10

p<- ggplot(data_for_plot1, aes(x=type, y=expo, color=type, size=count)) + geom_point()
p<- p+ scale_color_manual(values=c("F_up"="#E8ACB7", "F_dn"="#E16A86", "M_up"="#B2BAE5", "M_dn"="#768BE6", "shared_up"="#ADDFAD", "shared_dn"="#92a881", "inconsistent"="#DCDCDC"));
p<- p+ guides(color="none")
p<- p+ scale_size(limits=c(1, count_max))
p<- p+ theme_bw()
p<- p+ scale_x_discrete(guide = guide_axis(angle = 90))
p<- p+ ggtitle(paste0("ATAC motif liver adt"));
p<- p+ geom_vline(xintercept=c(2.5,4.5), color="black")
p<- p+ theme(panel.grid.major=element_blank())

pdf("dotplot_atac_lv_adt_motif.pdf", width=3.5, height=3);
print(p);
dev.off();

###### Piechart for uniq motif by sex
data_F<- unique(data_tissue2[grepl("_F_", data_tissue2$Group), c("Best_match_motif", "tag")])
data_M<- unique(data_tissue2[grepl("_M_", data_tissue2$Group), c("Best_match_motif", "tag")])

####
get_piechart<- function(my_data, my_sex){

data_freq<- data.frame(table(my_data$Best_match_motif))
multi_single_motif_vec<- data_freq[data_freq$Freq==1,"Var1"]
my_data2<- my_data[my_data$Best_match_motif %in% multi_single_motif_vec,]
data_for_plot<- data.frame(matrix(0, nrow=11, ncol=2));
colnames(data_for_plot)<- c("Freq", "label");
rownames(data_for_plot)<- c("As", "BPA10ug", "BPA10mg", "DEHP", "Pb", "PM2.5_BI", "PM2.5_MU", "TBT", "TCDD", "2", "3");
data_for_plot["2",1]<- nrow(data_freq[data_freq$Freq==2,])
data_for_plot["3",1]<- nrow(data_freq[data_freq$Freq==3,])
for(i in 1:9){
    data_for_plot[i, 1]<- nrow(my_data2[my_data2$tag== rownames(data_for_plot)[i],])   
   }
data_for_plot$Var1<- factor(rownames(data_for_plot), levels=c("As", "BPA10ug", "BPA10mg", "DEHP", "Pb", "PM2.5_BI", "PM2.5_MU", "TBT", "TCDD", "2", "3"));
data_for_plot$label<- paste0(data_for_plot$Var1, " n=", data_for_plot$Freq);
data_for_plot<- data_for_plot %>%
  mutate(csum = rev(cumsum(rev(Freq))),
         pos = Freq/2 + lead(csum, 1),
         pos = if_else(is.na(pos), Freq/2, pos))
data_for_plot[nrow(data_for_plot),"pos"]<- 1

p<- ggplot(data_for_plot, aes(x="", y=Freq, fill=Var1, color=Var1))+ geom_col(width=1, color="white")+ coord_polar(theta="y")
p<- p+ theme_void() + scale_fill_brewer(palette = "Set3")
p<- p+ ggtitle(my_sex);
p<- p+ geom_label_repel(aes(y=pos, label=label), size=4.5, nudge_x=1, show.legend=F, color="black")

return(p);
}
####

p1<- get_piechart(data_F, "uniq F");
p2<- get_piechart(data_M, "uniq M");

pdf("piechart_uniq_motif_bySex_atac_lv_adt.pdf", width=10, height=5);
grid.arrange(p1, p2, nrow=1);
dev.off();


