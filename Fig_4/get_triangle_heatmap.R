library(tidyr)
library(ggplot2)
library(RColorBrewer)

data_F<- read.table("GSEA_3stages_mh_F.txt", sep="\t", header=T)
data_M<- read.table("GSEA_3stages_mh_M.txt", sep="\t", header=T)

shared_F<- data_F[rownames(data_F) %in% rownames(data_M),]
shared_M<- data_M[rownames(data_M) %in% rownames(data_F),]

uniq_F<- data_F[!(rownames(data_F) %in% rownames(data_M)),]
uniq_M<- data_M[!(rownames(data_M) %in% rownames(data_F)),]
uniq_F<- uniq_F[, c(3:14,17:22,1:2,15:16,23,24)]

data_for_shared_F<- data.frame(shared_F %>% tibble::rownames_to_column(var="Row") %>% pivot_longer(-Row, names_to="cond", values_to="nes"))
data_for_shared_M<- data.frame(shared_M %>% tibble::rownames_to_column(var="Row") %>% pivot_longer(-Row, names_to="cond", values_to="nes"))
data_for_shared_F$x_pos<- match(data_for_shared_F$cond, colnames(uniq_F))
data_for_shared_F$y_pos<- as.numeric(as.factor(data_for_shared_F$Row)) +nrow(uniq_M);

list_out<- list();
list_count<- 1;
for(i in 1:nrow(data_for_shared_F)){
    list_out[[list_count]]<- c(data_for_shared_F[i, "Row"], data_for_shared_F[i, "cond"], data_for_shared_F[i, "nes"], "F", data_for_shared_F[i, "y_pos"]-0.5, data_for_shared_F[i, "x_pos"]-0.5);
    list_count<- list_count +1;
    list_out[[list_count]]<- c(data_for_shared_F[i, "Row"], data_for_shared_F[i, "cond"], data_for_shared_F[i, "nes"], "F", data_for_shared_F[i, "y_pos"]+0.5, data_for_shared_F[i, "x_pos"]-0.5);
    list_count<- list_count +1;
    list_out[[list_count]]<- c(data_for_shared_F[i, "Row"], data_for_shared_F[i, "cond"], data_for_shared_F[i, "nes"], "F", data_for_shared_F[i, "y_pos"]+0.5, data_for_shared_F[i, "x_pos"]+0.5);
    list_count<- list_count +1;
###
    list_out[[list_count]]<- c(data_for_shared_F[i, "Row"], data_for_shared_F[i, "cond"], shared_M[data_for_shared_F[i,"Row"], data_for_shared_F[i,"cond"]], "M", data_for_shared_F[i, "y_pos"]-0.5, data_for_shared_F[i, "x_pos"]-0.5);
    list_count<- list_count +1;
    list_out[[list_count]]<- c(data_for_shared_F[i, "Row"], data_for_shared_F[i, "cond"], shared_M[data_for_shared_F[i,"Row"], data_for_shared_F[i,"cond"]], "M", data_for_shared_F[i, "y_pos"]-0.5, data_for_shared_F[i, "x_pos"]+0.5);
    list_count<- list_count +1;
    list_out[[list_count]]<- c(data_for_shared_F[i, "Row"], data_for_shared_F[i, "cond"], shared_M[data_for_shared_F[i,"Row"], data_for_shared_F[i,"cond"]], "M", data_for_shared_F[i, "y_pos"]+0.5, data_for_shared_F[i, "x_pos"]+0.5);
    list_count<- list_count +1;
   }

data_for_plot_shared<- data.frame(do.call(rbind, list_out))
colnames(data_for_plot_shared)<- c("Row", "cond", "nes", "sex", "y_pos", "x_pos");
data_for_plot_shared$cond<- factor(data_for_plot_shared$cond, levels=c("BPA10ug_wl","BPA10ug_adt", "BPA10ug_aged", "BPA10mg_wl", "BPA10mg_adt", "BPA10mg_aged", "DEHP_wl", "DEHP_adt", "DEHP_aged", "Pb_wl", "Pb_adt", "Pb_aged", "PM2.5_Mutlu_wl", "PM2.5_Mutlu_adt", "PM2.5_Mutlu_aged","TBT_wl", "TBT_adt", "TBT_aged", "As_wl", "As_adt", "PM2.5_Biswal_wl", "PM2.5_Biswal_adt", "TCDD_wl", "TCDD_adt"))
data_for_plot_shared$nes<- as.numeric(as.character(data_for_plot_shared$nes));
data_for_plot_shared[data_for_plot_shared$nes==0,"nes"]<- NA;
data_for_plot_shared$x_pos<- as.numeric(as.character(data_for_plot_shared$x_pos));
data_for_plot_shared$y_pos<- as.numeric(as.character(data_for_plot_shared$y_pos));
data_for_plot_shared$Row<- factor(data_for_plot_shared$Row, levels=rownames(shared_F));

data_for_plot_F<- data.frame(uniq_F %>% tibble::rownames_to_column(var="Row") %>% pivot_longer(-Row, names_to="cond", values_to="nes"))
data_for_plot_F$sex<- "F";
data_for_plot_M<- data.frame(uniq_M %>% tibble::rownames_to_column(var="Row") %>% pivot_longer(-Row, names_to="cond", values_to="nes"))
data_for_plot_M$sex<- "M";

data_for_plot_F$y_pos<- as.numeric(as.factor(data_for_plot_F$Row)) + nrow(uniq_M) +4;
data_for_plot_F$x_pos<- match(data_for_plot_F$cond, colnames(uniq_F))

data_for_plot_M$y_pos<- as.numeric(as.factor(data_for_plot_M$Row))
data_for_plot_M$x_pos<- match(data_for_plot_M$cond, colnames(uniq_M))

data_for_plot_combined<- rbind.data.frame(data_for_plot_F, data_for_plot_M);
data_for_plot_combined$cond<- factor(data_for_plot_combined$cond, levels=colnames(uniq_F));
data_for_plot_combined[data_for_plot_combined$nes==0, "nes"]<- NA
data_for_plot_combined$nes<- as.numeric(as.character(data_for_plot_combined$nes))

data_for_plot_combined$cond<- factor(data_for_plot_combined$cond, levels=c("BPA10ug_wl","BPA10ug_adt", "BPA10ug_aged", "BPA10mg_wl", "BPA10mg_adt", "BPA10mg_aged", "DEHP_wl", "DEHP_adt", "DEHP_aged", "Pb_wl", "Pb_adt", "Pb_aged", "PM2.5_Mutlu_wl", "PM2.5_Mutlu_adt", "PM2.5_Mutlu_aged","TBT_wl", "TBT_adt", "TBT_aged", "As_wl", "As_adt", "PM2.5_Biswal_wl", "PM2.5_Biswal_adt", "TCDD_wl", "TCDD_adt"))

p<- ggplot(data_for_plot_combined, aes(y=y_pos, x=cond, fill=nes)) + geom_tile(color="darkgray", size=0.5) 
p<- p+ geom_polygon(data=data_for_plot_shared, aes(x=x_pos, y=y_pos, fill=nes, group=interaction(Row, cond,sex)), color="darkgray")
p<- p+ scale_fill_gradientn(colors=rev(brewer.pal(11, 'RdBu')), limits=c(-3.3, 3.3), na.value=brewer.pal(n = 3, name = "RdBu")[2])
p<- p+ theme_classic()
p<- p+ theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1))
p<- p+ geom_vline(xintercept=c(seq(3,15,3)+0.5, 20.5,22.5), color="white", size=1.2)
p<- p+ geom_vline(xintercept=c(18.5), color="black", size=1)
p<- p+ geom_hline(yintercept=c(5,9)+0.5, color="white", size=1.2)

dat_label<- unique(data_for_plot_combined[order(data_for_plot_combined$y_pos),c("Row", "y_pos")])
dat_label_shared<- unique(data_for_plot_shared[order(data_for_plot_shared$y_pos),c("Row", "y_pos")])
for_label<- c(dat_label[1:5,"Row"], as.character(unique(dat_label_shared$Row)), dat_label[6:nrow(dat_label),"Row"])
p<- p+ scale_y_continuous(expand=c(0,0), breaks=1:length(for_label), labels=for_label);
p<- p+ scale_x_discrete(expand=c(0,0));
p<- p+ ylab("M_uniq             shared                F_uniq");

pdf("GSEA_3stages_combined_sex_heatmap.pdf", width=9, height=6);
print(p);
dev.off();

