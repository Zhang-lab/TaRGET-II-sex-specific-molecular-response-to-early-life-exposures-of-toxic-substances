library(tidyr)
library(dplyr)
library(ggplot2)
library(scales)

'%!in%' <- function(x,y)!('%in%'(x,y))

all_wl_up<- read.table("../../FvsM_all_list/df_liver_wl_up", header=T, sep="\t", row.names=1);
all_wl_dn<- read.table("../../FvsM_all_list/df_liver_wl_dn", header=T, sep="\t", row.names=1);
all_adt_up<- read.table("../../FvsM_all_list/df_liver_adt_up", header=T, sep="\t", row.names=1);
all_adt_dn<- read.table("../../FvsM_all_list/df_liver_adt_dn", header=T, sep="\t", row.names=1);
all_aged_up<- read.table("../../FvsM_all_list/df_liver_aged_up", header=T, sep="\t", row.names=1);
all_aged_dn<- read.table("../../FvsM_all_list/df_liver_aged_dn", header=T, sep="\t", row.names=1);

combined_genes<- union(union(rownames(rbind.data.frame(all_wl_up, all_wl_dn)), rownames(rbind.data.frame(all_adt_up, all_adt_dn))), rownames(rbind.data.frame(all_aged_up, all_aged_dn) ) )

all_CPM<- data.frame(matrix(0, nrow=length(combined_genes), ncol=6));
rownames(all_CPM)<- combined_genes;
colnames(all_CPM)<- c("wl_F", "wl_M", "adt_F", "adt_M", "aged_F", "aged_M");

all_CPM[match(rownames(all_wl_up), rownames(all_CPM)),"wl_F"]<- all_wl_up$CPM_F
all_CPM[match(rownames(all_wl_up), rownames(all_CPM)),"wl_M"]<- all_wl_up$CPM_M
all_CPM[match(rownames(all_wl_dn), rownames(all_CPM)),"wl_F"]<- all_wl_dn$CPM_F
all_CPM[match(rownames(all_wl_dn), rownames(all_CPM)),"wl_M"]<- all_wl_dn$CPM_M

all_CPM[match(rownames(all_adt_up), rownames(all_CPM)),"adt_F"]<- all_adt_up$CPM_F
all_CPM[match(rownames(all_adt_up), rownames(all_CPM)),"adt_M"]<- all_adt_up$CPM_M
all_CPM[match(rownames(all_adt_dn), rownames(all_CPM)),"adt_F"]<- all_adt_dn$CPM_F
all_CPM[match(rownames(all_adt_dn), rownames(all_CPM)),"adt_M"]<- all_adt_dn$CPM_M

all_CPM[match(rownames(all_aged_up), rownames(all_CPM)),"aged_F"]<- all_aged_up$CPM_F
all_CPM[match(rownames(all_aged_up), rownames(all_CPM)),"aged_M"]<- all_aged_up$CPM_M
all_CPM[match(rownames(all_aged_dn), rownames(all_CPM)),"aged_F"]<- all_aged_dn$CPM_F
all_CPM[match(rownames(all_aged_dn), rownames(all_CPM)),"aged_M"]<- all_aged_dn$CPM_M

#######
data_wl_up<- read.table("../../FvsM_DEG_list_noChrXY/deg_liver_wl_up", header=T, sep="\t", row.names=1)
data_wl_dn<- read.table("../../FvsM_DEG_list_noChrXY/deg_liver_wl_dn", header=T, sep="\t", row.names=1)

data_adt_up<- read.table("../../FvsM_DEG_list_noChrXY/deg_liver_adt_up", header=T, sep="\t", row.names=1)
data_adt_dn<- read.table("../../FvsM_DEG_list_noChrXY/deg_liver_adt_dn", header=T, sep="\t", row.names=1)

data_aged_up<- read.table("../../FvsM_DEG_list_noChrXY/deg_liver_aged_up", header=T, sep="\t", row.names=1)
data_aged_dn<- read.table("../../FvsM_DEG_list_noChrXY/deg_liver_aged_dn", header=T, sep="\t", row.names=1)

deg_wl_up<- rownames(data_wl_up)
deg_adt_up<- rownames(data_adt_up)
deg_aged_up<- rownames(data_aged_up)

deg_wl_dn<- rownames(data_wl_dn)
deg_adt_dn<- rownames(data_adt_dn)
deg_aged_dn<- rownames(data_aged_dn)
########
all_CPM_adt_Fdom<- all_CPM[deg_adt_up,]
all_CPM_adt_Mdom<- all_CPM[deg_adt_dn,]
all_CPM_adt_Fdom<- all_CPM_adt_Fdom[order(all_CPM_adt_Fdom$adt_F),]
all_CPM_adt_Mdom<- all_CPM_adt_Mdom[order(all_CPM_adt_Mdom$adt_M),]

all_CPM_wl_Fdom<- all_CPM[deg_wl_up,]
all_CPM_wl_Mdom<- all_CPM[deg_wl_dn,]
all_CPM_wl_Fdom<- all_CPM_wl_Fdom[order(all_CPM_wl_Fdom$wl_F),]
all_CPM_wl_Mdom<- all_CPM_wl_Mdom[order(all_CPM_wl_Mdom$wl_M),]

all_CPM_aged_Fdom<- all_CPM[deg_aged_up,]
all_CPM_aged_Mdom<- all_CPM[deg_aged_dn,]
all_CPM_aged_Fdom<- all_CPM_aged_Fdom[order(all_CPM_aged_Fdom$aged_F),]
all_CPM_aged_Mdom<- all_CPM_aged_Mdom[order(all_CPM_aged_Mdom$aged_M),]

####################

CPM_lv_adt<- read.table("../../CPM_list/liver_adt_pseudo_corrected", header=T, row.names=1, sep="\t")
CPM_lv_wl<- read.table("../../CPM_list/liver_wl_pseudo_corrected", header=T, row.names=1, sep="\t")
CPM_lv_aged<- read.table("../../CPM_list/liver_aged_pseudo_corrected", header=T, row.names=1, sep="\t")

expo_files_adt<- list.files("../../treatment_DEG_list/adult_DEG/");
expo_files_adt_liver<- expo_files_adt[grepl("liver", expo_files_adt)]
expo_files_adt_blood<- expo_files_adt[grepl("blood", expo_files_adt)]
expo_files_adt_lung<- expo_files_adt[grepl("lung", expo_files_adt)]
expo_files_adt_brain<- expo_files_adt[grepl("brain", expo_files_adt)]

expo_files_wl<- list.files("../../treatment_DEG_list/weanling_DEG/");
expo_files_wl_liver<- expo_files_wl[grepl("liver", expo_files_wl)]
expo_files_wl_blood<- expo_files_wl[grepl("blood", expo_files_wl)]
expo_files_wl_lung<- expo_files_wl[grepl("lung", expo_files_wl)]
expo_files_wl_heart<- expo_files_wl[grepl("heart", expo_files_wl)]

expo_files_aged<- list.files("../../treatment_DEG_list/aged_DEG/");
expo_files_aged_liver<- expo_files_aged[grepl("liver", expo_files_aged)]
#######

#######
count_expo<- function(my_data, my_expo_files, my_path, my_pie_title){

expo_count<- rep(0, nrow(my_data));
expo_info<- rep("", nrow(my_data));
for(i in 1:length(my_expo_files)){

    data_expo<- read.table(paste0(my_path, my_expo_files[i]), header=T, sep="\t");
    expo_tag<- paste0(strsplit(my_expo_files[i], "_")[[1]][1], "_", strsplit(my_expo_files[i], "_")[[1]][7])
    for(j in 1:nrow(my_data)){
        if(rownames(my_data)[j] %in% data_expo$gene){
           expo_count[j]<- expo_count[j] +1;
           expo_info[j]<- paste0(expo_info[j], expo_tag, ",");
          }
       }
   }
overlap_detail<- cbind.data.frame(my_data, expo_count, expo_info)
write.table(data.frame("name"=rownames(overlap_detail), overlap_detail), file=paste0("info_", my_pie_title, ".txt"), quote=F, sep="\t", row.names=F, col.names=T);

out_name<- paste0("pie_", my_pie_title, ".pdf");
pdf(out_name);
print(pie(table(expo_count), main=paste0(my_pie_title, " overlap to expo"), labels=paste0(names(table(expo_count)), " n=", as.character(table(expo_count)))));
dev.off();
expo_count[expo_count>3]<- 3
return(expo_count);
}
##########
expo_count_adt_lv_Fdom_F<- count_expo(all_CPM_adt_Fdom, expo_files_adt_liver[grep("_Female_", expo_files_adt_liver)],"../../treatment_DEG_list/adult_DEG/","lv_adt_DEG_Fdom_Fexpo");
expo_count_adt_lv_Fdom_M<- count_expo(all_CPM_adt_Fdom, expo_files_adt_liver[grep("_Male_", expo_files_adt_liver)], "../../treatment_DEG_list/adult_DEG/", "lv_adt_DEG_Fdom_Mexpo");

expo_count_adt_lv_Mdom_F<- count_expo(all_CPM_adt_Mdom, expo_files_adt_liver[grep("_Female_", expo_files_adt_liver)],"../../treatment_DEG_list/adult_DEG/","lv_adt_DEG_Mdom_Fexpo");
expo_count_adt_lv_Mdom_M<- count_expo(all_CPM_adt_Mdom, expo_files_adt_liver[grep("_Male_", expo_files_adt_liver)],"../../treatment_DEG_list/adult_DEG/", "lv_adt_DEG_Mdom_Mexpo");

expo_count_wl_lv_Fdom_F<- count_expo(all_CPM_wl_Fdom, expo_files_wl_liver[grep("_Female_", expo_files_wl_liver)], "../../treatment_DEG_list/weanling_DEG/", "lv_wl_DEG_Fdom_F");
expo_count_wl_lv_Fdom_M<- count_expo(all_CPM_wl_Fdom, expo_files_wl_liver[grep("_Male_", expo_files_wl_liver)], "../../treatment_DEG_list/weanling_DEG/", "lv_wl_DEG_Fdom_M");

expo_count_wl_lv_Mdom_F<- count_expo(all_CPM_wl_Mdom, expo_files_wl_liver[grep("_Female_", expo_files_wl_liver)], "../../treatment_DEG_list/weanling_DEG/", "lv_wl_DEG_Mdom_F");
expo_count_wl_lv_Mdom_M<- count_expo(all_CPM_wl_Mdom, expo_files_wl_liver[grep("_Male_", expo_files_wl_liver)], "../../treatment_DEG_list/weanling_DEG/", "lv_wl_DEG_Mdom_M");

expo_count_aged_lv_Fdom_F<-count_expo(all_CPM_aged_Fdom,expo_files_aged_liver[grep("_Female_", expo_files_aged_liver)], "../../treatment_DEG_list/aged_DEG/", "lv_aged_DEG_Fdom_F");
expo_count_aged_lv_Fdom_M<-count_expo(all_CPM_aged_Fdom,expo_files_aged_liver[grep("_Male_", expo_files_aged_liver)], "../../treatment_DEG_list/aged_DEG/", "lv_aged_DEG_Fdom_M");

expo_count_aged_lv_Mdom_F<-count_expo(all_CPM_aged_Mdom, expo_files_aged_liver[grep("_Female_", expo_files_aged_liver)], "../../treatment_DEG_list/aged_DEG/", "lv_aged_DEG_Mdom_F");
expo_count_aged_lv_Mdom_M<-count_expo(all_CPM_aged_Mdom, expo_files_aged_liver[grep("_Male_", expo_files_aged_liver)], "../../treatment_DEG_list/aged_DEG/", "lv_aged_DEG_Mdom_M");


##########
get_ctrl_plot<-function(data_for_plot,my_tag1, my_tag2){
data_for_plot$CPM<- as.numeric(as.character(data_for_plot$CPM));

mycol_1<- "red";
mycol_2<- "blue";
myshape_1<- 16;
myshape_2<- 15;
if(grepl("M", my_tag1)){
   mycol_1<- "blue";
   mycol_2<- "red";
   myshape_1<- 15;
   myshape_2<- 16;
  }

p<- ggplot(data_for_plot[data_for_plot$type==my_tag1,], aes(x=gene, y=log2(1+CPM), group=type)) + geom_line(col=mycol_1, size=0.5) + geom_point(shape=myshape_1, size=0.7, aes(col=as.character(data_for_plot[data_for_plot$type==paste0("expo_",strsplit(my_tag1, "_")[[1]][2]) ,"CPM"])));

p<- p+ theme_classic()
p<- p+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=5))
p<- p+ labs(y="log2CPM", x="Gene", title=paste0("RNAseq ", my_tag1," sexDEG n=", length(unique(data_for_plot$gene))), color="Overlap to expo");

p<- p+ geom_point(data=data_for_plot[data_for_plot$type==my_tag2,], aes(x=gene, y=log2(1+CPM), col=as.character(data_for_plot[data_for_plot$type==paste0("expo_",strsplit(my_tag2, "_")[[1]][2]) ,"CPM"])), size=0.7, shape=myshape_2);

hex<- hue_pal()(4)
p<- p+ scale_fill_manual(values=c("0"=hex[1], "1"=hex[2], "2"=hex[3], "3"=hex[4]));

out_name<- paste0("lineplot_sexDEG_", my_tag1, ".pdf");
pdf(out_name, height=6, width=9);
print(p);
dev.off();
}
##########

info_adt_Fdom<- cbind.data.frame(all_CPM_adt_Fdom, expo_count_adt_lv_Fdom_F, expo_count_adt_lv_Fdom_M)
colnames(info_adt_Fdom)[7:8]<- c("expo_F", "expo_M");
data_for_plot_adt_Fdom<- data.frame(info_adt_Fdom[,c(3:4,7:8)] %>% as_tibble(rownames="gene") %>% pivot_longer(cols=-gene, names_to="type", values_to="CPM") %>% group_by(type) )
data_for_plot_adt_Fdom$gene<- factor(data_for_plot_adt_Fdom$gene, levels=rownames(all_CPM_adt_Fdom));

info_adt_Mdom<- cbind.data.frame(all_CPM_adt_Mdom, expo_count_adt_lv_Mdom_F, expo_count_adt_lv_Mdom_M)
colnames(info_adt_Mdom)[7:8]<- c("expo_F", "expo_M");
data_for_plot_adt_Mdom<- data.frame(info_adt_Mdom[,c(3:4,7:8)] %>% as_tibble(rownames="gene") %>% pivot_longer(cols=-gene, names_to="type", values_to="CPM") %>% group_by(type) )
data_for_plot_adt_Mdom$gene<- factor(data_for_plot_adt_Mdom$gene, levels=rownames(all_CPM_adt_Mdom));

get_ctrl_plot(data_for_plot_adt_Fdom, "adt_F", "adt_M");
get_ctrl_plot(data_for_plot_adt_Mdom, "adt_M", "adt_F");
######
info_wl_Fdom<- cbind.data.frame(all_CPM_wl_Fdom, expo_count_wl_lv_Fdom_F, expo_count_wl_lv_Fdom_M)
colnames(info_wl_Fdom)[7:8]<- c("expo_F", "expo_M");
data_for_plot_wl_Fdom<- data.frame(info_wl_Fdom[,c(1:2,7:8)] %>% as_tibble(rownames="gene") %>% pivot_longer(cols=-gene, names_to="type", values_to="CPM") %>% group_by(type) )
data_for_plot_wl_Fdom$gene<- factor(data_for_plot_wl_Fdom$gene, levels=rownames(all_CPM_wl_Fdom));

info_wl_Mdom<- cbind.data.frame(all_CPM_wl_Mdom, expo_count_wl_lv_Mdom_F, expo_count_wl_lv_Mdom_M)
colnames(info_wl_Mdom)[7:8]<- c("expo_F", "expo_M");
data_for_plot_wl_Mdom<- data.frame(info_wl_Mdom[,c(1:2,7:8)] %>% as_tibble(rownames="gene") %>% pivot_longer(cols=-gene, names_to="type", values_to="CPM") %>% group_by(type) )
data_for_plot_wl_Mdom$gene<- factor(data_for_plot_wl_Mdom$gene, levels=rownames(all_CPM_wl_Mdom));

get_ctrl_plot(data_for_plot_wl_Fdom, "wl_F", "wl_M");
get_ctrl_plot(data_for_plot_wl_Mdom, "wl_M", "wl_F");
######
info_aged_Fdom<- cbind.data.frame(all_CPM_aged_Fdom, expo_count_aged_lv_Fdom_F, expo_count_aged_lv_Fdom_M)
colnames(info_aged_Fdom)[7:8]<- c("expo_F", "expo_M");
data_for_plot_aged_Fdom<- data.frame(info_aged_Fdom[,c(5:6,7:8)] %>% as_tibble(rownames="gene") %>% pivot_longer(cols=-gene, names_to="type", values_to="CPM") %>% group_by(type) )
data_for_plot_aged_Fdom$gene<- factor(data_for_plot_aged_Fdom$gene, levels=rownames(all_CPM_aged_Fdom));

info_aged_Mdom<- cbind.data.frame(all_CPM_aged_Mdom, expo_count_aged_lv_Mdom_F, expo_count_aged_lv_Mdom_M)
colnames(info_aged_Mdom)[7:8]<- c("expo_F", "expo_M");
data_for_plot_aged_Mdom<- data.frame(info_aged_Mdom[,c(5:6,7:8)] %>% as_tibble(rownames="gene") %>% pivot_longer(cols=-gene, names_to="type", values_to="CPM") %>% group_by(type) )
data_for_plot_aged_Mdom$gene<- factor(data_for_plot_aged_Mdom$gene, levels=rownames(all_CPM_aged_Mdom));

get_ctrl_plot(data_for_plot_aged_Fdom, "aged_F", "aged_M");
get_ctrl_plot(data_for_plot_aged_Mdom, "aged_M", "aged_F");
###
















#######
convert_lab<- function(my_lab){
if(my_lab=="Zhibin"){
   my_lab<- "ZB";
  }
if(my_lab=="Bartolomei"){
   my_lab<- "BA";
  }
if(my_lab=="Dolinoy"){
   my_lab<- "DO";
  }
if(my_lab=="Biswal"){
   my_lab<- "BI";
  }
if(my_lab=="Mutlu"){
   my_lab<- "MU";
  }
if(my_lab=="Walker"){
   my_lab<- "WK";
  }
if(my_lab=="Aylor"){
   my_lab<- "AL";
  }
return(my_lab);
}
#######
add_expo<- function(my_expo_list, my_data, my_tag1, my_tag2, my_stage, my_CPM){

expo_meta<- data.frame(t(data.frame(strsplit(my_expo_list, "_"))))
colnames(expo_meta)<- c("expo", "sex1", "stage1", "ctrl", "sex2",  "stage2", "lab", "tissue", "direction");

expo_tags<- unique(paste0(expo_meta[,"expo"],"_", expo_meta[,"lab"]))

for(i in 1:length(expo_tags)){

    expo<- strsplit(expo_tags[i], "_")[[1]][1]
    lab<- strsplit(expo_tags[i], "_")[[1]][2]

    sub_expo_list<- my_expo_list[which(expo_meta$expo==expo & expo_meta$lab==lab)]
    sub_expo_list_F<- sub_expo_list[grep("Female", sub_expo_list)]
    sub_expo_list_M<- sub_expo_list[grep("Male", sub_expo_list)]
    lab<- convert_lab(lab);

    expo_F<- rep(NA, nrow(my_data));
    expo_M<- rep(NA, nrow(my_data));

    deg_F<- c();
    deg_M<- c();
    for(j in 1:length(sub_expo_list_F)){
        data_F<- read.table(paste0("../../treatment_DEG_list/", my_stage, "_DEG/", sub_expo_list_F[j]), header=T,  sep="\t");
        deg_F<- c(deg_F, data_F$gene);
       }
    deg_F_sexDEG<- deg_F[deg_F %in% rownames(my_data)]
    expo_F[match(deg_F_sexDEG, rownames(my_data))]<- my_CPM[deg_F_sexDEG, intersect(grep(paste0(lab, "_", expo), colnames(my_CPM)) , grep("_F$", colnames(my_CPM)))]

    for(j in 1:length(sub_expo_list_M)){
        data_M<- read.table(paste0("../../treatment_DEG_list/", my_stage, "_DEG/", sub_expo_list_M[j]), header=T,  sep="\t");
        deg_M<- c(deg_M, data_M$gene);
       }
    deg_M_sexDEG<- deg_M[deg_M %in% rownames(my_data)]
    expo_M[match(deg_M_sexDEG, rownames(my_data))]<- my_CPM[deg_M_sexDEG, intersect(grep(paste0(lab, "_", expo), colnames(my_CPM)) , grep("_M$", colnames(my_CPM)))]

    my_data1<- cbind.data.frame(my_data, expo_F, expo_M)
### Only expo DEGs overlapped to sex DEGs
    my_data2<- my_data1[!(is.na(expo_F) & is.na(expo_M)),]

    if(my_stage=="weanling"){
       data_for_plot<- cbind.data.frame(rep(rownames(my_data2),2), c(my_data2$wl_F, my_data2$wl_M), c(rep(c("Fdom", "Mdom"), each=nrow(my_data2))));
      }
    if(my_stage=="adult"){
       data_for_plot<- cbind.data.frame(rep(rownames(my_data2),2), c(my_data2$adt_F, my_data2$adt_M), c(rep(c("Fdom", "Mdom"), each=nrow(my_data2))));
      }
    if(my_stage=="aged"){
       data_for_plot<- cbind.data.frame(rep(rownames(my_data2),2), c(my_data2$aged_F, my_data2$aged_M), c(rep(c("Fdom", "Mdom"), each=nrow(my_data2))));
      }

    colnames(data_for_plot)<- c("gene", "CPM", "type");
    data_for_plot$gene<- factor(data_for_plot$gene, levels=rownames(my_data2));
    data_for_plot$CPM<- as.numeric(as.character(data_for_plot$CPM));

    col_sex1<- "brown";
    col_sex2<- "blue";
    Col_sex1<- "brown";
    Col_sex2<- "blue";
    if(my_tag1=="Mdom"){
       col_sex1<- "blue";
       col_sex2<- "brown";
       Col_sex1<- "blue";
       Col_sex2<- "brown";
      }
    p<- ggplot(data_for_plot[data_for_plot$type==my_tag1,], aes(x=gene, y=log2(1+CPM), group=type)) + geom_line(col=col_sex1, size=1.2) + geom_point(col=col_sex1, size=2, shape=1);

    p<- p+ theme_classic()
    p<- p+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=5))
    p<- p+ labs(y="log2CPM", x="Gene", title=paste0(my_stage, " ", my_tag1," in ", expo_tags[i]," sexDEG n=", nrow(my_data2)));

    p<- p+ geom_point(data=data_for_plot[data_for_plot$type==my_tag2,], aes(x=gene, y=log2(1+CPM)), col=col_sex2, size=1.4, shape=1);

    segment_data<- data.frame(matrix(ncol=4, nrow=2*nrow(my_data2)));
    colnames(segment_data)<- c("gene", "type", "ystart", "yend");
    segment_data$gene<- rep(rownames(my_data2), 2);
    segment_data$type<- rep(c("Fdom", "Mdom"), each=nrow(my_data2));
    if(my_stage== "weanling"){
       segment_data$ystart<- c(my_data2$wl_F, my_data2$wl_M);
      }
    if(my_stage== "adult"){
       segment_data$ystart<- c(my_data2$adt_F, my_data2$adt_M);
      }
    if(my_stage== "aged"){
       segment_data$ystart<- c(my_data2$aged_F, my_data2$aged_M);
      }
    segment_data$yend<- c(my_data2$expo_F, my_data2$expo_M);
    segment_data<- segment_data[!is.na(segment_data$yend),]
    segment_data$ystart<- as.numeric(as.character(segment_data$ystart));
    segment_data$yend<- as.numeric(as.character(segment_data$yend));
    q<- p+ geom_segment(data=segment_data, aes(x=gene, xend=gene, y=log2(ystart+1), yend=log2(yend+1), col=type),size=1.2) 
    q<- q+ scale_color_manual(values=c("Fdom"="brown", "Mdom"="blue"));
    q<- q+ geom_point(data=segment_data[segment_data$type==my_tag1 & segment_data$ystart -segment_data$yend<0,], mapping=aes(x=gene, y=log2(yend+1)), size=2.2, fill=Col_sex1, col=Col_sex1, shape=24);
    q<- q+ geom_point(data=segment_data[segment_data$type==my_tag1 & segment_data$ystart -segment_data$yend>0,], mapping=aes(x=gene, y=log2(yend+1)), size=2.2, fill=Col_sex1, col=Col_sex1,  shape=25);

    q<- q+ geom_point(data=segment_data[segment_data$type==my_tag2 & segment_data$ystart -segment_data$yend<0,], mapping=aes(x=gene, y=log2(yend+1)), size=2.2, fill=Col_sex2, col=Col_sex2, shape=24);
    q<- q+ geom_point(data=segment_data[segment_data$type==my_tag2 & segment_data$ystart -segment_data$yend>0,], mapping=aes(x=gene, y=log2(yend+1)), size=2.2, fill=Col_sex2, col=Col_sex2, shape=25);

    q<- q+ ggtitle(paste0(my_stage, " ", my_tag1, " ", expo, " ", lab, " sexDEG in expoDEG ", substr(my_tag1,1,1), "/", substr(my_tag2,1,1),"=", nrow(segment_data[segment_data$type==my_tag1,]), " ", nrow(segment_data[segment_data$type==my_tag2,])));


    out_name<- paste0("dotline_lv_", my_stage, "_", my_tag1, "_", expo, "_", lab, ".pdf");
    pdf(out_name, width=9, height=6);
    print(q);
    dev.off();
   }
}
#######

add_expo(expo_files_adt_liver, all_CPM_adt_Fdom, "Fdom", "Mdom", "adult", CPM_lv_adt);
add_expo(expo_files_adt_liver, all_CPM_adt_Mdom, "Mdom", "Fdom", "adult", CPM_lv_adt);

add_expo(expo_files_wl_liver, all_CPM_wl_Fdom, "Fdom", "Mdom", "weanling", CPM_lv_wl);
add_expo(expo_files_wl_liver, all_CPM_wl_Mdom, "Mdom", "Fdom", "weanling", CPM_lv_wl);

add_expo(expo_files_aged_liver, all_CPM_aged_Fdom, "Fdom", "Mdom", "aged", CPM_lv_aged);
add_expo(expo_files_aged_liver, all_CPM_aged_Mdom, "Mdom", "Fdom", "aged", CPM_lv_aged);







