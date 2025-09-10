# 2024-01-11 GO for F-/M-specific, shared
library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)
library(pheatmap)

#https://www.littlemissdata.com/blog/set-analysis

path_wl<- "/BRC/shuhua/target/RNAseq/DEGs/treatment_DEG_list/GO_R_combined_direction/weanling_GO_R/"
path_adt<- "/BRC/shuhua/target/RNAseq/DEGs/treatment_DEG_list/GO_R_combined_direction/adult_GO_R/"
path_aged<- "/BRC/shuhua/target/RNAseq/DEGs/treatment_DEG_list/GO_R_combined_direction/aged_GO_R/"

list_wl<- list.files(path_wl)
list_wl<- list_wl[grep("liver", list_wl)]
list_adt<- list.files(path_adt)
list_adt<- list_adt[grep("liver", list_adt)]
list_aged<- list.files(path_aged)
list_aged<- list_aged[grep("liver", list_aged)]

data_info<- unique(data.frame(t(data.frame(strsplit(list_wl, "_"))))[,c(2,8)])

GO_shared_data<- data.frame(matrix(ncol=6, nrow=0));
colnames(GO_shared_data)<- c("GO", "F_qval", "M_qval", "expo", "lab", "stage");

#####
get_dotplot<- function(my_data, my_sex){

colnames(my_data)<- c("wl", "adt", "aged");
data_for_plot<- data.frame(my_data %>% tibble::rownames_to_column(var="tag") %>% pivot_longer(-tag))
data_for_plot$name<- factor(data_for_plot$name, levels=c("wl", "adt", "aged"));
data_for_plot$tag<- factor(data_for_plot$tag, levels=rev(c("As", "BPA10ug", "BPA10mg", "DEHP", "Pb", "PM2.5_Biswal", "PM2.5_Mutlu", "TBT", "TCDD")))

p<- ggplot(data_for_plot, aes(x=name, y=tag, color=value, size=value)) + geom_point()
p<- p+ scale_color_gradient(low = "lightblue", high = "darkblue", limits=c(1, 300))
p<- p+ scale_size(limits=c(1, 300))
p<- p+ theme_bw()
p<- p+ ggtitle(paste0("GO ", my_sex));
return(p);
}
######
get_shared_GO_dotplot<- function(my_data, my_stage, my_expo, my_lab){

data_for_plot1<- data.frame(my_data[,1:2] %>% tibble::rownames_to_column(var="tag") %>% pivot_longer(-tag, values_to="qval"))
data_for_plot2<- data.frame(my_data[,3:4] %>% tibble::rownames_to_column(var="tag") %>% pivot_longer(-tag, values_to="count"))
data_for_plot<- cbind.data.frame(data_for_plot1, data_for_plot2$count);
colnames(data_for_plot)[4]<- "ratio";
data_for_plot$name<- ifelse(grepl("F_", data_for_plot$name), "F", "M")
data_for_plot$ratio<- as.numeric(data_for_plot$ratio)
data_for_plot$qval<- as.numeric(data_for_plot$qval)
data_for_plot$tag<- factor(data_for_plot$tag, levels=rev(unique(data_for_plot1[order(data_for_plot1$qval),"tag"])))

p<- ggplot(data_for_plot, aes(x=name, y=tag, color=-log10(qval), size=ratio)) + geom_point()
p<- p+ scale_color_gradient(low = "lightblue", high = "darkblue") 
p<- p+ theme_bw()
p<- p+ ggtitle(paste0("GO ", my_stage, " ", my_expo, " ", my_lab, " n=", nrow(my_data)));

pdf(paste0("dotplot_shared_GO_", my_stage, "_", my_expo, "_", my_lab, ".pdf"));
print(p);
dev.off();
}
######
get_venn_plot<- function(my_expo, my_lab){

sub_list_wl<- list_wl[grepl(my_lab, list_wl) & grepl(my_expo, list_wl) ]
sub_list_adt<- list_adt[grepl(my_lab, list_adt) & grepl(my_expo, list_adt) ]
sub_list_aged<- list_aged[grepl(my_lab, list_aged) & grepl(my_expo, list_aged) ]

gene_cutoff<- 5;
qval_cutoff<- 0.05
data_wl_F<- read.table(paste0(path_wl, sub_list_wl[grepl("Female", sub_list_wl)]), header=T, sep="\t", quote="")
data_wl_F<- data_wl_F[as.numeric(data_wl_F$Count)>= gene_cutoff & as.numeric(data_wl_F$qvalue)< qval_cutoff,]
data_wl_M<- read.table(paste0(path_wl, sub_list_wl[grepl("Male", sub_list_wl)]), header=T, sep="\t", quote="")
data_wl_M<- data_wl_M[as.numeric(data_wl_M$Count)>= gene_cutoff & as.numeric(data_wl_M$qvalue)< qval_cutoff,]

data_adt_F<- read.table(paste0(path_adt, sub_list_adt[grepl("Female", sub_list_adt)]), header=T, sep="\t", quote="")
data_adt_F<- data_adt_F[as.numeric(data_adt_F$Count)>= gene_cutoff & as.numeric(data_adt_F$qvalue)< qval_cutoff,]
data_adt_M<- read.table(paste0(path_adt, sub_list_adt[grepl("Male", sub_list_adt)]), header=T, sep="\t", quote="")
data_adt_M<- data_adt_M[as.numeric(data_adt_M$Count)>= gene_cutoff & as.numeric(data_adt_M$qvalue)< qval_cutoff,]

count_wl_shared<- length(data_wl_F$ID[data_wl_F$ID %in% data_wl_M$ID])
count_wl_F<- length(data_wl_F$ID[!(data_wl_F$ID %in% data_wl_M$ID)])
count_wl_M<- length(data_wl_M$ID[!(data_wl_M$ID %in% data_wl_F$ID)])

count_adt_shared<- length(data_adt_F$ID[data_adt_F$ID %in% data_adt_M$ID])
count_adt_F<- length(data_adt_F$ID[!(data_adt_F$ID %in% data_adt_M$ID)])
count_adt_M<- length(data_adt_M$ID[!(data_adt_M$ID %in% data_adt_F$ID)])

if(count_wl_shared>0){
   data_shared_wl<- data.frame(matrix(0, ncol=4, nrow=count_wl_shared));
   colnames(data_shared_wl)<- c("F_qval", "M_qval", "F_ratio", "M_ratio");
   rownames(data_shared_wl)<- data_wl_F[data_wl_F$ID %in% data_wl_M$ID,"Description"]
   data_shared_wl$F_qval<- data_wl_F[data_wl_F$ID %in% data_wl_M$ID,"qvalue"]
   data_shared_wl$F_ratio<- data_wl_F[data_wl_F$ID %in% data_wl_M$ID,"Count"]
   data_shared_wl$M_qval<- data_wl_M[match(rownames(data_shared_wl), data_wl_M$Description),"qvalue"]
   data_shared_wl$M_ratio<- data_wl_M[match(rownames(data_shared_wl), data_wl_M$Description),"Count"]
   gene_size<- as.numeric(strsplit(data_wl_F[1,"GeneRatio"], "/")[[1]][2])
   data_shared_wl$F_ratio<- as.numeric(data_shared_wl$F_ratio)/gene_size
   data_shared_wl$M_ratio<- as.numeric(data_shared_wl$M_ratio)/gene_size
   data_GO_expo<- cbind.data.frame(row.names(data_shared_wl), data_shared_wl$F_qval, data_shared_wl$M_qval, my_expo, my_lab, "wl")
   colnames(data_GO_expo)<- c("GO", "F_qval", "M_qval", "expo", "lab", "stage");
   GO_shared_data<<- rbind.data.frame(GO_shared_data, data_GO_expo)
   get_shared_GO_dotplot(data_shared_wl, "wl", my_expo, my_lab);
  }
if(count_adt_shared>0){
   data_shared_adt<- data.frame(matrix(0, ncol=4, nrow=count_adt_shared));
   colnames(data_shared_adt)<- c("F_qval", "M_qval", "F_ratio", "M_ratio");
   rownames(data_shared_adt)<- data_adt_F[data_adt_F$ID %in% data_adt_M$ID,"Description"]
   data_shared_adt$F_qval<- data_adt_F[data_adt_F$ID %in% data_adt_M$ID,"qvalue"]
   data_shared_adt$F_ratio<- data_adt_F[data_adt_F$ID %in% data_adt_M$ID,"Count"]
   data_shared_adt$M_qval<- data_adt_M[match(rownames(data_shared_adt), data_adt_M$Description),"qvalue"]
   data_shared_adt$M_ratio<- data_adt_M[match(rownames(data_shared_adt), data_adt_M$Description),"Count"]
   gene_size<- as.numeric(strsplit(data_adt_F[1,"GeneRatio"], "/")[[1]][2])
   data_shared_adt$F_ratio<- as.numeric(data_shared_adt$F_ratio)/gene_size
   data_shared_adt$M_ratio<- as.numeric(data_shared_adt$M_ratio)/gene_size
   data_GO_expo<- cbind.data.frame(row.names(data_shared_adt), data_shared_adt$F_qval, data_shared_adt$M_qval, my_expo, my_lab, "adt")
   colnames(data_GO_expo)<- c("GO", "F_qval", "M_qval", "expo", "lab", "stage");
   GO_shared_data<<- rbind.data.frame(GO_shared_data, data_GO_expo)
   get_shared_GO_dotplot(data_shared_adt, "adt", my_expo, my_lab);
  }

count_aged_shared<- 0;
count_aged_F<- 0;
count_aged_M<- 0;

if(length(sub_list_aged)==0){
  temp<- "tmp"
  }else{
   data_aged_F<- read.table(paste0(path_aged, sub_list_aged[grepl("Female", sub_list_aged)]), header=T, sep="\t", quote="")
   data_aged_F<- data_aged_F[as.numeric(data_aged_F$Count)>= gene_cutoff & as.numeric(data_aged_F$qvalue)< qval_cutoff,]
   data_aged_M<- read.table(paste0(path_aged, sub_list_aged[grepl("Male", sub_list_aged)]), header=T, sep="\t", quote="")
   data_aged_M<- data_aged_M[as.numeric(data_aged_M$Count)>= gene_cutoff & as.numeric(data_aged_M$qvalue)< qval_cutoff,]

   count_aged_shared<- length(data_aged_F$ID[data_aged_F$ID %in% data_aged_M$ID])
   count_aged_F<- length(data_aged_F$ID[!(data_aged_F$ID %in% data_aged_M$ID)])
   count_aged_M<- length(data_aged_M$ID[!(data_aged_M$ID %in% data_aged_F$ID)])
  }
if(count_aged_shared>0){
   data_shared_aged<- data.frame(matrix(0, ncol=4, nrow=count_aged_shared));
   colnames(data_shared_aged)<- c("F_qval", "M_qval", "F_ratio", "M_ratio");
   rownames(data_shared_aged)<- data_aged_F[data_aged_F$ID %in% data_aged_M$ID,"Description"]
   data_shared_aged$F_qval<- data_aged_F[data_aged_F$ID %in% data_aged_M$ID,"qvalue"]
   data_shared_aged$F_ratio<- data_aged_F[data_aged_F$ID %in% data_aged_M$ID,"Count"]
   data_shared_aged$M_qval<- data_aged_M[match(rownames(data_shared_aged), data_aged_M$Description),"qvalue"]
   data_shared_aged$M_ratio<- data_aged_M[match(rownames(data_shared_aged), data_aged_M$Description),"Count"]
   gene_size<- as.numeric(strsplit(data_aged_F[1,"GeneRatio"], "/")[[1]][2])
   data_shared_aged$F_ratio<- as.numeric(data_shared_aged$F_ratio)/gene_size
   data_shared_aged$M_ratio<- as.numeric(data_shared_aged$M_ratio)/gene_size
   data_GO_expo<- cbind.data.frame(row.names(data_shared_aged), data_shared_aged$F_qval, data_shared_aged$M_qval, my_expo, my_lab, "aged")
   colnames(data_GO_expo)<- c("GO", "F_qval", "M_qval", "expo", "lab", "stage");
   GO_shared_data<<- rbind.data.frame(GO_shared_data, data_GO_expo)
   get_shared_GO_dotplot(data_shared_aged, "aged", my_expo, my_lab);
  }

return(c(count_wl_F, count_wl_shared, count_wl_M, count_adt_F, count_adt_shared, count_adt_M, count_aged_F, count_aged_shared, count_aged_M));
}
######

data_to_plot<- data.frame(matrix(0, nrow=nrow(data_info), ncol=9));
data_info<- data_info[c(1,3,2,4:9),]
rownames(data_to_plot)<- ifelse(data_info[,1]=="PM2.5", paste0(data_info[,1], "_", data_info[,2]), data_info[,1])
colnames(data_to_plot)<- c("wl_F", "wl_shared", "wl_M",  "adt_F", "adt_shared", "adt_M",  "aged_F", "aged_shared", "aged_M")

for(i in 1:nrow(data_info)){
    data_to_plot[i,]<- get_venn_plot(data_info[i,1], data_info[i, 2]) 
   }

####

p<- get_dotplot(data_to_plot[,grep("_F", colnames(data_to_plot))], "F_specific");
q<- get_dotplot(data_to_plot[,grep("_M", colnames(data_to_plot))], "M_specific");
s<- get_dotplot(data_to_plot[,grep("shared", colnames(data_to_plot))], "shared");

pdf("GO_count_across_stages.pdf", height=3, width=8);
grid.arrange(p, q, s, ncol=3);
dev.off();

######
get_heatmap<- function(my_stage){

GO_shared_stage<- GO_shared_data[GO_shared_data$stage==my_stage,];
GO_shared_stage$tag<- paste0(GO_shared_stage$expo, "_", GO_shared_stage$lab);
GO_vec<- unique(GO_shared_stage$GO);
tag_vec<- unique(GO_shared_stage$tag);

data_heatmap<- data.frame(matrix(1, nrow=length(GO_vec), ncol= 2*length(tag_vec)));
rownames(data_heatmap)<- GO_vec;
colnames(data_heatmap)<- paste0(rep(tag_vec, each=2), c("_F", "_M"));

for(i in 1:ncol(data_heatmap)){
    tag<- gsub("_F$|_M$", "", colnames(data_heatmap)[i])
    GO_shared_expo<- GO_shared_stage[GO_shared_stage$tag==tag,]
    if(grepl("_F", colnames(data_heatmap)[i])){
       data_heatmap[,i]<- as.numeric(GO_shared_expo[match(rownames(data_heatmap), GO_shared_expo$GO),"F_qval"]);
      }else{
       data_heatmap[,i]<- as.numeric(GO_shared_expo[match(rownames(data_heatmap), GO_shared_expo$GO),"M_qval"]);
      }
   }
data_heatmap[is.na(data_heatmap)]<- 1 
data_heatmap_log<- (-1)*log10(data_heatmap)
colnames(data_heatmap_log)<- gsub("_Walker|_Zhibin|_Bartolomei|_Dolinoy|_Aylor", "", colnames(data_heatmap_log))

paletteLength <- 31
myColor <- colorRampPalette(c("navy", "white", "firebrick3"))(paletteLength)
myBreaks <- c(seq(min(data_heatmap_log), 0, length.out=ceiling(paletteLength/2) + 1),
          seq(max(data_heatmap_log)/paletteLength, max(data_heatmap_log), length.out=floor(paletteLength/2)))

exposure<- gsub("_F$|_M$", "", colnames(data_heatmap_log))
sex<- ifelse(grepl("_F", colnames(data_heatmap)), "F", "M")
ann<- cbind.data.frame(exposure, sex);
rownames(ann)<-  colnames(data_heatmap_log)
ann$exposure<- factor(ann$exposure, levels=c("As", "BPA10ug", "BPA10mg", "DEHP", "Pb", "PM2.5_Biswal", "PM2.5_Mutlu", "TBT", "TCDD"));

pdf(paste0("heatmap_shared_expoGO_LFC_", my_stage, ".pdf"), width=12, height=10);
print(pheatmap(data_heatmap_log, cluster_cols=F, cluster_rows=T,  annotation_col=ann, main=paste0("-log10Q shared expoGO  n=", nrow(data_heatmap_log)), gaps_col=seq(2, ncol(data_heatmap_log),2) , color = colorRampPalette(c( "#FFFFFF", "#3d85c6" ))(30)   ))
dev.off();

write.table(data.frame("qval"=rownames(data_heatmap_log), data_heatmap_log), file=paste0("data_shared_expoGO_LFC_", my_stage, ".txt"), quote=F, sep="\t", row.names=F, col.names=T);

}
######

get_heatmap("wl");
get_heatmap("adt");
get_heatmap("aged");

