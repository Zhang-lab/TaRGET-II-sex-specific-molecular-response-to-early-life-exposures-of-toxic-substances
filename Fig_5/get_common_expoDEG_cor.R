#2024-01-25 LFC for non-expoDEG is not 0
library(tidyr)
library(ggplot2)
library(pheatmap)
library(corrplot)

cor.mtest <- function(mat, ...) {
    mat <- as.matrix(mat)
    n <- ncol(mat)
    p.mat<- matrix(NA, n, n)
    diag(p.mat) <- 0
    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            tmp <- cor.test(mat[, i], mat[, j], ...)
            p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
        }
    }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

file_path<- "../treatment_DEG_list_noChrXY/"

#######
get_data<- function(my_stage, my_sex){

file_path_stage<- paste0(file_path, my_stage, "_DEG/")
file_list<- list.files(file_path_stage)

file_list<- file_list[intersect(grep("liver", file_list), grep(my_sex, file_list))]

data_info<- data.frame(t(data.frame(strsplit(file_list, "_"))))
expo_vec<- unique(ifelse(data_info[,7] %in% c("Biswal", "Mutlu"), paste0(data_info[,1], "_", data_info[,7]), data_info[,1]))
expo_vec<- factor(expo_vec, levels=c("As", "BPA10ug", "BPA10mg", "DEHP", "Pb", "PM2.5_Biswal", "PM2.5_Mutlu", "TBT", "TCDD"));
expo_vec<- expo_vec[order(expo_vec)]

list_out<- list();
list_count<- 1;

for(i in 1:length(file_list)){
    expo<- strsplit(file_list[i], "_")[[1]][1]
    lab<- strsplit(file_list[i], "_")[[1]][7]
    data<- read.table(paste0(file_path_stage, file_list[i]), header=T, sep="\t")
    for(j in 1:nrow(data)){
        list_out[[list_count]]<- c(data[j, "gene"], data[j, "log2FoldChange"], expo, lab); 
        list_count<- list_count +1;
       }
   }
data_out<- data.frame(do.call(rbind, list_out))
colnames(data_out)<- c("gene", "LFC", "expo", "lab");
data_out$expo<- ifelse(data_out$expo=="PM2.5", paste0(data_out$expo, "_", data_out$lab), data_out$expo)

gene_freq<- data.frame(table(data_out$gene))
gene_freq<- gene_freq[order(-gene_freq$Freq),]
data_freq<- data.frame(matrix(nrow=length(expo_vec), ncol= max(gene_freq$Freq)));
colnames(data_freq)<- 1:max(gene_freq$Freq)
rownames(data_freq)<- expo_vec;

for(i in 1:nrow(data_freq)){
    freq_vec<- gene_freq[gene_freq$Var1 %in% data_out[data_out$expo== expo_vec[i],"gene"], "Freq"];
    for(j in 1:ncol(data_freq)){
        data_freq[i, j]<-length(which(freq_vec==as.numeric(colnames(data_freq)[j])))
       }
   }
data_for_plot<- data.frame(data_freq %>% tibble::rownames_to_column(var="expo") %>% pivot_longer(-expo, names_to="n_expo", values_to="count") )
data_freq$gene_sum<- rowSums(data_freq)

data_for_plot$expo<- factor(data_for_plot$expo, levels=rev(c("As", "BPA10ug", "BPA10mg", "DEHP", "Pb", "PM2.5_Biswal", "PM2.5_Mutlu", "TBT", "TCDD")))
data_for_plot$n_expo<- ifelse(data_for_plot$n_expo=="1", "expo_specific", data_for_plot$n_expo)
data_for_plot$count<- as.numeric(data_for_plot$count)
data_for_plot$n_expo<- factor(data_for_plot$n_expo, levels=c(max(gene_freq$Freq):2, "expo_specific"));

p<- ggplot(data_for_plot, aes(x=expo, y=count, fill=n_expo))+geom_bar(position="fill", stat="identity")
p<- p + scale_fill_brewer(palette = "Set3")
p<- p+ theme_classic()
p<- p+ coord_flip()
p<- p+ scale_y_continuous(expand = c(0, 0))
p<- p+ ggtitle(paste0(my_stage, " ", my_sex));
p<- p+ labs(fill="Common by", y="Ratio", x="Exposure");
p<- p+ geom_text(data=data_freq, aes(rownames(data_freq), 0.2, label=gene_sum, fill=NULL))

pdf(paste0("across_expo_plots/stacked_barplot_", my_stage, "_", my_sex, ".pdf"), height=5);
print(p);
dev.off();

freq_cutoff<- max(gene_freq$Freq) -2;
if(freq_cutoff<= 2){
   freq_cutoff<- 3;
  }
genes_selected<- as.character(gene_freq[gene_freq$Freq >= freq_cutoff, "Var1"])

data_LFC<- data.frame(matrix(0, nrow=nrow(gene_freq), ncol=length(expo_vec)));
rownames(data_LFC)<- gene_freq[, "Var1"]; 
colnames(data_LFC)<- expo_vec;
for(i in 1:nrow(data_out)){
    data_LFC[data_out[i,"gene"],data_out[i, "expo"]]<-as.numeric(data_out[i, "LFC"])
   }
genes_selected<- rownames(data_LFC[rowSums(data_LFC[,c("TBT", "TCDD")]==0)<2,])

data_LFC_expo_ge2<- data_LFC[rowSums(data_LFC!=0)>=2,]
data_heatmap<- data_LFC[genes_selected,]

### Check LFC to replace 0 for non-expoDEGs
path_genes<- "../treatment_all_list/"
list_genes<- list.files(paste0(path_genes, my_stage, "_all/"))
list_genes<- list_genes[intersect(grep("liver", list_genes), grep(my_stage, list_genes))]
list_genes<- list_genes[grep(my_sex, list_genes)]
info_all<- data.frame(t(data.frame(strsplit(list_genes, "_"))))
info_all$tag<- ifelse(info_all$X1== "PM2.5", paste0(info_all$X1, "_", info_all$X7),info_all[,1])
for(i in 1:ncol(data_LFC_expo_ge2)){
    file_n<- which(info_all$tag==colnames(data_LFC_expo_ge2)[i])
    for(j in 1:length(file_n)){
        data<- read.table(paste0(path_genes, my_stage, "_all/", list_genes[file_n[j]]), header=T, sep="\t");
        for(k in 1:nrow(data_LFC_expo_ge2)){
            if(data_LFC_expo_ge2[k, i]==0 &  length(which(data$gene==rownames(data_LFC_expo_ge2)[k]) )>0 ){
               data_LFC_expo_ge2[k, i]<- data[data$gene==rownames(data_LFC_expo_ge2)[k], "log2FoldChange"];
              }
           }
       }
   }


paletteLength <- 31
myColor <- colorRampPalette(c("navy", "white", "firebrick3"))(paletteLength)
myBreaks <- c(seq(min(data_heatmap), 0, length.out=ceiling(paletteLength/2) + 1),
          seq(max(data_heatmap)/paletteLength, max(data_heatmap), length.out=floor(paletteLength/2)))

exposure<- gsub("_F$|_M$", "", colnames(data_heatmap))

ann<- cbind.data.frame(exposure);
rownames(ann)<-  colnames(data_heatmap)
ann$exposure<- factor(ann$exposure, levels=c("As", "BPA10ug", "BPA10mg", "DEHP", "Pb", "PM2.5_Biswal", "PM2.5_Mutlu", "TBT", "TCDD"));

frequency<- data.frame(gene_freq[match(genes_selected, gene_freq$Var1),"Freq"])
rownames(frequency)<-  rownames(data_heatmap)
colnames(frequency)<- ("Frequency")

freq_data<- data.frame(table(frequency$Frequency))
for_gap_row<- cumsum(rev(freq_data$Freq))
for_gap_row<- for_gap_row[1:(length(for_gap_row)-1)]

pdf(paste0("across_expo_plots/heatmap_common_expoDEG_LFC_", my_stage, "_", my_sex, ".pdf"), width=4, height=10);
pheatmap(data_heatmap, cluster_cols=F, cluster_rows=F, color=myColor, annotation_col=ann, annotation_row=frequency, main=paste0("LFC of frequent DEGs n=", nrow(data_heatmap)), gaps_row= for_gap_row, breaks=myBreaks)
dev.off();

corr_position<- ifelse(my_sex=="Female", "upper", "lower");
M<- cor(data_LFC_expo_ge2)
p.mat<- cor.mtest(data_LFC_expo_ge2)
pdf(paste0("across_expo_plots/corr_common_expoGO_ge2_", my_stage, "_",  my_sex, ".pdf"))
col <- colorRampPalette(rev(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA")))
print(corrplot(M, method="color", col=col(200),
         type=corr_position, order="original",
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = p.mat, sig.level = 0.01, insig="pch",
         # hide correlation coefficient on the principal diagonal
         diag=FALSE
         ));
dev.off();

### scatter plot for TBT vs TCDD in male adt
if(my_stage=="adult" & my_sex=="Male"){
   data_scatter<- data_LFC[,c("TBT", "TCDD")]
   data_scatter<- data_scatter[rowSums(data_scatter==0)<2,]
  }
}
#######

get_data("weanling", "Female");
get_data("weanling", "Male");

get_data("adult", "Female");
get_data("adult", "Male");

get_data("aged", "Female");
get_data("aged", "Male");

