#2023-04-24 color different omics w/ different colors
library(ggplot2)
library(gridExtra)

data_RNA<- read.table("/BRC/shuhua/target/RNAseq/DEGs/sex_overlap_expo/sexDEG_overlap_expoDEG_adt", header=T)
data_ATAC<- read.table("/BRC/shuhua/target/ATAC/sex_overlap_expo/sexDAR_overlap_liver_adt_expoDAR", header=T)
data_WGBS<- read.table("/BRC/shuhua/target/sex_WGBS/sex_overlap_expo/expoDMR_overlap_liver_sexDMR", header=T)

expo_RNA<- sapply(strsplit(data_RNA[,1], "_"), head, 1)
expo_RNA[intersect(which(grepl("Biswal", data_RNA[,1])), which(grepl("PM2.5", data_RNA[,1])))]<- "PM2.5_Biswal";
expo_RNA[intersect(which(grepl("Mutlu", data_RNA[,1])), which(grepl("PM2.5", data_RNA[,1])))]<- "PM2.5_Mutlu";
data_RNA[,1]<- expo_RNA
expo_RNA<- factor(expo_RNA, levels=c("As", "BPA10ug", "BPA10mg", "DEHP", "Pb", "PM2.5_Biswal","PM2.5_Mutlu", "TBT", "TCDD"))
expo_RNA<- expo_RNA[order(expo_RNA)]

expo_ATAC<- as.character(data.frame(strsplit(data_ATAC[,1], "_"))[3,])
expo_ATAC[intersect(which(grepl("^BI_", data_ATAC[,1])), which(grepl("PM2.5", data_ATAC[,1])))]<- "PM2.5_Biswal";
expo_ATAC[intersect(which(grepl("^MU_", data_ATAC[,1])), which(grepl("PM2.5", data_ATAC[,1])))]<- "PM2.5_Mutlu";
data_ATAC[,1]<- expo_ATAC

expo_WGBS<- data_WGBS[,1]
expo_WGBS[expo_WGBS=="Lead"]<- "Pb"
expo_WGBS[expo_WGBS=="LOW_BPA"]<- "BPA10ug"
expo_WGBS[expo_WGBS=="UP_BPA"]<- "BPA10mg"
expo_WGBS[expo_WGBS=="PM25_BI"]<- "PM2.5_Biswal"
expo_WGBS[expo_WGBS=="PM25_UC"]<- "PM2.5_Mutlu"
data_WGBS[,1]<- expo_WGBS

data_for_plot_F_expoF<- data.frame(matrix(nrow=54, ncol=4));
data_for_plot_F_expoM<- data.frame(matrix(nrow=54, ncol=4));
data_for_plot_M_expoF<- data.frame(matrix(nrow=54, ncol=4));
data_for_plot_M_expoM<- data.frame(matrix(nrow=54, ncol=4));
colnames(data_for_plot_F_expoF)<- c("expo", "percent", "direction", "count");
colnames(data_for_plot_F_expoM)<- c("expo", "percent", "direction", "count");
colnames(data_for_plot_M_expoF)<- c("expo", "percent", "direction", "count");
colnames(data_for_plot_M_expoM)<- c("expo", "percent", "direction", "count");

for(i in 1:length(expo_RNA)){

    RNA_F<- as.numeric(strsplit(colnames(data_RNA)[1], "_")[[1]][2])
    RNA_M<- as.numeric(strsplit(colnames(data_RNA)[1], "_")[[1]][4])

    per_RNA_F_Fexpo_up<- round(100*data_RNA[data_RNA[,1]==expo_RNA[i],][2]/RNA_F, 2)
    per_RNA_F_Fexpo_dn<- round(100*data_RNA[data_RNA[,1]==expo_RNA[i],][3]/RNA_F, 2)

    per_RNA_F_Mexpo_up<- round(100*data_RNA[data_RNA[,1]==expo_RNA[i],][4]/RNA_F, 2)
    per_RNA_F_Mexpo_dn<- round(100*data_RNA[data_RNA[,1]==expo_RNA[i],][5]/RNA_F, 2)

    per_RNA_M_Fexpo_up<- round(100*data_RNA[data_RNA[,1]==expo_RNA[i],][6]/RNA_M, 2)
    per_RNA_M_Fexpo_dn<- round(100*data_RNA[data_RNA[,1]==expo_RNA[i],][7]/RNA_M, 2)

    per_RNA_M_Mexpo_up<- round(100*data_RNA[data_RNA[,1]==expo_RNA[i],][8]/RNA_M, 2)
    per_RNA_M_Mexpo_dn<- round(100*data_RNA[data_RNA[,1]==expo_RNA[i],][9]/RNA_M, 2)

    data_for_plot_F_expoF[6*i-5,]<- c(paste0(expo_RNA[i], "_RNA"), per_RNA_F_Fexpo_up, "up", as.numeric(data_RNA[data_RNA[,1]==expo_RNA[i],][2]));
    data_for_plot_F_expoF[6*i-4,]<- c(paste0(expo_RNA[i], "_RNA"), per_RNA_F_Fexpo_dn, "dn", as.numeric(data_RNA[data_RNA[,1]==expo_RNA[i],][3]));

    data_for_plot_F_expoM[6*i-5,]<- c(paste0(expo_RNA[i], "_RNA"), per_RNA_F_Mexpo_up, "up", as.numeric(data_RNA[data_RNA[,1]==expo_RNA[i],][4]));
    data_for_plot_F_expoM[6*i-4,]<- c(paste0(expo_RNA[i], "_RNA"), per_RNA_F_Mexpo_dn, "dn", as.numeric(data_RNA[data_RNA[,1]==expo_RNA[i],][5]));

    data_for_plot_M_expoF[6*i-5,]<- c(paste0(expo_RNA[i], "_RNA"), per_RNA_M_Fexpo_up, "up", as.numeric(data_RNA[data_RNA[,1]==expo_RNA[i],][6]));
    data_for_plot_M_expoF[6*i-4,]<- c(paste0(expo_RNA[i], "_RNA"), per_RNA_M_Fexpo_dn, "dn", as.numeric(data_RNA[data_RNA[,1]==expo_RNA[i],][7]));

    data_for_plot_M_expoM[6*i-5,]<- c(paste0(expo_RNA[i], "_RNA"), per_RNA_M_Mexpo_up, "up", as.numeric(data_RNA[data_RNA[,1]==expo_RNA[i],][8]));
    data_for_plot_M_expoM[6*i-4,]<- c(paste0(expo_RNA[i], "_RNA"), per_RNA_M_Mexpo_dn, "dn", as.numeric(data_RNA[data_RNA[,1]==expo_RNA[i],][9]));

#######
    ATAC_F<- as.numeric(strsplit(colnames(data_ATAC)[1], "_")[[1]][2])
    ATAC_M<- as.numeric(strsplit(colnames(data_ATAC)[1], "_")[[1]][4])

    per_ATAC_F_Fexpo_up<- round(100*data_ATAC[data_ATAC[,1]==expo_RNA[i],][2]/ATAC_F, 2)
    per_ATAC_F_Fexpo_dn<- round(100*data_ATAC[data_ATAC[,1]==expo_RNA[i],][3]/ATAC_F, 2)

    per_ATAC_F_Mexpo_up<- round(100*data_ATAC[data_ATAC[,1]==expo_RNA[i],][4]/ATAC_F, 2)
    per_ATAC_F_Mexpo_dn<- round(100*data_ATAC[data_ATAC[,1]==expo_RNA[i],][5]/ATAC_F, 2)

    per_ATAC_M_Fexpo_up<- round(100*data_ATAC[data_ATAC[,1]==expo_RNA[i],][6]/ATAC_M, 2)
    per_ATAC_M_Fexpo_dn<- round(100*data_ATAC[data_ATAC[,1]==expo_RNA[i],][7]/ATAC_M, 2)

    per_ATAC_M_Mexpo_up<- round(100*data_ATAC[data_ATAC[,1]==expo_RNA[i],][8]/ATAC_M, 2)
    per_ATAC_M_Mexpo_dn<- round(100*data_ATAC[data_ATAC[,1]==expo_RNA[i],][9]/ATAC_M, 2)

    data_for_plot_F_expoF[6*i-3,]<- c(paste0(expo_RNA[i], "_ATAC"), per_ATAC_F_Fexpo_up, "up", as.numeric(data_ATAC[data_ATAC[,1]==expo_RNA[i],][2]));
    data_for_plot_F_expoF[6*i-2,]<- c(paste0(expo_RNA[i], "_ATAC"), per_ATAC_F_Fexpo_dn, "dn", as.numeric(data_ATAC[data_ATAC[,1]==expo_RNA[i],][3]));

    data_for_plot_F_expoM[6*i-3,]<- c(paste0(expo_RNA[i], "_ATAC"), per_ATAC_F_Mexpo_up, "up", as.numeric(data_ATAC[data_ATAC[,1]==expo_RNA[i],][4]));
    data_for_plot_F_expoM[6*i-2,]<- c(paste0(expo_RNA[i], "_ATAC"), per_ATAC_F_Mexpo_dn, "dn", as.numeric(data_ATAC[data_ATAC[,1]==expo_RNA[i],][5]));

    data_for_plot_M_expoF[6*i-3,]<- c(paste0(expo_RNA[i], "_ATAC"), per_ATAC_M_Fexpo_up, "up", as.numeric(data_ATAC[data_ATAC[,1]==expo_RNA[i],][6]));
    data_for_plot_M_expoF[6*i-2,]<- c(paste0(expo_RNA[i], "_ATAC"), per_ATAC_M_Fexpo_dn, "dn", as.numeric(data_ATAC[data_ATAC[,1]==expo_RNA[i],][7]));

    data_for_plot_M_expoM[6*i-3,]<- c(paste0(expo_RNA[i], "_ATAC"), per_ATAC_M_Mexpo_up, "up", as.numeric(data_ATAC[data_ATAC[,1]==expo_RNA[i],][8]));
    data_for_plot_M_expoM[6*i-2,]<- c(paste0(expo_RNA[i], "_ATAC"), per_ATAC_M_Mexpo_dn, "dn", as.numeric(data_ATAC[data_ATAC[,1]==expo_RNA[i],][9]));

#######
    WGBS_F<- as.numeric(strsplit(colnames(data_WGBS)[1], "_")[[1]][2])
    WGBS_M<- as.numeric(strsplit(colnames(data_WGBS)[1], "_")[[1]][4])

    per_WGBS_F_Fexpo_up<- round(100*data_WGBS[data_WGBS[,1]==expo_RNA[i],][2]/WGBS_F, 2)
    per_WGBS_F_Fexpo_dn<- round(100*data_WGBS[data_WGBS[,1]==expo_RNA[i],][3]/WGBS_F, 2)

    per_WGBS_F_Mexpo_up<- round(100*data_WGBS[data_WGBS[,1]==expo_RNA[i],][4]/WGBS_F, 2)
    per_WGBS_F_Mexpo_dn<- round(100*data_WGBS[data_WGBS[,1]==expo_RNA[i],][5]/WGBS_F, 2)

    per_WGBS_M_Fexpo_up<- round(100*data_WGBS[data_WGBS[,1]==expo_RNA[i],][6]/WGBS_M, 2)
    per_WGBS_M_Fexpo_dn<- round(100*data_WGBS[data_WGBS[,1]==expo_RNA[i],][7]/WGBS_M, 2)

    per_WGBS_M_Mexpo_up<- round(100*data_WGBS[data_WGBS[,1]==expo_RNA[i],][8]/WGBS_M, 2)
    per_WGBS_M_Mexpo_dn<- round(100*data_WGBS[data_WGBS[,1]==expo_RNA[i],][9]/WGBS_M, 2)

    data_for_plot_F_expoF[6*i-1,]<- c(paste0(expo_RNA[i], "_WGBS"), per_WGBS_F_Fexpo_up, "up", as.numeric(data_WGBS[data_WGBS[,1]==expo_RNA[i],][2]));
    data_for_plot_F_expoF[6*i,]<- c(paste0(expo_RNA[i], "_WGBS"), per_WGBS_F_Fexpo_dn, "dn", as.numeric(data_WGBS[data_WGBS[,1]==expo_RNA[i],][3]));

    data_for_plot_F_expoM[6*i-1,]<- c(paste0(expo_RNA[i], "_WGBS"), per_WGBS_F_Mexpo_up, "up", as.numeric(data_WGBS[data_WGBS[,1]==expo_RNA[i],][4]));
    data_for_plot_F_expoM[6*i,]<- c(paste0(expo_RNA[i], "_WGBS"), per_WGBS_F_Mexpo_dn, "dn", as.numeric(data_WGBS[data_WGBS[,1]==expo_RNA[i],][5]));

    data_for_plot_M_expoF[6*i-1,]<- c(paste0(expo_RNA[i], "_WGBS"), per_WGBS_M_Fexpo_up, "up", as.numeric(data_WGBS[data_WGBS[,1]==expo_RNA[i],][6]));
    data_for_plot_M_expoF[6*i,]<- c(paste0(expo_RNA[i], "_WGBS"), per_WGBS_M_Fexpo_dn, "dn", as.numeric(data_WGBS[data_WGBS[,1]==expo_RNA[i],][7]));

    data_for_plot_M_expoM[6*i-1,]<- c(paste0(expo_RNA[i], "_WGBS"), per_WGBS_M_Mexpo_up, "up", as.numeric(data_WGBS[data_WGBS[,1]==expo_RNA[i],][8]));
    data_for_plot_M_expoM[6*i,]<- c(paste0(expo_RNA[i], "_WGBS"), per_WGBS_M_Mexpo_dn, "dn", as.numeric(data_WGBS[data_WGBS[,1]==expo_RNA[i],][9]));
   }

data_for_plot_F_expoF$type<- paste0(sapply(strsplit(data_for_plot_F_expoF$expo, "_"), tail,1), "_", data_for_plot_F_expoF$direction)
data_for_plot_F_expoM$type<- paste0(sapply(strsplit(data_for_plot_F_expoM$expo, "_"), tail,1), "_", data_for_plot_F_expoM$direction)

data_for_plot_M_expoF$type<- paste0(sapply(strsplit(data_for_plot_M_expoF$expo, "_"), tail,1), "_", data_for_plot_M_expoF$direction)
data_for_plot_M_expoM$type<- paste0(sapply(strsplit(data_for_plot_M_expoM$expo, "_"), tail,1), "_", data_for_plot_M_expoM$direction)

data_for_plot_F_expoF$expo<- factor(data_for_plot_F_expoF$expo, levels=paste0(rep(expo_RNA, each=3), c("_RNA", "_ATAC", "_WGBS")))
data_for_plot_F_expoM$expo<- factor(data_for_plot_F_expoM$expo, levels=paste0(rep(expo_RNA, each=3), c("_RNA", "_ATAC", "_WGBS")))

data_for_plot_M_expoF$expo<- factor(data_for_plot_M_expoF$expo, levels=paste0(rep(expo_RNA, each=3), c("_RNA", "_ATAC", "_WGBS")))
data_for_plot_M_expoM$expo<- factor(data_for_plot_M_expoM$expo, levels=paste0(rep(expo_RNA, each=3), c("_RNA", "_ATAC", "_WGBS")))


data_for_plot_F_expoF$direction<- factor(data_for_plot_F_expoF$direction, levels=c("up", "dn"))
data_for_plot_F_expoM$direction<- factor(data_for_plot_F_expoM$direction, levels=c("up", "dn"))
data_for_plot_M_expoF$direction<- factor(data_for_plot_M_expoF$direction, levels=c("up", "dn"))
data_for_plot_M_expoM$direction<- factor(data_for_plot_M_expoM$direction, levels=c("up", "dn"))

p<- ggplot(data_for_plot_F_expoF, aes(x= expo, y= ifelse(direction=="up", percent, -percent), fill=type)) + geom_bar(stat="identity", position="identity");
p<- p+ scale_y_continuous(limits=c(-21,21));
p<- p+ theme_classic();
p<- p+ theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10), plot.title = element_text(hjust = 0.5))
p<- p+ geom_vline(size=0.3, color="gray", xintercept=0.5+seq(3,24,3))
p<- p+ labs(x="F Exposure", y="Percent of expo signature among sex signature", title="Liver adult F sexDEG")
p<- p+ scale_fill_manual(values=c("RNA_up"="#E8ACB7","RNA_dn"= "#E16A86", "ATAC_up"="#fdf1df", "ATAC_dn"="#ddaa56", "WGBS_up"="#fedebe", "WGBS_dn"="#fe6e00"));
p<- p+ geom_text(aes(label=count, y=ifelse(direction=="up",percent+1, (-1)*percent-1)), position=position_dodge(0.1), vjust=0.5, angle=90) 
####
q<- ggplot(data_for_plot_F_expoM, aes(x= expo, y= ifelse(direction=="up", percent, -percent), fill=type)) + geom_bar(stat="identity", position="identity");
q<- q+ scale_y_continuous(limits=c(-21,21));
q<- q+ theme_classic();
q<- q+ theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10), plot.title = element_text(hjust = 0.5))
q<- q+ geom_vline(size=0.3, color="gray", xintercept=0.5+seq(3,24,3))
q<- q+ labs(x="M Exposure", y="Percent of expo signature among sex signature", title="Liver adult F sexDEG")
q<- q+ scale_fill_manual(values=c("RNA_up"="#E8ACB7","RNA_dn"= "#E16A86", "ATAC_up"="#fdf1df", "ATAC_dn"="#ddaa56", "WGBS_up"="#fedebe", "WGBS_dn"="#fe6e00"));
q<- q+ geom_text(aes(label=count, y=ifelse(direction=="up",percent+1, (-1)*percent-0.1)), position=position_dodge(0.1), vjust=0.5, angle=90) 
####
r<- ggplot(data_for_plot_M_expoF, aes(x= expo, y= ifelse(direction=="up", percent, -percent), fill=type)) + geom_bar(stat="identity", position="identity");
r<- r+ scale_y_continuous(limits=c(-21,21));
r<- r+ theme_classic();
r<- r+ theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10), plot.title = element_text(hjust = 0.5))
r<- r+ geom_vline(size=0.3, color="gray", xintercept=0.5+seq(3,24,3))
r<- r+ labs(x="F Exposure", y="Percent of expo signature among sex signature", title="Liver adult M sexDEG")
r<- r+ scale_fill_manual(values=c("RNA_up"="#e6eff2", "RNA_dn"="#68a0b4", "ATAC_up"="#EEE6EA", "ATAC_dn"="#AC8295", "WGBS_up"="#B2BAE5","WGBS_dn"= "#768BE6"));
r<- r+ geom_text(aes(label=count, y=ifelse(direction=="up",percent+1, (-1)*percent-1)), position=position_dodge(0.1), vjust=0.5, angle=90)
####
s<- ggplot(data_for_plot_M_expoM, aes(x= expo, y= ifelse(direction=="up", percent, -percent), fill=type)) + geom_bar(stat="identity", position="identity");
s<- s+ scale_y_continuous(limits=c(-21,21));
s<- s+ theme_classic();
s<- s+ theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10), plot.title = element_text(hjust = 0.5))
s<- s+ geom_vline(size=0.3, color="gray", xintercept=0.5+seq(3,24,3))
s<- s+ labs(x="M Exposure", y="Percent of expo signature among sex signature", title="Liver adult M sexDEG")
s<- s+ scale_fill_manual(values=c("RNA_up"="#e6eff2", "RNA_dn"="#68a0b4", "ATAC_up"="#EEE6EA", "ATAC_dn"="#AC8295", "WGBS_up"="#B2BAE5","WGBS_dn"= "#768BE6"));
s<- s+ geom_text(aes(label=count, y=ifelse(direction=="up",percent+1, (-1)*percent-0.1)), position=position_dodge(0.1), vjust=0.5, angle=90)

pdf("bidirectional_barplot_sex_overlap_expo_omics.pdf", width=16, height=14);
grid.arrange(p, q, r, s, ncol=2)
dev.off();

