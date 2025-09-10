#Only hpc2 successful 2022-06-09 4:59pm
#.libPaths("/home/shuhua/R/x86_64-pc-linux-gnu-library/4.0")
#.libPaths("/usr/lib/R/site-library/") #hpc3


library(ggalluvial)
library(tidyr)
#library(rtracklayer)

#z<- import("/home/Resource/Genome/mm10/gencode.vM20.annotation.gtf")
#dat_z<- data.frame(z)

#dat_gene<- dat_z[dat_z$type=="gene",]
#chrXY_genes<- unique(dat_gene[dat_gene$seqnames %in% c("chrX", "chrY"),"gene_name"])
chrXY_data<- read.table("/BRC/shuhua/target/RNAseq/DEGs/mm10.v20_chrXY_genes", header=F);
chrXY_genes<- chrXY_data[,1]

data_liver_wl_up<-read.table("DEGs/FvsM_DEG_list/deg_liver_wl_up", header=T, row.names=1)
data_liver_wl_dn<-read.table("DEGs/FvsM_DEG_list/deg_liver_wl_dn", header=T, row.names=1)

data_liver_adt_up<-read.table("DEGs/FvsM_DEG_list/deg_liver_adt_up", header=T, row.names=1)
data_liver_adt_dn<-read.table("DEGs/FvsM_DEG_list/deg_liver_adt_dn", header=T, row.names=1)

data_liver_aged_up<-read.table("DEGs/FvsM_DEG_list/deg_liver_aged_up", header=T, row.names=1)
data_liver_aged_dn<-read.table("DEGs/FvsM_DEG_list/deg_liver_aged_dn", header=T, row.names=1)

deg_liver_wl_up<- rownames(data_liver_wl_up)
deg_liver_wl_dn<- rownames(data_liver_wl_dn)

deg_liver_adt_up<- rownames(data_liver_adt_up)
deg_liver_adt_dn<- rownames(data_liver_adt_dn)

deg_liver_aged_up<- rownames(data_liver_aged_up)
deg_liver_aged_dn<- rownames(data_liver_aged_dn)

degs<- unique(c(deg_liver_wl_up, deg_liver_wl_dn, deg_liver_adt_up, deg_liver_adt_dn, deg_liver_aged_up, deg_liver_aged_dn))
################# Remove chrXY
degs<- degs[!(degs %in% chrXY_genes)]

wl_stage<- rep("Stably expressed", length(degs) ); 
wl_stage[degs %in% deg_liver_wl_up]<- "Up-regulated"
wl_stage[degs %in% deg_liver_wl_dn]<- "Down-regulated"

adt_stage<- rep("Stably expressed", length(degs) );
adt_stage[degs %in% deg_liver_adt_up]<- "Up-regulated"
adt_stage[degs %in% deg_liver_adt_dn]<- "Down-regulated"

aged_stage<- rep("Stably expressed", length(degs) );
aged_stage[degs %in% deg_liver_aged_up]<- "Up-regulated"
aged_stage[degs %in% deg_liver_aged_dn]<- "Down-regulated"

data_stage<- cbind.data.frame(wl_stage, adt_stage, aged_stage)
colnames(data_stage)<- c("Weaning", "Adult", "Aged");
rownames(data_stage)<- degs;

data_to_plot<- data.frame(data_stage %>% tibble::rownames_to_column(var="Row") %>% pivot_longer(-Row, names_to="Stage", values_to=("Expression") ))

data_to_plot$Expression<- factor(data_to_plot$Expression, levels=c("Up-regulated", "Down-regulated", "Stably expressed"))
data_to_plot$Stage<- factor(data_to_plot$Stage, levels=c("Weaning", "Adult", "Aged") );

p<- ggplot(data_to_plot, aes(x=Stage, stratum=Expression, alluvium=Row, label=Expression, fill=Expression))
p<- p+ scale_x_discrete(expand=c(0.05,0.05))
#q<- q+ geom_alluvium(aes(fill=Expression, width=1/8))

p<- p+ geom_flow(stat="alluvium", width=0.1)  + geom_stratum(alpha= 0.5)
p<- p+ scale_y_continuous(expand=c(0,0) )



p<- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), plot.title=element_text(hjust=0.5, size=22, face="bold"), axis.text=element_text(size=15), axis.title=element_text(size=20), legend.text=element_text(size=12), legend.title=element_text(size=18) )

p<- p+ labs(y="Gene count", x="Stage", title=paste0("Gene expression across stages in liver\nn=", length(degs)) );

p<- p+ scale_fill_manual(values=c("#E69F00", "#56B4E9","#999999"))


pdf("sex_deg_liver_stages_alluvial.pdf");
p;
dev.off();

###############

for(i in 1:length(degs)){
    d<- data_to_plot[data_to_plot$Row==degs[i],];
    if(length(unique(d$Expression))==1 & d[1,"Expression"]=="Stably expressed"){
       print(i); 
      }
  }






###### get plot for weaning & adult tissue
get_plot<- function(my_tissue){

data_wl_up<-read.table(paste0("DEGs/FvsM_DEG_list/deg_", my_tissue,"_wl_up"), header=T, row.names=1)
data_wl_dn<-read.table(paste0("DEGs/FvsM_DEG_list/deg_", my_tissue,"_wl_dn"), header=T, row.names=1)

data_adt_up<-read.table(paste0("DEGs/FvsM_DEG_list/deg_", my_tissue,"_adt_up"), header=T, row.names=1)
data_adt_dn<-read.table(paste0("DEGs/FvsM_DEG_list/deg_", my_tissue,"_adt_dn"), header=T, row.names=1)

deg_wl_up<- rownames(data_wl_up);
deg_wl_dn<- rownames(data_wl_dn);

deg_adt_up<- rownames(data_adt_up);
deg_adt_dn<- rownames(data_adt_dn);

DEGs<- unique(c(deg_wl_up, deg_wl_dn, deg_adt_up, deg_adt_dn));
##### Remove chrXY
DEGs<- DEGs[!(DEGs %in% chrXY_genes)]

Weaning<-  rep("Stably expressed", length(DEGs));
Weaning[DEGs %in% deg_wl_up]<- "Up-regulated";
Weaning[DEGs %in% deg_wl_dn]<- "Down-regulated";

Adult<-  rep("Stably expressed", length(DEGs));
Adult[DEGs %in% deg_adt_up]<- "Up-regulated";
Adult[DEGs %in% deg_adt_dn]<- "Down-regulated";

Data_stage<- cbind.data.frame(Weaning, Adult);
rownames(Data_stage)<- DEGs;

Data_to_plot<- data.frame(Data_stage %>% tibble::rownames_to_column(var="Row") %>% pivot_longer(-Row, names_to="Stage", values_to=("Expression") ))

Data_to_plot$Expression<- factor(Data_to_plot$Expression, levels=c("Up-regulated", "Down-regulated", "Stably expressed"))
Data_to_plot$Stage<- factor(Data_to_plot$Stage, levels=c("Weaning", "Adult") );

q<- ggplot(Data_to_plot, aes(x=Stage, stratum=Expression, alluvium=Row, label=Expression, fill=Expression))
q<- q+ scale_x_discrete(expand=c(0.05,0.05))
#q<- q+ geom_alluvium(aes(fill=Expression, width=1/8))

q<- q+ geom_flow(stat="alluvium", width=0.1)  + geom_stratum(alpha= 0.5) 
q<- q+ scale_y_continuous(expand=c(0,0) ) 
q<- q + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), plot.title=element_text(hjust=0.5, size=22, face="bold"), axis.text=element_text(size=15), axis.title=element_text(size=20), legend.text=element_text(size=12), legend.title=element_text(size=18) )

q<- q+ labs(y="Gene count", x="Stage", title=paste0("Gene expression across stages in ", my_tissue,"\nn=", length(DEGs)) );
q<- q+ scale_fill_manual(values=c("#E69F00", "#56B4E9","#999999"))

pdf(paste0("sex_deg_", my_tissue, "_stages_alluvial.pdf") ); 
print(q);
dev.off();
}
###########

get_plot("lung");
get_plot("blood");




