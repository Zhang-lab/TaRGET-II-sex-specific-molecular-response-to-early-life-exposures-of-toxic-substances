library(ggalluvial)
library(tidyr)
library(gridExtra)

path_wl<- "../../treatment_DEG_list_noChrXY/weanling_DEG/";
path_adt<- "../../treatment_DEG_list_noChrXY/adult_DEG/";
path_aged<- "../../treatment_DEG_list_noChrXY/aged_DEG/";

list_wl<- list.files(path_wl);
list_wl<- list_wl[grep("liver", list_wl)]
list_adt<- list.files(path_adt);
list_adt<- list_adt[grep("liver", list_adt)]
list_aged<- list.files(path_aged);
list_aged<- list_aged[grep("liver", list_aged)]

data_info<- unique(data.frame(t(data.frame(strsplit(list_wl, "_"))))[,c(1,7)])

######
get_plot<- function(my_expo, my_lab, my_sex){

sub_list_wl<- list_wl[grepl(my_lab, list_wl) & grepl(my_expo, list_wl) & grepl(my_sex, list_wl) ]
sub_list_adt<- list_adt[grepl(my_lab, list_adt) & grepl(my_expo, list_adt) & grepl(my_sex, list_adt) ]
sub_list_aged<- list_aged[grepl(my_lab, list_aged) & grepl(my_expo, list_aged) & grepl(my_sex, list_aged) ]

out_list<- list();
list_count<- 1;

degs<- c();
for(i in 1:length(sub_list_wl)){
    direction<- sapply(strsplit(sub_list_wl[i], "_"), tail, 1)
    data<- read.table(paste0(path_wl, sub_list_wl[i]), header=T, sep="\t");
    for(j in 1:nrow(data)){
        out_list[[list_count]]<- c(data[j, "gene"], direction, "wl");
        list_count<- list_count +1;
       }
   }
for(i in 1:length(sub_list_adt)){
    direction<- sapply(strsplit(sub_list_adt[i], "_"), tail, 1)
    data<- read.table(paste0(path_adt, sub_list_adt[i]), header=T, sep="\t");
    for(j in 1:nrow(data)){
        out_list[[list_count]]<- c(data[j, "gene"], direction, "adt");
        list_count<- list_count +1;
       }
   }
if(length(sub_list_aged)>0){
   for(i in 1:length(sub_list_aged)){
       direction<- sapply(strsplit(sub_list_aged[i], "_"), tail, 1)
       data<-read.table(paste0(path_aged,sub_list_aged[i]), header=T, sep="\t");
       for(j in 1:nrow(data)){
           out_list[[list_count]]<- c(data[j, "gene"], direction, "aged");
           list_count<- list_count +1;
          }
      }
  }
data_for_plot<- data.frame(do.call(rbind, out_list))
colnames(data_for_plot)<- c("gene", "direction", "stage");

wl_up_adt_up<- length(which(duplicated(data_for_plot[data_for_plot$direction=="up" & data_for_plot$stage %in% c("wl", "adt"),"gene"])))
wl_dn_adt_dn<- length(which(duplicated(data_for_plot[data_for_plot$direction=="down" & data_for_plot$stage %in% c("wl", "adt"),"gene"])))

wl_up_adt_dn<- length(which(duplicated(data_for_plot[(data_for_plot$direction=="up" & data_for_plot$stage=="wl") | (data_for_plot$direction=="down" & data_for_plot$stage=="adt")  ,"gene"])))
wl_dn_adt_up<- length(which(duplicated(data_for_plot[(data_for_plot$direction=="down" & data_for_plot$stage=="wl") | (data_for_plot$direction=="up" & data_for_plot$stage=="adt")  ,"gene"])))

aged_up_adt_up<- length(which(duplicated(data_for_plot[data_for_plot$direction=="up" & data_for_plot$stage %in% c("aged", "adt"),"gene"])))
aged_dn_adt_dn<- length(which(duplicated(data_for_plot[data_for_plot$direction=="down" & data_for_plot$stage %in% c("aged", "adt"),"gene"])))

aged_up_adt_dn<- length(which(duplicated(data_for_plot[(data_for_plot$direction=="up" & data_for_plot$stage=="aged") | (data_for_plot$direction=="down" & data_for_plot$stage=="adt")  ,"gene"])))
aged_dn_adt_up<- length(which(duplicated(data_for_plot[(data_for_plot$direction=="down" & data_for_plot$stage=="aged") | (data_for_plot$direction=="up" & data_for_plot$stage=="adt")  ,"gene"])))



data_for_plot$stage<-factor(data_for_plot$stage, levels=c("wl", "adt", "aged"));
data_for_plot$direction<-factor(data_for_plot$direction,levels=c("up", "down"));

#p<- ggplot(data_for_plot, aes(x=stage, stratum=direction, alluvium=gene, label=direction, fill=direction))
p<- ggplot(data_for_plot, aes(x=stage, stratum=direction, alluvium=gene, label=direction, fill=direction))
p<- p+ scale_x_discrete(expand=c(0.05,0.05))

#p<- p+ geom_flow(stat="alluvium", width=0.1)  + geom_stratum(alpha= 0.5)
p<- p+ geom_flow(aes(fill=direction),  width=0.1)  + geom_stratum(alpha= 0.5)
p<- p+ scale_y_continuous(expand=c(0,0) )

p<- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), plot.title=element_text(hjust=0.5, size=11, face="bold"), axis.text=element_text(size=15), axis.title=element_text(size=18), legend.text=element_text(size=12), legend.title=element_text(size=16) )

p<- p+ labs(y="Gene count", x="Stage", title=paste0("DEGs across stages liver\n", my_expo, " ", my_lab, " n=", length(unique(data_for_plot$gene)), "\nwl_up_adt_up=", wl_up_adt_up, " wl_dn_adt_dn=",wl_dn_adt_dn,"\nwl_up_adt_dn=", wl_up_adt_dn," wl_dn_adt_up=", wl_dn_adt_up,"\naged_up_adt_up=", aged_up_adt_up," aged_dn_adt_dn=", aged_dn_adt_dn,"\naged_up_adt_dn=", aged_up_adt_dn," aged_dn_adt_up=", aged_dn_adt_up ) );
#p<- p+ scale_fill_manual(values=c("#E69F00", "#56B4E9","#999999"))
return(p);
}
######

data_info<- data_info[c(1,3,2,4:9),]
pList<- list();
qList<- list();
for(i in 1:nrow(data_info)){
    pList[[i]]<- get_plot(data_info[i, 1], data_info[i, 2], "Female");
    qList[[i]]<- get_plot(data_info[i, 1], data_info[i, 2], "Male");
   }
pdf("alluvial_expoDEG_across_stages_F.pdf", width=38, height=5);
grid.arrange(grobs=pList, ncol=9)
dev.off();

pdf("alluvial_expoDEG_across_stages_M.pdf", width=38, height=5);
grid.arrange(grobs=qList, ncol=9)
dev.off();

