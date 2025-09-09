library(ggplot2)
library(ggrepel)
library(ggpubr)
library(gridExtra)

shared_files<- list.files(".", pattern="^shared_*")
shared_files_wl<- shared_files[grep("_wl_", shared_files)]
shared_files_adt<- shared_files[grep("_adt_", shared_files)]
shared_files_aged<- shared_files[grep("_aged_", shared_files)]

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

########
get_plot<- function(my_shared_files, my_stage){

data1<- read.table(my_shared_files[1], header=T, sep="\t")
colnames(data1)<- c("name", "F_LFC", "M_LFC")
data1$expo<- gsub("shared_", "", my_shared_files[1]);
data1$expo<- gsub(paste0(my_stage, "_"), "", data1$expo);

for(i in 2:length(my_shared_files)){
    data<- read.table(my_shared_files[i], header=T, sep="\t")
    if(nrow(data)>0){
       colnames(data)<- c("name", "F_LFC", "M_LFC")
       data$expo<- gsub("shared_", "", my_shared_files[i]);
       data$expo<- gsub(paste0(my_stage, "_"), "", data$expo);
       data1<- rbind.data.frame(data1, data);
      }
   }

data1$F_LFC<- as.numeric(as.character(data1$F_LFC))
data1$M_LFC<- as.numeric(as.character(data1$M_LFC))

set_cutoff<- 4
data1[data1[,2]< (-1)*set_cutoff,2]<-  (-1)*set_cutoff
data1[data1[,2]>set_cutoff,2]<- set_cutoff
data1[data1[,3]<(-1)*set_cutoff, 3]<- (-1)* set_cutoff
data1[data1[,3]>set_cutoff, 3]<- set_cutoff
data1$col<- ifelse(data1$F_LFC * data1$M_LFC >0, "same", "diff")

max_val<-  max(abs(data1[,2:3]) )+ 0.0005

data1$expo<- factor(data1$expo, levels=c("As", "BPA10ug", "BPA10mg", "DEHP", "Pb", "PM2.5_Biswal", "PM2.5_Mutlu", "TBT", "TCDD"))
p<- ggplot(data1, aes(x=F_LFC, y=M_LFC)) + geom_point(fill=NA, colour=NA, size=0.7)
p<- p+ theme_classic()
p<- p+ xlim((-1)* max_val, max_val) + ylim((-1)* max_val, max_val)
p<- p+ geom_smooth(method="lm", aes(color=expo))
p<- p+ stat_cor(data=data1[data1$expo!="BPA10ug",],aes(color=expo))
p<- p+ scale_color_manual(values=c("As"="#F8766D", "BPA10ug"="#D39200", "BPA10mg"="#93AA00", "DEHP"="#00BA38", "Pb"="#00C19F", "PM2.5_Biswal"="#00B9E3", "PM2.5_Mutlu"="#619CFF", "TBT"="#DB72FB", "TCDD"="#FF61C3"))
p<- p+ ggtitle(paste0(my_stage," expoDEG shared by sex"));

return(p);
}
########

p_wl<- get_plot(shared_files_wl, "wl")
p_adt<- get_plot(shared_files_adt, "adt")
p_aged<- get_plot(shared_files_aged, "aged")

pdf("trend_lines_shared_expoDEG_stages.pdf", width=18, height=5);
grid.arrange(p_wl, p_adt, p_aged, nrow=1);
dev.off();

