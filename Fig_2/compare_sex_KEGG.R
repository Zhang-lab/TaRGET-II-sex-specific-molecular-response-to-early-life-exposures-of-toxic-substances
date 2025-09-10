### 2025-04-29 add two terms: Lipid and atheroslerosis and AGE-RAGE
library(ggplot2)
library(stringr)
library(gridExtra)
library(pheatmap)

total_out<- list();
total_count<- 1;

shared_out<- list();
shared_count<- 1;

####
compare_sex<- function(my_stage){

file_path<- paste0("../", my_stage, "_KEGG/")
file_list<- list.files(file_path, pattern="^KEGG_*")
file_list_F<- file_list[grep("Female", file_list)]

qval_cutoff<- 0.1
count_cutoff<- 5;

list_out<- list();
list_count<- 1;

for(i in 1:length(file_list_F)){
    expo<- strsplit(file_list_F[i], "_")[[1]][2]
    lab<- strsplit(file_list_F[i], "_")[[1]][8]
    if(lab %in% c("Biswal", "Mutlu")){
       expo<- paste0(expo, "_", lab);
      }
    data_F<- read.table(paste0(file_path, file_list_F[i]), header=T, sep="\t");
    data_F<- data_F[as.numeric(data_F$count)>=count_cutoff&as.numeric(data_F$lowest_p)<qval_cutoff,]
    data_M<- read.table(paste0(file_path, gsub("Female", "Male", file_list_F[i])), header=T, sep="\t");
    data_M<- data_M[as.numeric(data_M$count)>=count_cutoff&as.numeric(data_M$lowest_p)<qval_cutoff,]
    data_F<- data_F[order(as.numeric(data_F$lowest_p)),]
    data_M<- data_M[order(as.numeric(data_M$lowest_p)),]
    if(nrow(data_F)>0){
       for(j in 1:nrow(data_F)){
           total_out[[total_count]]<<- c(my_stage, "F", expo, data_F[j, "Term_Description"], data_F[j, "lowest_p"], data_F[j, "count"]); 
           total_count<<- total_count +1;
          }
      }
    if(nrow(data_M)>0){
       for(j in 1:nrow(data_M)){
           total_out[[total_count]]<<- c(my_stage, "M", expo, data_M[j, "Term_Description"], data_M[j, "lowest_p"], data_M[j, "count"]);
           total_count<<- total_count +1;
          }
      }   
    shared_terms<- intersect(data_F$Term_Description, data_M$Term_Description)
    if(length(shared_terms)>0){
       for(j in 1:length(shared_terms)){
           shared_out[[shared_count]]<<- c(my_stage, expo, shared_terms[j], as.numeric(data_F[data_F$Term_Description== shared_terms[j],c("lowest_p", "count")]), as.numeric(data_M[data_M$Term_Description== shared_terms[j],c("lowest_p", "count")]))
           shared_count<<- shared_count +1;
          }
      }
    F_terms<- data_F[!(data_F$ID %in% data_M$ID), "Term_Description"]
    M_terms<- data_M[!(data_M$ID %in% data_F$ID), "Term_Description"]
    list_out[[list_count]]<- c(expo, my_stage, "shared", length(shared_terms));
    list_count<- list_count +1;
    list_out[[list_count]]<- c(expo, my_stage, "Fspecific", length(F_terms));
    list_count<- list_count +1;
    list_out[[list_count]]<- c(expo, my_stage, "Mspecific", length(M_terms));
    list_count<- list_count +1;

    term_cutoff<- 8;
    uniq_F<- data_F[data_F$Term_Description %in% F_terms,]; 
    uniq_M<- data_M[data_M$Term_Description %in% M_terms,]; 
    if(nrow(uniq_F)>term_cutoff){
       uniq_F<- uniq_F[1:term_cutoff,]
      }
    if(nrow(uniq_M)>term_cutoff){
       uniq_M<- uniq_M[1:term_cutoff,]
      }
    if(nrow(uniq_F)==0 & nrow(uniq_M)==0){
       next;
      }
    data_for_plot<- rbind.data.frame(uniq_F[,c("Term_Description", "lowest_p", "count")], uniq_M[,c("Term_Description", "lowest_p", "count")])
    data_for_plot$sex<- c(rep("F", nrow(uniq_F)), rep("M", nrow(uniq_M)));
    data_for_plot$logQ<- (-1)*log10(as.numeric(data_for_plot$lowest_p))
    data_for_plot$count<- as.numeric(data_for_plot$count)
    data_for_plot$Term_Description<- factor(data_for_plot$Term_Description, levels=rev(data_for_plot$Term_Description))
print(i); 
    p<- ggplot(data_for_plot, aes(x=sex, y=Term_Description, color=logQ, size=count) ) + geom_point() + scale_color_gradient(low = "lightblue", high = "darkblue")
   p<- p + theme_bw() + ylab("") + xlab("") + ggtitle(paste0("KEGG uniq sex ", my_stage," ", expo));
   p<- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
   p<- p + scale_y_discrete(labels = function(n) str_wrap(n, width = 40))
   pdf(paste0("uniq_KEGG_plots/KEGG_uniq_sex_top", term_cutoff, "_", my_stage, "_", expo, ".pdf"), height=5.5, width=5);
   print(p);
   dev.off();
   }
data_out<- data.frame(do.call(rbind, list_out))
colnames(data_out)<- c("expo", "stage", "type", "count");
return(data_out);
}
####

data_wl<- compare_sex("weanling");
data_adt<- compare_sex("adult");
data_aged<- compare_sex("aged");

data_combined<- rbind.data.frame(data_wl, data_adt, data_aged);
data_combined$count<- as.numeric(data_combined$count)
count_max<- ceiling(max(data_combined$count)/10)*10

#####
get_plot<- function(my_data, my_type){

my_data<- my_data[my_data$type== my_type,]
my_data$stage<- factor(my_data$stage, levels=c("weanling", "adult", "aged"));
my_data$expo<- factor(my_data$expo, levels=rev(c("As", "BPA10ug", "BPA10mg", "DEHP", "Pb", "PM2.5_Biswal", "PM2.5_Mutlu", "TBT", "TCDD")))

p<- ggplot(my_data, aes(x=stage, y=expo, color=count, size=count)) + geom_point()
p<- p+ scale_color_gradient(low = "lightblue", high = "darkblue", limits=c(1, count_max))
p<- p+ scale_size(limits=c(1, count_max))
p<- p+ theme_bw()
p<- p+ ggtitle(paste0("KEGG ", my_type));

return(p);
}
#####

p1<- get_plot(data_combined, "Fspecific")
p2<- get_plot(data_combined, "Mspecific")
p3<- get_plot(data_combined, "shared")

pdf("KEGG_count_across_stages.pdf", width=8, height=3);
print(grid.arrange(p1,p2, p3, nrow=1));
dev.off();
############
total_data<- data.frame(do.call(rbind, total_out))
colnames(total_data)<- c("stage", "sex", "expo", "term", "pval", "count")

shared_data<- data.frame(do.call(rbind, shared_out));
colnames(shared_data)<- c("stage", "expo", "term", "F_p", "F_count", "M_p", "M_count");
shared_freq<- data.frame(table(shared_data$term))
shared_freq<- shared_freq[order(-shared_freq$Freq),]

freq_cutoff<- 3;
shared_terms<- shared_freq[shared_freq$Freq>= freq_cutoff, 1]

shared_cond<- unique(shared_data[shared_data$term %in% shared_terms,c(1,2)])
data_heatmap<- data.frame(matrix(0, nrow=length(shared_terms), ncol= 2*nrow(shared_cond)))
colnames(data_heatmap)<- paste0(rep(paste0(shared_cond[,1], "_", shared_cond[,2]),each=2), c("_F", "_M"))
rownames(data_heatmap)<- shared_terms;

for(i in 1:nrow(data_heatmap)){
    for(j in 1:nrow(shared_cond)){
        sub_data<- total_data[total_data$term== rownames(data_heatmap)[i] & total_data$stage== shared_cond[j, "stage"] & total_data$expo== shared_cond[j, "expo"],]
        if(nrow(sub_data)==0){
           next;
          }
        for(k in 1:nrow(sub_data)){
            tag<- paste(sub_data[k, "stage"], sub_data[k, "expo"], sub_data[k, "sex"], sep="_")
            data_heatmap[i, tag]<- (-1)* log10(as.numeric(sub_data[k, "pval"]));
           }
       }
   }
data_heatmap[data_heatmap>10]<- 10;

row_ann<- data.frame(shared_freq[shared_freq$Freq>= freq_cutoff,"Freq"])
colnames(row_ann)<- "freq"
rownames(row_ann)<- shared_terms;

freq_count<- data.frame(table(shared_freq[shared_freq$Freq>= freq_cutoff,"Freq"]))
for_gap_row<- cumsum(rev(freq_count$Freq))
for_gap_row<- for_gap_row[1:(length(for_gap_row)-1)]

pdf("heatmap_KEGG_shared_by_sex.pdf", width=9);
pheatmap(data_heatmap, cluster_cols=F, cluster_rows=F, scale="none", annotation_row=row_ann, gaps_col=2*c(1:(nrow(shared_cond)-1)), gaps_row= for_gap_row, main=paste0("-logadjP shared_terms n=", nrow(data_heatmap)), 
color=colorRampPalette(c("white", "red"))(100) )
dev.off();
#####

total_data_F<- total_data[total_data$sex=="F",]
total_data_M<- total_data[total_data$sex=="M",]

#### 
get_heatmap_by_sex<- function(my_sex, my_freq_cutoff){

my_data<- total_data[total_data$sex== my_sex,]
term_freq<- data.frame(table(my_data$term))
term_freq<- term_freq[order(-term_freq$Freq),]
term_freq<- term_freq[term_freq$Freq>= my_freq_cutoff,]

data_for_plot<- data.frame(matrix(0, ncol=24, nrow=nrow(term_freq)));
rownames(data_for_plot)<- term_freq[,1]

colnames(data_for_plot)<- c("As_weanling", "As_adult", "BPA10ug_weanling", "BPA10ug_adult","BPA10ug_aged", "BPA10mg_weanling", "BPA10mg_adult", "BPA10mg_aged", "DEHP_weanling", "DEHP_adult", "DEHP_aged", "Pb_weanling", "Pb_adult", "Pb_aged", "PM2.5_Biswal_weanling", "PM2.5_Biswal_adult", "PM2.5_Mutlu_weanling", "PM2.5_Mutlu_adult","PM2.5_Mutlu_aged", "TBT_weanling", "TBT_adult","TBT_aged", "TCDD_weanling", "TCDD_adult")

for(i in 1:nrow(my_data)){

    if(my_data[i, "term"] %in% rownames(data_for_plot)){
       tag<- paste0(my_data[i, "expo"],"_", my_data[i, "stage"])
       data_for_plot[my_data[i, "term"], tag]<-  (-1)*log10(as.numeric(my_data[i,"pval"]))
      }
   }
data_for_plot[data_for_plot>10]<- 10;

row_ann<- data.frame(term_freq[,"Freq"])
colnames(row_ann)<- "freq"
rownames(row_ann)<- rownames(data_for_plot);

freq_count<- data.frame(table(term_freq[,"Freq"]))
for_gap_row<- cumsum(rev(freq_count$Freq))
for_gap_row<- for_gap_row[1:(length(for_gap_row)-1)]
for_gaps_col<- c(2, 5, 8, 11, 14, 16, 19, 22)

pdf(paste0("heatmap_freq_KEGG_", my_sex, ".pdf"), width=8, height=10);
print(pheatmap(data_for_plot, cluster_cols=F, cluster_rows=F, scale="none", annotation_row=row_ann, gaps_col=for_gaps_col, gaps_row= for_gap_row, main=paste0(my_sex," -logadjP terms n=", nrow(data_for_plot)),
color=colorRampPalette(c("white", "red"))(100) ) )
dev.off();
}
#### 

get_heatmap_by_sex("F", 6)
get_heatmap_by_sex("M", 7)

####################
get_freq_term_by_stage<- function(my_data, my_stage, my_cutoff){

my_data_stage<- my_data[my_data$stage== my_stage,]
stage_term_freq<- data.frame(table(my_data_stage$term))
stage_term_freq<- stage_term_freq[order(-stage_term_freq$Freq),]
selected_terms<- as.character(stage_term_freq[stage_term_freq$Freq>= my_cutoff,"Var1"])
excluded_terms<- c("Breast cancer", "Protein processing in endoplasmic reticulum", "Regulation of actin cytoskeleton", "Progesterone-mediated oocyte maturation", "Oocyte meiosis", "Kaposi sarcoma-associated herpesvirus infection", "Measles");
selected_terms<- selected_terms[!(selected_terms %in% excluded_terms)]
return(selected_terms);
}
####

terms_wl<- get_freq_term_by_stage(shared_data, "weanling", 2);
terms_adt<- get_freq_term_by_stage(shared_data, "adult", 2);
terms_aged<- get_freq_term_by_stage(shared_data, "aged", 3);
#### add two KEGG pathways
terms_wl<- c(terms_wl, "Lipid and atherosclerosis", "AGE-RAGE signaling pathway in diabetic complications");
terms_adt<- c(terms_adt, "Lipid and atherosclerosis", "AGE-RAGE signaling pathway in diabetic complications");
terms_aged<- c(terms_aged, "Lipid and atherosclerosis", "AGE-RAGE signaling pathway in diabetic complications");


terms_combined<- c(terms_wl, terms_adt, terms_aged)

data_combined<- data.frame(matrix(0, nrow=length(terms_combined), ncol=9));
rownames(data_combined)<- c(paste0(terms_wl, "_wl"), paste0(terms_adt, "_adt"), paste0(terms_aged, "_aged"))
colnames(data_combined)<- c("As", "BPA10ug", "BPA10mg", "DEHP", "Pb", "PM2.5_Biswal", "PM2.5_Mutlu", "TBT", "TCDD")

for(i in 1:nrow(data_combined)){
    stage<- sapply(strsplit(rownames(data_combined)[i],"_"), tail, 1)
    term<- sapply(strsplit(rownames(data_combined)[i],"_"), head, 1)
    if(stage=="wl"){
       stage<- "weanling"
      }
    if(stage=="adt"){
       stage<- "adult"
      }
    for(j in 1:ncol(data_combined)){
        sub_data<- total_data[total_data$stage==stage & total_data$term== term & total_data$expo== colnames(data_combined)[j],]
        if(nrow(sub_data)==2){
           data_combined[i, j]<- 3;
          }
        if(nrow(sub_data)==1 & sub_data[1,"sex"]== "F"){
           data_combined[i, j]<- 1;
          }
        if(nrow(sub_data)==1 & sub_data[1,"sex"]== "M"){
           data_combined[i, j]<- 2;
          }
       }
   }
data_combined_adt<- data_combined[grep("_adt", rownames(data_combined)),]
data_combined_aged<- data_combined[grep("_aged", rownames(data_combined)),]
data_combined_aged[,c("As", "PM2.5_Biswal", "TCDD")]<- (-1);
data_combined1<- rbind.data.frame(data_combined[order(-rowSums(data_combined[1:length(terms_wl),]!=0)),], 
data_combined_adt[order(-rowSums(data_combined_adt!=0)),],
data_combined_aged[order(-rowSums(data_combined_aged!=0)),])


myColors<- colorRampPalette(c("#DCDCDC", "white", "#E16A86", "#768BE6", "#92a881"))(5);
myBreaks<- c(-1.5, -0.5, 0.5, 1.5, 2.5, 3.5)
gaps_for_row<- c(length(terms_wl), length(terms_wl)+length(terms_adt));

row_ann<- data.frame(c(rep("wl", length(terms_wl)), rep("adt",length(terms_adt)), rep("aged", length(terms_aged))));
colnames(row_ann)<- "stage"
rownames(row_ann)<- rownames(data_combined1);
row_ann$stage<- factor(row_ann$stage, levels=c("wl", "adt", "aged"));
### Remove some terms before plotting
data_combined1<- data_combined1[!grepl("Epstein-Barr virus infection_aged", rownames(data_combined1)),]


pdf("heatmap_KEGG_freq_shared_by_sex.pdf", width=7, height=6.5);
pheatmap(data_combined1, cluster_rows=F, cluster_cols=F, breaks=myBreaks, color= myColors, scale="none", gaps_row=gaps_for_row, annotation_row=row_ann, main=paste0("Freq shared KEGG by sex n=", nrow(data_combined1)));
dev.off();














