## On hpc5 https://gist.github.com/jokergoo/ac65ed4061c83fe26a4a3e182c9cb6f9
library(circlize)
library(ComplexHeatmap)

list_wl<- list.files("../GSEA_output/out_wl/")
list_adt<- list.files("../GSEA_output/out_adt/")
list_aged<- list.files("../GSEA_output/out_aged/")

#######
extract_data<- function(my_list, my_path, my_stage){

list_out<- list();
list_count<- 1;

for(i in 1:length(my_list)){

    enplot_list<- list.files(paste0(my_path, my_list[i]), pattern="^enplot_CYP_HS_*")  
    if(length(enplot_list)>=0){
print(i);
       report_list<- list.files(paste0(my_path, my_list[i]), pattern="^gsea_*")
       report_list<- report_list[grepl(".tsv", report_list)]
       tag<- gsub("m2_edited_", "", gsub(paste0(".Gsea.", sapply(strsplit(gsub("m2_edited_","",my_list[i]), "[.]"), tail, 1)), "", my_list[i]))
       for(j in 1:length(report_list)){
           data<- read.table(paste0(paste0(my_path, my_list[i]), "/", report_list[j]), header=T, sep="\t")
           if("CYP_HS" %in% data$NAME){
print(j); 
              ES<- data[data$NAME=="CYP_HS","ES"]
              NES<- data[data$NAME=="CYP_HS","NES"]
              pval<- data[data$NAME=="CYP_HS", "NOM.p.val"]
              pval<- data[data$NAME=="CYP_HS", "NOM.p.val"]
              FDR<- data[data$NAME=="CYP_HS", "FDR.q.val"]
              FWER<- data[data$NAME=="CYP_HS", "FWER.p.val"]
              rank<- data[data$NAME=="CYP_HS", "RANK.AT.MAX"]
              size<- data[data$NAME=="CYP_HS", "SIZE"]
              group<- gsub("gsea_report_for_", "", gsub(paste0("_",sapply(strsplit(report_list[j], "_"),tail,1)), "", report_list[j]))
              list_out[[list_count]]<- c(my_stage, tag, group, ES, NES,  pval, FDR, FWER,rank, size);
              list_count<- list_count +1;
             }
          }
      }
   }
data_out<- data.frame(do.call(rbind, list_out))
colnames(data_out)<- c("stage", "compare", "group", "ES", "NES", "pval", "FDR", "FWER", "rank", "size");

return(data_out);
}
#######

out_wl<- extract_data(list_wl, "../GSEA_output/out_wl/", "wl");
out_adt<- extract_data(list_adt, "../GSEA_output/out_adt/", "adt");
out_aged<- extract_data(list_aged, "../GSEA_output/out_aged/", "aged");

data_out<- rbind.data.frame(out_wl, out_adt, out_aged)

data_out$sex<- ifelse(grepl("_M_", data_out$compare), "M", "F")
data_out$expo<- gsub("_M$|_F$", "", sapply(strsplit(data_out$compare, "_vs"), head, 1))

data_out$tag<- paste0(data_out$stage, "_", data_out$expo)

data_out_F<- data_out[data_out$sex=="F",]
data_out_M<- data_out[data_out$sex=="M",]

#######
get_plot<- function(my_data, my_sex){

data_for_plot_NES<- data.frame(matrix(NA, nrow=3, ncol=9))

colnames(data_for_plot_NES)<- c("As", "BPA10ug", "BPA10mg", "DEHP", "Pb", "PM2.5_Biswal", "PM2.5_Mutlu", "TBT","TCDD");
rownames(data_for_plot_NES)<- c("wl", "adt", "aged");
data_for_plot_sig<- data_for_plot_NES;

data_wl<- my_data[my_data$stage=="wl", ]
data_adt<- my_data[my_data$stage=="adt", ]
data_aged<- my_data[my_data$stage=="aged", ]

data_for_plot_NES["wl",]<- as.numeric(data_wl[match(colnames(data_for_plot_NES), data_wl$expo), "NES"])
data_for_plot_NES["adt",]<- as.numeric(data_adt[match(colnames(data_for_plot_NES), data_adt$expo), "NES"])
data_for_plot_NES["aged",]<- as.numeric(data_aged[match(colnames(data_for_plot_NES), data_aged$expo), "NES"])

if(my_sex=="F"){
   data_for_plot_sig["adt", "PM2.5_Mutlu"]<- 1;
   data_for_plot_sig["aged", "DEHP"]<- 1;
   data_for_plot_sig["aged", "Pb"]<- 1;
  }
if(my_sex=="M"){
   data_for_plot_sig["wl", "As"]<- 1;
   data_for_plot_sig["wl", "BPA10mg"]<- 1;
   data_for_plot_sig["adt", "PM2.5_Biswal"]<- 1;
   data_for_plot_sig["adt", "PM2.5_Mutlu"]<- 1;
   data_for_plot_sig["aged", "BPA10mg"]<- 1;
   data_for_plot_sig["aged", "TBT"]<- 1;
  }
data_for_plot_sig[is.na(data_for_plot_sig)]<- 0

col_fun = colorRamp2(c(-2.5, 0, 2.5), c("blue", "white", "red"))

pdf(paste0("heatmap_CYP_NES_", my_sex, ".pdf"), height=2.5, width=4);

print(Heatmap(as.matrix(data_for_plot_NES),  na_col="lightgray", name="NES", cluster_columns=FALSE, cluster_rows=FALSE, column_title=paste0("CYP ", my_sex), col = col_fun, cell_fun = function(j, i, x, y, w, h, fill ) {
        if(data_for_plot_sig[i, j] ==1) {
                grid.text("*", x, y)
        }
}))

dev.off();
}
########


get_plot(data_out_F, "F");
get_plot(data_out_M, "M");

