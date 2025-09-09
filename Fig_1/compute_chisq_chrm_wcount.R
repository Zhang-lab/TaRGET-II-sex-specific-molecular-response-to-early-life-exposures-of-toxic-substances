.libPaths("/home/shuhua/R/x86_64-pc-linux-gnu-library/4.0/")

library(data.table)
suppressMessages(library(foreach) )
suppressMessages(library(doParallel) )

numCores<- 20;
print(paste("nCPUs:",numCores, sep=" "));
registerDoParallel(numCores)

args=commandArgs(T)

expo_count_file<- args[1];
expo_ratio_file<- args[2];
expo_meth_file<- args[3];
expo_unmeth_file<- args[4];

ctrl_ratio_file<- args[5];
ctrl_meth_file<- args[6];
ctrl_unmeth_file<- args[7];

sex<- args[8]
expo<- args[9]

chrm<- args[10]

default_settting<- function(){

expo_count_file<- "count_liver_TCDD_F"
expo_ratio_file<- "ratio_liver_TCDD_F"
expo_meth_file<- "meth_liver_TCDD_F"
expo_unmeth_file<- "unmeth_liver_TCDD_F"

ctrl_ratio_file<- "/BRC/shuhua/target/WGBS/reseq_new_analysis/processed_reseq_new/DMR_analysis_metilene/chisq_test/process_ctrl/ratio_ctrl_liver_F";
ctrl_meth_file<- "/BRC/shuhua/target/WGBS/reseq_new_analysis/processed_reseq_new/DMR_analysis_metilene/chisq_test/process_ctrl/meth_ctrl_liver_F";
ctrl_unmeth_file<- "/BRC/shuhua/target/WGBS/reseq_new_analysis/processed_reseq_new/DMR_analysis_metilene/chisq_test/process_ctrl/unmeth_ctrl_liver_F";

sex<- "F";
expo<- "TCDD";
chrm<- "chr1"
}


expo_count_data<- data.frame(fread(expo_count_file,header=T), row.names=1);
expo_ratio_data<- data.frame(fread(expo_ratio_file,header=T), row.names=1);
expo_meth_data<- data.frame(fread(expo_meth_file,header=T), row.names=1);
expo_unmeth_data<- data.frame(fread(expo_unmeth_file,header=T), row.names=1);

ctrl_ratio_data<- data.frame(fread(ctrl_ratio_file,header=T), row.names=1);
ctrl_meth_data<- data.frame(fread(ctrl_meth_file,header=T), row.names=1);
ctrl_unmeth_data<- data.frame(fread(ctrl_unmeth_file,header=T), row.names=1);


shared_wins<- rownames(expo_ratio_data)[rownames(expo_ratio_data) %in% rownames(ctrl_ratio_data)]
#For 1 chrm
shared_wins<- shared_wins[grep(paste0(chrm, "_"), shared_wins)]

print(paste0("Shard wins: ", length(shared_wins)));


expo_ratio_data1<- expo_ratio_data[shared_wins,]
expo_meth_data1<- expo_meth_data[shared_wins,]
expo_unmeth_data1<- expo_unmeth_data[shared_wins,]

ctrl_ratio_data1<- ctrl_ratio_data[shared_wins,]
ctrl_meth_data1<- ctrl_meth_data[shared_wins,]
ctrl_unmeth_data1<- ctrl_unmeth_data[shared_wins,]

########################################
compute_chisq<- function(my_i){

C_Me_num<- 0;
C_Un_num<- 0;
C_meth_vals<- c();

w1<- ctrl_ratio_data1[my_i,]
w2<- ctrl_meth_data1[my_i,]
w3<- ctrl_unmeth_data1[my_i,]

if(nrow(w1)==1){
   w1<- unlist(w1);
   w2<- unlist(w2);
   w3<- unlist(w3);
   w1<- w1[!is.na(w1)]
   w2<- w2[!is.na(w2)]
   w3<- w3[!is.na(w3)]
   C_meth_vals<- as.numeric(w1);
   C_Me_num<- sum(as.numeric(w2));
   C_Un_num<- sum(as.numeric(w3));
  }

E_Me_num<- 0;
E_Un_num<- 0;
E_meth_vals<- c();

v1<- expo_ratio_data1[my_i,]
v2<- expo_meth_data1[my_i,]
v3<- expo_unmeth_data1[my_i,]

if(nrow(v1)==1){
   v1<- unlist(v1);
   v2<- unlist(v2);
   v3<- unlist(v3);
   v1<- v1[!is.na(v1)]
   v2<- v2[!is.na(v2)]
   v3<- v3[!is.na(v3)]
   E_meth_vals<- as.numeric(v1);
   E_Me_num<- sum(as.numeric(v2));
   E_Un_num<- sum(as.numeric(v3));
  }

#Add pseudo count as 1 to "0"
if(E_Un_num==0){
   E_Un_num<- 1;
  }
if(C_Un_num==0){
   C_Un_num<- 1;
  }

pval<- chisq.test(matrix(c(E_Me_num, E_Un_num, C_Me_num, C_Un_num), nrow=2))$p.value;
ratio_diff<- mean(E_meth_vals) -mean(C_meth_vals);
return(c(shared_wins[my_i], pval, ratio_diff) )
}
#####****************************** End function compute_chisq

system.time({t<- foreach(i=1:length(shared_wins), combine=rbind) %dopar% { print(i); compute_chisq(i) }  })


data_out<- data.frame(matrix(unlist(t), ncol=3, byrow=T))

cpgNum<- rowMeans(expo_count_data[match(data_out[,1], rownames(expo_count_data)),])
data_out<- cbind.data.frame(data_out, cpgNum);

colnames(data_out)<- c("pos", "pval", "delta", "cpgNum");

data_out$delta<- as.numeric(as.character(data_out$delta));
data_out$pval<- as.numeric(as.character(data_out$pval));

data_out1<- data_out[ which((!is.na(data_out$pval)) & abs(data_out$delta) >0.1),]

qval_vec<- p.adjust(data_out1$pval)

data_out2<- cbind.data.frame(data_out1, qval_vec);

write.table(data_out, file=paste0("unfiltered_region_", sex, "_", expo, "_", chrm), quote=F, sep="\t", col.names=T, row.names=F);
write.table(data_out1, file=paste0("pval_region_", sex, "_", expo, "_", chrm), quote=F, sep="\t", col.names=T, row.names=F);
write.table(data_out2, file=paste0("qval_region_", sex, "_", expo, "_", chrm), quote=F, sep="\t", col.names=T, row.names=F);




