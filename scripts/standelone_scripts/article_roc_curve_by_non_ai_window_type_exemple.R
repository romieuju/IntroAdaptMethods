###################################################ROC ALL#############################################################################
#Script used to obtain graphs representing the ROC curves and mccf1 curves for a test set (adjacente, independante neutral introgression and second chromosome windows) 
library(extraDistr, quietly=TRUE)
library(ini, quietly = TRUE)
library(mccf1, quietly = TRUE)
library(ROCR, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(spaMM, quietly = TRUE)
library(dplyr, quietly = TRUE)

#path of introadapt.R scripts
source("/home/romieuju/ABC_Introadapt/ABC_IntroAdapt/scripts/introadapt.R")
#path of comparison folder (folder created after end_evaluation rule in snakemake file)
eval_dir            <- "/home/romieuju/results/classification/human_2chro_s_inter/comparison"
setwd(eval_dir)
#path to option file of the project
options_file        <- "/home/romieuju/results/classification/human_2chro_s_inter./project_options.ini"
#path to project prediction file with all score value by windows by methods after end_evaluation rule execution
predictiontable_all <- read.csv("/home/romieuju/results/classification/human_2chro_s_inter/comparison/predictions_human_2chro_s_inter..csv", sep="")
#type of non-AI window in the test set
analyse_type        <- "CHRO2" # for second chromosome windows or "NI" (for Independant neutral introgression window) or "ADJ" (for adjavent windows)
#path to neutral (NI) project prediction file with all score value by windows by methods after end_evaluation rule execution
file_NI             <- "/home/romieuju/results/human_2chro_neutral/comparison/predictions_human_2chro_neutral.csv"

#####to calculate volcano max########
predictiontable_IN <- read.csv(file_NI, sep="")
max_vol_LR <- max(predictiontable_all[predictiontable_all$method=="vol_rel",]$predictions_no_norm, predictiontable_IN[predictiontable_IN$method=="vol_rel",]$predictions_no_norm)

if (analyse_type=="ADJ"){
  predictiontable <- predictiontable_all
  predictiontable[predictiontable$method=="vol_rel",]$predictions<-predictiontable[predictiontable$method=="vol_rel",]$predictions_no_norm/max_vol_LR
  predictiontable <- predictiontable[predictiontable$start<1000000,]
} else if (analyse_type=="NI"){
  predictiontable    <- predictiontable_all
  predictiontable[predictiontable$method=="vol_rel",]$predictions <- predictiontable[predictiontable$method=="vol_rel",]$predictions_no_norm/max_vol_LR
  predictiontable    <- predictiontable[predictiontable$classifier==1,]
  predictiontable_IN <- read.csv(file_NI, sep="")
  predictiontable_IN <- predictiontable_IN[predictiontable_IN$classifier==0,]
  predictiontable_IN[predictiontable_IN$method=="vol_rel",]$predictions<-predictiontable_IN[predictiontable_IN$method=="vol_rel",]$predictions_no_norm/max_vol_LR
  predictiontable <- rbind(predictiontable, predictiontable_IN)
}else if (analyse_type=="CHRO2"){
  predictiontable    <- predictiontable_all
  predictiontable[predictiontable$method=="vol_rel",]$predictions <- predictiontable[predictiontable$method=="vol_rel",]$predictions_no_norm/max_vol_LR
  predictiontable_2chro <- predictiontable[predictiontable$end>1000000,]
  predictiontable <- predictiontable[predictiontable$classifier==1,]
  predictiontable <- rbind(predictiontable, predictiontable_2chro)
}

opts     = get_project_options(options_file)
Settings = opts$Settings
Inference = opts$Inference
Genome = opts$Genome
analysis     <-Settings$analysis
project     <-Settings$project
project_dir <-Settings$project_dir
simul_type <-Settings$simul_type
maladapt_switch   <- Inference$maladapt_switch
volcano_switch    <- Inference$volcano_switch
genomatnn_switch  <- Inference$genomatnn_switch
volcano_threshold   <- Inference$volcano_threshold
genomatnn_threshold <- Inference$genomatnn_threshold
maladapt_threshold  <- Inference$maladapt_threshold
q95_threshold       <- Inference$q95_threshold
mal_feat        <- unlist(strsplit(Inference$maladapt_features," "))
mal_method_list <- paste("mal_",mal_feat,sep="")
gnn_mods        <- unlist(strsplit(Inference$genomatnn_trained_mod," "))
gnn_method_list <- paste("gnn_",1:length(gnn_mods),sep="")
method_list<- list()
cbPalette<-list()
cbPoint<-list()
cbLinetype<-list()
cblabels<-list()

if (genomatnn_switch == "On"){
  method_list<-c(unlist(method_list),gnn_method_list)
  if (length(gnn_method_list)==1){
    if(gnn_mods=="Nea_to_CEU_af-0.25_2918410235"){
      cbPalette<-c(unlist(cbPalette), "#E69F00")
      cbPoint<-c(unlist(cbPoint), 2)
      cbLinetype<-c(unlist(cbLinetype), "twodash")
      cblabels<-c(unlist(cblabels), "genomatnn_0.25")
    }else if(gnn_mods=="Nea_to_CEU_af-0.05_2250018620"){
      cbPalette<-c(unlist(cbPalette), "#F0E442")
      cbPoint<-c(unlist(cbPoint), 6)
      cbLinetype<-c(unlist(cbLinetype), "dotdash")
      cblabels<-c(unlist(cblabels), "genomatnn_0.05")
    }
  }else{
    cbPalette<-c(unlist(cbPalette), "#E69F00", "#F0E442")
    cbPoint<-c(unlist(cbPoint), 2, 6)
    cbLinetype<-c(unlist(cbLinetype), "twodash", "dotdash")
    cblabels<-c(unlist(cblabels), "genomatnn_0.25", "genomatnn_0.05")
  }
}
if (maladapt_switch == "On"){
  method_list<-c(unlist(method_list),mal_method_list)
  if (length(mal_method_list)==1){
    if(mal_method_list=="mal_noenoZnS"){
      cbPalette<-c(unlist(cbPalette), "#56B4E9")
      cbPoint<-c(unlist(cbPoint), 0)
      cbLinetype<-c(unlist(cbLinetype), "longdash")
      cblabels<-c(unlist(cblabels), "MaLAdapt_Q95")
    }else if(mal_method_list=="mal_noQnoenoZnS"){
      cbPalette<-c(unlist(cbPalette), "#0072B2")
      cbPoint<-c(unlist(cbPoint), 5)
      cbLinetype<-c(unlist(cbLinetype), "dashed")
      cblabels<-c(unlist(cblabels), "MaLAdapt_noQ95")
    }
  }else{
    cbPalette<-c(unlist(cbPalette), "#56B4E9", "#0072B2")
    cbPoint<-c(unlist(cbPoint), 0, 5)
    cbLinetype<-c(unlist(cbLinetype), "longdash", "dashed")
    cblabels<-c(unlist(cblabels), "MaLAdapt_Q95", "MaLAdapt_noQ95")
  }
}
method_list<-c(unlist(method_list),"Q_1_100_q95")
cbPalette<-c(unlist(cbPalette), "#009E73")
cbPoint<-c(unlist(cbPoint), 1)
cbLinetype<-c(unlist(cbLinetype), "3D")
cblabels<-c(unlist(cblabels),"Q95(1%,100%)")
if (volcano_switch == "On"){
  method_list<-c(unlist(method_list),"vol_rel")
  cbPalette<-c(unlist(cbPalette), "#D55E00")
  cbPoint<-c(unlist(cbPoint), 8)
  cbLinetype<-c(unlist(cbLinetype), "solid")
  cblabels<-c(unlist(cblabels),"VolcanoFinder")
}
perfectpalette <- c("black","#0072B2","#56B4E9","#00A300","#228B22","#00CC00","#F0E442","#E69F00", "#D55E00")
currentthreshold <- NaN
predictiontable_end <- data.frame()
graphtable<-data.frame()
summarytable<-data.frame()
######Neutral################
for(method in method_list){
  currentprediction<-predictiontable[predictiontable$method==method,]
  pred<-prediction(currentprediction$predictions, currentprediction$classifier)
  currentsummary<-data.frame(method=method,apriori_threshold=currentthreshold)
  if(is.na(currentthreshold)){
    actual_threshold<-NA
    currentsummary$tp<-NA
    currentsummary$tn<-NA
    currentsummary$fp<-NA
    currentsummary$fn<-NA
    currentsummary$accuracy<-NA
    currentsummary$sensitivity<-NA
    currentsummary$specificity<-NA
    currentsummary$precision<-NA
    currentsummary$npv<-NA
    currentprediction$bin_pred<-NA
  }else{
    currentprediction<-currentprediction%>%
      mutate(bin_pred=as.numeric(predictions>=currentthreshold))
    #confusion matrix info
    currentsummary$tp<-sum(currentprediction$bin_pred & currentprediction$classifier)
    currentsummary$tn<-sum(!currentprediction$bin_pred & !currentprediction$classifier)
    currentsummary$fp<-sum(currentprediction$bin_pred & !currentprediction$classifier)
    currentsummary$fn<-sum(!currentprediction$bin_pred & currentprediction$classifier)
    #metrics
    currentsummary$accuracy<-(currentsummary$tp+currentsummary$tn)/(currentsummary$tp+currentsummary$tn+currentsummary$fp+currentsummary$fn)
    currentsummary$sensitivity<-(currentsummary$tp)/(currentsummary$tp+currentsummary$fn)
    currentsummary$specificity<-(currentsummary$tn)/(currentsummary$tn+currentsummary$fp)
    currentsummary$precision<-(currentsummary$tp)/(currentsummary$tp+currentsummary$fp)
    currentsummary$npv<-(currentsummary$tn)/(currentsummary$tn+currentsummary$fn)
    #takes the lowest value among the data above the threshold (give "Inf" if there are none, which should also be recognised)
    actual_threshold<-min(currentprediction$predictions[as.logical(currentprediction$bin_pred)])
  }
  currentsummary$auroc<-unlist(performance(pred, "auc")@y.values)#area under ROC curve
  currentsummary$aupr<-unlist(performance(pred, "aucpr")@y.values)#area under PR curve
  #gives coordinates for the mcc-f1 curve
  mccf1out<-mccf1(response=currentprediction$classifier,currentprediction$predictions)
  #gives the mccf1 metric and the "best threshold" to use, and saves them in the dataframe
  mccf1sum<-summary(mccf1out)
  best_threshold<-mccf1sum$best_threshold
  currentsummary$best_threshold<-best_threshold
  currentsummary$mccf1_metric<-mccf1sum$mccf1_metric
  #adds mcc-f1 curve coordinates to the relevant table
  currentgraph<-data.frame(method=method,ideal=(best_threshold==mccf1out$thresholds),focus=(actual_threshold==mccf1out$thresholds),graph="mccf1",x=mccf1out$f1,y=mccf1out$normalized_mcc,thresholds=mccf1out$thresholds)
  graphtable<-rbind(graphtable,currentgraph)
  #roc
  roc<-performance(pred,"tpr","fpr")
  currentgraph<-data.frame(method=method,ideal=(best_threshold==unlist(roc@alpha.values)),focus=(actual_threshold==unlist(roc@alpha.values)),graph="roc",x=unlist(roc@x.values),y=unlist(roc@y.values),thresholds=unlist(roc@alpha.values))
  graphtable<-rbind(graphtable,currentgraph)
  currentprediction<-currentprediction%>%
    mutate(best_bin_pred=as.numeric(predictions>=best_threshold))
  #confusion matrix info
  currentsummary$best_tp<-unlist(pred@tp)[unlist(pred@cutoffs)==best_threshold]#currentsummary$best_tp<-sum(currentprediction$best_bin_pred & currentprediction$classifier)
  currentsummary$best_tn<-unlist(pred@tn)[unlist(pred@cutoffs)==best_threshold]#currentsummary$best_tn<-sum(!currentprediction$best_bin_pred & !currentprediction$classifier)
  currentsummary$best_fp<-unlist(pred@fp)[unlist(pred@cutoffs)==best_threshold]#currentsummary$best_fp<-sum(currentprediction$best_bin_pred & !currentprediction$classifier)
  currentsummary$best_fn<-unlist(pred@fn)[unlist(pred@cutoffs)==best_threshold]#currentsummary$best_fn<-sum(!currentprediction$best_bin_pred & currentprediction$classifier)
  #metrics
  currentsummary$best_accuracy<-(currentsummary$best_tp+currentsummary$best_tn)/(currentsummary$best_tp+currentsummary$best_tn+currentsummary$best_fp+currentsummary$best_fn)
  currentsummary$best_sensitivity<-(currentsummary$best_tp)/(currentsummary$best_tp+currentsummary$best_fn)
  currentsummary$best_specificity<-(currentsummary$best_tn)/(currentsummary$best_tn+currentsummary$best_fp)
  currentsummary$best_precision<-(currentsummary$best_tp)/(currentsummary$best_tp+currentsummary$best_fp)
  currentsummary$best_npv<-(currentsummary$best_tn)/(currentsummary$best_tn+currentsummary$best_fn)
  fpr<-performance( pred, "fpr" )
  #if(method == "Q_1_100_q95"){
  #  best_q95_fpr<-unlist(fpr@y.values)[unlist(fpr@x.values)==best_threshold]
  #}
  #currentsummary$q95like_threshold<-min(unlist(fpr@x.values)[unlist(fpr@y.values)<=best_q95_fpr])
  #currentsummary$fpr0.05_threshold<-min(unlist(fpr@x.values)[unlist(fpr@y.values)<=best_threshold])
  currentsummary$fpr0.05_threshold<-min(unlist(fpr@x.values)[unlist(fpr@y.values)<=0.05])
  #currentprediction<-currentprediction%>%
  #  mutate(q95like_bin_pred=as.numeric(predictions>=currentsummary$q95like_threshold))
  currentprediction<-currentprediction%>%
    mutate(fpr0.05_bin_pred=as.numeric(predictions>=currentsummary$fpr0.05_threshold))
  #confusion matrix info
  currentsummary$fpr0.05_tp<-unlist(pred@tp)[unlist(pred@cutoffs)==currentsummary$fpr0.05_threshold]#currentsummary$best_tp<-sum(currentprediction$best_bin_pred & currentprediction$classifier)
  currentsummary$fpr0.05_tn<-unlist(pred@tn)[unlist(pred@cutoffs)==currentsummary$fpr0.05_threshold]#currentsummary$best_tn<-sum(!currentprediction$best_bin_pred & !currentprediction$classifier)
  currentsummary$fpr0.05_fp<-unlist(pred@fp)[unlist(pred@cutoffs)==currentsummary$fpr0.05_threshold]#currentsummary$best_fp<-sum(currentprediction$best_bin_pred & !currentprediction$classifier)
  currentsummary$fpr0.05_fn<-unlist(pred@fn)[unlist(pred@cutoffs)==currentsummary$fpr0.05_threshold]#currentsummary$best_fn<-sum(!currentprediction$best_bin_pred & currentprediction$classifier)
  #metrics
  currentsummary$fpr0.05_accuracy<-(currentsummary$fpr0.05_tp+currentsummary$fpr0.05_tn)/(currentsummary$fpr0.05_tp+currentsummary$fpr0.05_tn+currentsummary$fpr0.05_fp+currentsummary$fpr0.05_fn)
  currentsummary$fpr0.05_sensitivity<-(currentsummary$fpr0.05_tp)/(currentsummary$fpr0.05_tp+currentsummary$fpr0.05_fn)
  currentsummary$fpr0.05_specificity<-(currentsummary$fpr0.05_tn)/(currentsummary$fpr0.05_tn+currentsummary$fpr0.05_fp)
  currentsummary$fpr0.05_precision<-(currentsummary$fpr0.05_tp)/(currentsummary$fpr0.05_tp+currentsummary$fpr0.05_fp)
  currentsummary$fpr0.05_npv<-(currentsummary$fpr0.05_tn)/(currentsummary$fpr0.05_tn+currentsummary$fpr0.05_fn)
  summarytable<-rbind(summarytable,currentsummary)
  predictiontable_end<-rbind(predictiontable_end,currentprediction)
}
cbLinetype<-c(88,43,42,61,111,95)
cbPalette <- c("#D55E00", "#F0E442","#0072B2","#56B4E9","#E69F00","black")
cbPoint<-c(88,43,42,61,111,95)

graphtable$method<-factor(graphtable$method, levels=c("vol_rel", "Q_1_100_q95", "mal_noenoZnS","mal_noQnoenoZnS", "gnn_1","gnn_2"))
cbPalette <- c("#D55E00" ,"#009E73", "#56B4E9","#0072B2", "#E69F00", "#F0E442")
cblabels<- c("VolcanoFinder","Q95(1%,100%)", "MaLAdapt_Q95", "MaLAdapt_noQ95", "genomatnn_0.25", "genomatnn_0.05")
cbLinetype<-c("solid","twodash","dotdash","longdash","dashed","3D")
graphtablebis <-graphtable #graphtable[graphtable$method=="mal_noQnoenoZnS",]
print(ggplot(graphtablebis[graphtablebis$graph=="roc" & graphtable$x!=0 & graphtable$y!=0 & is.finite(graphtable$thresholds),],aes(x=x,y=y,group=method,shape=method))+
        geom_line(aes(color=method, linetype=method),  linewidth = 1.0, alpha=1.0)+
        scale_shape_manual(values = cbPoint, labels=cblabels)+
        scale_colour_manual(values=cbPalette, labels=cblabels)+
        scale_linetype_manual(values = cbLinetype, labels = cblabels)+
        xlim(0,1)+
        ylim(0,1)+
        geom_abline(aes(intercept=0,slope=1), linetype="dotted")+
        labs(x="False positive rate",y="True positive rate",title="")+
        theme(panel.background = element_rect(fill = "white"), panel.grid = element_line(color = "#B5DFED",
                                                                                         linewidth = 0.2,
                                                                                         linetype = 2))+
        theme(text=element_text(size=20))+
        theme(plot.title = element_text(hjust = 0, vjust = -4.5)))

rocplot <- ggplot(graphtablebis[graphtablebis$graph=="roc" & graphtable$x!=0 & graphtable$y!=0 & is.finite(graphtable$thresholds),],aes(x=x,y=y,group=method,shape=method))+
  geom_line(aes(color=method, linetype=method),  linewidth = 1.3, alpha=1.0)+
  scale_shape_manual(values = cbPoint, labels=cblabels)+
  scale_colour_manual(values=cbPalette, labels=cblabels)+
  scale_linetype_manual(values = cbLinetype, labels = cblabels)+
  xlim(0,1)+
  ylim(0,1)+
  geom_abline(aes(intercept=0,slope=1), linetype="dotted")+
  labs(x="False positive rate",y="True positive rate",title="")+
  theme(panel.background = element_rect(fill = "white"), panel.grid = element_line(color = "#B5DFED",
                                                                                   linewidth = 0.2,
                                                                                   linetype = 2))+
  theme(text=element_text(size=20))+
  theme(plot.title = element_text(hjust = 0, vjust = -4.5))
ggsave(paste("ROC_AI_VS_",analyse_type,"_all_",project,".pdf",sep=""), rocplot, width = 11.69, height = 8.27, units = "in", dpi = 300, device = "pdf")


cbPointmccf1<-c(8,0,5,2,6,4)
graphtable$method<-factor(graphtable$method, levels=c("vol_rel", "mal_noenoZnS","mal_noQnoenoZnS", "gnn_1","gnn_2", "Q_1_100_q95"))
cbPalette <- c("#D55E00" , "#56B4E9","#0072B2", "#E69F00", "#F0E442", "#009E73")
cblabels<- c("VolcanoFinder", "MaLAdapt_Q95", "MaLAdapt_noQ95", "genomatnn_0.25", "genomatnn_0.05", "Q95(1%,100%)")
mccf1plot<-ggplot(graphtable[graphtable$graph=="mccf1",],aes(x=x,y=y,group=method,color=method,shape=method, alpha = 0.7))+
  geom_path()+
  geom_point(data = subset(graphtable[graphtable$graph=="mccf1",], ideal), size=3, alpha = 1.0)+
  geom_hline(aes(yintercept=0.5), linetype="dotted")+
  xlim(0,1)+
  ylim(0,1)+
  labs(x="F1",y="Normalized MCC",title="")+
  scale_shape_manual(values = cbPointmccf1, labels=cblabels)+
  scale_colour_manual(values=cbPalette, labels=cblabels)+
  scale_linetype_manual(values = cbLinetype, labels = cblabels)+
  theme(panel.background = element_rect(fill = "white"), panel.grid = element_line(color = "#B5DFED",
                                                                                   linewidth = 0.2,
                                                                                   linetype = 2))+
  theme(text=element_text(size=20))+
  theme(plot.title = element_text(hjust = 0, vjust = -4.5))
ggsave(paste("MCCF1_AI_VS_",analyse_type,"_all_",project,".pdf",sep=""), mccf1plot, width = 11.69, height = 8.27, units = "in", dpi = 300, device = "pdf")

#Save AUROC results and confusion matrix + mccf1 score (summarytable file)
auroc_all<-data.frame("method"=c("gnn_1","gnn_2","mal_noenoZnS","mal_noQnoenoZnS","Q_1_100_q95","vol_rel"))
auroc_table<-summarytable[,c("method", "auroc")]
auroc_all <- merge(auroc_all,auroc_table,by="method")
write.table(auroc_table, file =paste("auroc_AI_VS_",analyse_type,"_all_",project,".csv",sep=""), col.names = TRUE, row.names = FALSE)
write.table(summarytable, file =paste("summarytable_AI_VS_",analyse_type,"_all_",project,".csv",sep=""), col.names = TRUE, row.names = FALSE)
