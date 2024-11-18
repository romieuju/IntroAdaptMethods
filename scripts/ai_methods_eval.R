#    IntroAdapt: joint inference of introgression and selection
#    Copyright (C) 2024 Ghislain Camarata and Jules Romieu, ISEM/CBGP, CNRS/INRAE. 
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.

# This script is used to calculate performance statistics and create associated figures
# for Racimo's Q95, MaLAdapt, genomatnn and VolcanoFinder

# R package used for this script
library(extraDistr, quietly=TRUE)
library(ini, quietly = TRUE)
source("scripts/timeadapt.R")
library(mccf1, quietly = TRUE)
library(ROCR, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(spaMM, quietly = TRUE)

#timer
compare_start <- Sys.time()

# read command line arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0){
  stop("Two positional arguments necessary in the command (RDS options file and batch number)")
  quit(save="no")
} else {
  options_file = args[1]
}

#Save the ini information on variables :
opts      <- get_project_options(options_file)
Settings  <- opts$Settings
Inference <- opts$Inference
Genome    <- opts$Genome

#where are your results ?
analysis    <- Settings$analysis
project     <- Settings$project
project_dir <- Settings$project_dir
simul_type  <- Settings$simul_type
eval_dir    <- paste(project_dir, "comparison",sep="/")
dir.create(eval_dir, showWarnings = FALSE)

#what are you comparing ?
maladapt_switch   <- Inference$maladapt_switch
volcano_switch    <- Inference$volcano_switch
genomatnn_switch  <- Inference$genomatnn_switch

#which values of the criteria are sufficient to be classified as AI ?
volcano_threshold   <- Inference$volcano_threshold
genomatnn_threshold <- Inference$genomatnn_threshold
maladapt_threshold  <- Inference$maladapt_threshold
q95_threshold       <- Inference$q95_threshold

#which models do you use ?
mal_feat        <- unlist(strsplit(Inference$maladapt_features," "))
mal_method_list <- paste("mal_",mal_feat,sep="")
gnn_mods        <- unlist(strsplit(Inference$genomatnn_trained_mod," "))
gnn_method_list <- paste("gnn_",1:length(gnn_mods),sep="")

#loads table containing maladapt outputs (conditionnal) and inputs (mandatory)
if (maladapt_switch == "On"){
  basetable<-read.csv(paste(project_dir,'/maladapt_out_',project,'.csv',sep=""))
}else{
  basetable<-read.csv(paste(project_dir,'/reftable_wind_',analysis,"_",project,'.csv',sep=""),sep=" ")  %>%
  rename_at('num_sim', ~'sim_id') %>%
  rename_at('IA', ~'classifier')
}

if(simul_type!="AI" & max(basetable$classifier)==1){
  stop("Windows tagged as AI in a non AI project")
}

if(simul_type=="AI" & max(basetable$classifier)==0){
  stop("No AI windows in an AI project")
}

if((simul_type=="NI") & max(basetable$mean_p1)==0){
  stop("No introgression in a neutral introgression project")
}

#filters out the simulations with no introgression if there should be neutral introgression
if(simul_type=="NI"){
  basetable<-basetable  %>%
  filter(mean_p1>0)
}

#filters out the simulations with no AI if there should be some
if(simul_type=="AI"){
  basetable<-basetable %>%
  group_by(sim_id) %>%
  mutate(any_AI=max(classifier)) %>%
  ungroup() %>%
  filter(any_AI==1) %>%
  data.frame()
}

simul_vector<-unique(basetable$sim_id)

gnntotal<-data.frame()
#For each simulation that passes the filters, load the relevant inference : 
for (sim in simul_vector){
  # find simulation path
  cat(project,sim,"\n")
  job_dir     <-paste0(analysis,"_",project,"_sim_", sim, sep="")
  job_dir_path<-paste(project_dir, job_dir,sep="/")
  if (Settings$verbose >= 10) write(paste("Job folder:", job_dir), stdout())

  #adds the beneficial allele frequency to the table
  basetable$freq_benef[basetable$sim_id==sim] <- 0
  benef_freq                                  <- read.csv(paste(job_dir_path,'/forwsim_freq_mut_don_in_rec_',analysis,'_',project,'_sim_',sim,'.txt',sep=""),sep=",")
  if (is.nan(benef_freq$position[1])==FALSE){
    for(mut in 1:nrow(benef_freq)){
      cur_pos  <- benef_freq$position[mut]
      cur_freq <- benef_freq$frequency[mut]
      basetable$freq_benef[(basetable$start<=cur_pos & basetable$end>=cur_pos)&basetable$sim_id==sim] <- max(cur_freq,basetable$freq_benef[(basetable$start>=cur_pos & basetable$end<=cur_pos)&basetable$sim_id==sim])
    }
  }

  if (volcano_switch == "On"){
    #groups volcanofinder's output for every chromosome of the simulation
    volfull <- data.frame()
    for (chr in seq_len(Genome$nchr)) {
      volout <- read.csv(paste(job_dir_path,'/VolcanoFinder/volcanotest_',sim,'_chr_',chr,'.out',sep=""),sep="\t")
      volfull <- rbind(volfull,volout)
    }
    
    max_scores    <- c()
    max_positions <- c()
    
    # for each window, finds the highest volcanofinder score and the matching position
    for (i in 1:sum(basetable$sim_id==sim)) {
      positions     <- volfull$location[volfull$location > basetable$start[basetable$sim_id==sim][i] & volfull$location <= basetable$end[basetable$sim_id==sim][i]]
      max_scores    <- c(max_scores,max(volfull$LR[volfull$location %in% positions]))
      max_positions <- c(max_positions,positions[which.max(volfull$LR[volfull$location %in% positions])])
    }
    
    # adds the results to the maladapt table
    basetable$vol_loc[basetable$sim_id==sim] <- max_positions
    basetable$vol_LR[basetable$sim_id==sim]  <- max_scores
  }
  
  if (genomatnn_switch == "On"){
    gnntable <- data.frame()
    for (mod_id in 1:length(gnn_mods)){
      #loads genomatnn outputs and determines wether there was AI or not using the maladapt table
      mod_table        <- read.csv(paste(job_dir_path,'/genomatnn/',gnn_mods[mod_id],'/predictions.txt',sep=""),sep="\t")
      mod_table$method <- gnn_method_list[mod_id]
      gnntable         <- rbind(gnntable,mod_table)
    }
    for (i in 1:nrow(gnntable)) {
      if(sum(basetable$classifier[basetable$sim_id==sim & (basetable$start==(gnntable$start[i]-1)|basetable$end==(gnntable$end[i]-1))])==0){
        gnntable$classifier[i] <- 0
        gnntable$freq_benef[i] <- 0
      }else{
        gnntable$classifier[i] <- 1
        gnntable$freq_benef[i] <- max(basetable$freq_benef[basetable$sim_id==sim & (basetable$start==(gnntable$start[i]-1)|basetable$end==(gnntable$end[i]-1))])
      }
    }
    gnntable$sim_id            <- sim
    gnntable$fixation_g_r_back <- unique(basetable$fixation_g_r_back[basetable$sim_id==sim])
    gnntable$mean_p1           <- unique(basetable$mean_p1[basetable$sim_id==sim])
    gnntotal                   <- rbind(gnntotal,gnntable)
  }
}

if (volcano_switch == "On"){
#gives the relative distribution of VolcanoFinder, to make comparizons with other methods easier
basetable<-basetable %>%
  group_by(sim_id) %>%
  mutate(vol_rel= vol_LR / max(vol_LR)) %>%
  ungroup() %>%
  data.frame()

}

predictiontable <- data.frame()
summarytable    <- data.frame()
graphtable      <- data.frame()

method_list<-c("Q_1_100_q95")
if (volcano_switch == "On"){method_list<-c(method_list,"vol_rel")}
if (genomatnn_switch == "On"){method_list<-c(method_list,gnn_method_list)}
if (maladapt_switch == "On"){method_list<-c(method_list,mal_method_list)}

#for each method to compare 
for(method in method_list){
  if(method %in% gnn_method_list){
    currentprediction <- data.frame(sim_id=gnntotal$sim_id[gnntotal$method==method],start=gnntotal$start[gnntotal$method==method],end=gnntotal$end[gnntotal$method==method],method=method,classifier=gnntotal$classifier[gnntotal$method==method],predictions=gnntotal$Pr.AI.[gnntotal$method==method],freq_benef=gnntotal$freq_benef[gnntotal$method==method],fixation_g_r_back=gnntotal$fixation_g_r_back[gnntotal$method==method],mean_p1=gnntotal$mean_p1[gnntotal$method==method])
  }else{
    currentprediction <- data.frame(sim_id=basetable$sim_id,start=basetable$start,end=basetable$end,method=method,classifier=basetable$classifier,predictions=basetable[[method]],freq_benef=basetable$freq_benef,fixation_g_r_back=basetable$fixation_g_r_back,mean_p1=basetable$mean_p1)
  }

  #uses the relevant threshold for this method
  if(method=="Q_1_100_q95"){
    currentthreshold <- q95_threshold
  }else if(method=="vol_rel"){
    currentthreshold <- volcano_threshold
  }else if(method %in% gnn_method_list){
    currentthreshold <- genomatnn_threshold
  }else if(method %in% mal_method_list){
    currentthreshold <- maladapt_threshold
  }

  #creates the dataframe that will hold all the metrics for this method
  currentsummary <- data.frame(method=method,apriori_threshold=currentthreshold)

  #uses the threshold set a priori if there is one in order to create the confusion matrix and the corresponding metrics
  if(is.na(currentthreshold)){
    actual_threshold           <- NA
    currentsummary$tp          <- NA
    currentsummary$tn          <- NA
    currentsummary$fp          <- NA
    currentsummary$fn          <- NA
    currentsummary$accuracy    <- NA
    currentsummary$sensitivity <- NA
    currentsummary$specificity <- NA
    currentsummary$precision   <- NA
    currentsummary$npv         <- NA
    currentprediction$bin_pred <- NA
  }else{
    currentprediction<-currentprediction%>%
      mutate(bin_pred=as.numeric(predictions>=currentthreshold))
      
    #confusion matrix info
    currentsummary$tp <- sum(currentprediction$bin_pred & currentprediction$classifier)
    currentsummary$tn <- sum(!currentprediction$bin_pred & !currentprediction$classifier)
    currentsummary$fp <- sum(currentprediction$bin_pred & !currentprediction$classifier)
    currentsummary$fn <- sum(!currentprediction$bin_pred & currentprediction$classifier)

    #metrics
    currentsummary$accuracy    <- (currentsummary$tp+currentsummary$tn)/(currentsummary$tp+currentsummary$tn+currentsummary$fp+currentsummary$fn)
    currentsummary$sensitivity <- (currentsummary$tp)/(currentsummary$tp+currentsummary$fn)
    currentsummary$specificity <- (currentsummary$tn)/(currentsummary$tn+currentsummary$fp)
    currentsummary$precision   <- (currentsummary$tp)/(currentsummary$tp+currentsummary$fp)
    currentsummary$npv         <- (currentsummary$tn)/(currentsummary$tn+currentsummary$fn)

    #takes the lowest value among the data above the threshold (give "Inf" if there are none, which should also be recognised)
    actual_threshold           <- min(currentprediction$predictions[as.logical(currentprediction$bin_pred)])
  }

  if(length(unique(currentprediction$classifier))==2 & length(unique(currentprediction$predictions))>1){  
    #puts the data in a ROCR friendly format
    pred                          <-prediction(currentprediction$predictions, currentprediction$classifier)
    currentsummary$auroc          <- unlist(performance(pred, "auc")@y.values)#area under ROC curve
    currentsummary$aupr           <- unlist(performance(pred, "aucpr")@y.values)#area under PR curve
  
    #gives coordinates for the mcc-f1 curve 
    mccf1out                      <- mccf1(response=currentprediction$classifier,currentprediction$predictions)

    #gives the mccf1 metric and the "best threshold" to use, and saves them in the dataframe
    mccf1sum                      <- summary(mccf1out)
    best_threshold                <- mccf1sum$best_threshold
    currentsummary$best_threshold <- best_threshold
    currentsummary$mccf1_metric   <- mccf1sum$mccf1_metric

    #adds mcc-f1 curve coordinates to the relevant table
    currentgraph                  <- data.frame(method=method,ideal=(best_threshold==mccf1out$thresholds),focus=(actual_threshold==mccf1out$thresholds),graph="mccf1",x=mccf1out$f1,y=mccf1out$normalized_mcc,thresholds=mccf1out$thresholds)
    graphtable                    <- rbind(graphtable,currentgraph)

    #roc
    roc                           <- performance(pred,"tpr","fpr")
    currentgraph                  <- data.frame(method=method,ideal=(best_threshold==unlist(roc@alpha.values)),focus=(actual_threshold==unlist(roc@alpha.values)),graph="roc",x=unlist(roc@x.values),y=unlist(roc@y.values),thresholds=unlist(roc@alpha.values))
    graphtable                    <- rbind(graphtable,currentgraph)

    #precision/recall curve
    precrec                       <- performance(pred, "prec", "rec")
    currentgraph                  <- data.frame(method=method,ideal=(best_threshold==unlist(precrec@alpha.values)),focus=(actual_threshold==unlist(precrec@alpha.values)),graph="precrec",x=unlist(precrec@x.values),y=unlist(precrec@y.values),thresholds=unlist(precrec@alpha.values))
    graphtable                    <- rbind(graphtable,currentgraph)

    #sensitivity/specificity curve 
    senspe                        <- performance(pred, "sens", "spec")
    currentgraph                  <- data.frame(method=method,ideal=(best_threshold==unlist(senspe@alpha.values)),focus=(actual_threshold==unlist(senspe@alpha.values)),graph="senspe",x=unlist(senspe@x.values),y=unlist(senspe@y.values),thresholds=unlist(senspe@alpha.values))
    graphtable                    <- rbind(graphtable,currentgraph)

    currentprediction<-currentprediction%>%
      mutate(best_bin_pred=as.numeric(predictions>=best_threshold))
    
    #confusion matrix info
    currentsummary$best_tp <- unlist(pred@tp)[unlist(pred@cutoffs)==best_threshold]#currentsummary$best_tp<-sum(currentprediction$best_bin_pred & currentprediction$classifier)
    currentsummary$best_tn <- unlist(pred@tn)[unlist(pred@cutoffs)==best_threshold]#currentsummary$best_tn<-sum(!currentprediction$best_bin_pred & !currentprediction$classifier)
    currentsummary$best_fp <- unlist(pred@fp)[unlist(pred@cutoffs)==best_threshold]#currentsummary$best_fp<-sum(currentprediction$best_bin_pred & !currentprediction$classifier)
    currentsummary$best_fn <- unlist(pred@fn)[unlist(pred@cutoffs)==best_threshold]#currentsummary$best_fn<-sum(!currentprediction$best_bin_pred & currentprediction$classifier)
        
    #metrics
    currentsummary$best_accuracy<-(currentsummary$best_tp+currentsummary$best_tn)/(currentsummary$best_tp+currentsummary$best_tn+currentsummary$best_fp+currentsummary$best_fn)
    currentsummary$best_sensitivity<-(currentsummary$best_tp)/(currentsummary$best_tp+currentsummary$best_fn)
    currentsummary$best_specificity<-(currentsummary$best_tn)/(currentsummary$best_tn+currentsummary$best_fp)
    currentsummary$best_precision<-(currentsummary$best_tp)/(currentsummary$best_tp+currentsummary$best_fp)
    currentsummary$best_npv<-(currentsummary$best_tn)/(currentsummary$best_tn+currentsummary$best_fn)

    #how do you get a 5% error rate and what does everything else look like in these conditions ?
    fpr<-performance( pred, "fpr" )
    #if(method == "Q_1_100_q95"){
    #  best_q95_fpr<-unlist(fpr@y.values)[unlist(fpr@x.values)==best_threshold]
    #}

    #currentsummary$q95like_threshold<-min(unlist(fpr@x.values)[unlist(fpr@y.values)<=best_q95_fpr])
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
  }else{
    currentsummary$auroc<-NA
    currentsummary$aupr<-NA
    best_threshold<-NA
    currentsummary$best_threshold<-best_threshold
    currentsummary$mccf1_metric<-NA
    currentsummary$best_tp<-NA
    currentsummary$best_tn<-NA
    currentsummary$best_fp<-NA
    currentsummary$best_fn<-NA
    currentsummary$best_accuracy<-NA
    currentsummary$best_sensitivity<-NA
    currentsummary$best_specificity<-NA
    currentsummary$best_precision<-NA
    currentsummary$best_npv<-NA
    currentprediction$best_bin_pred<-NA
    #currentsummary$q95like_threshold<-NA
    currentsummary$fpr0.05_threshold<-NA
    currentsummary$fpr0.05_tp<-NA
    currentsummary$fpr0.05_tn<-NA
    currentsummary$fpr0.05_fp<-NA
    currentsummary$fpr0.05_fn<-NA
    currentsummary$fpr0.05_accuracy<-NA
    currentsummary$fpr0.05_sensitivity<-NA
    currentsummary$fpr0.05_specificity<-NA
    currentsummary$fpr0.05_precision<-NA
    currentsummary$fpr0.05_npv<-NA
    #currentprediction$q95like_bin_pred<-NA
    currentprediction$fpr0.05_bin_pred<-NA
  }

  #adds the results in the same data.frame, following the same format
  predictiontable<-rbind(predictiontable,currentprediction)

  summarytable<-rbind(summarytable,currentsummary)
}

perfectpalette <- c("black","#0072B2","#56B4E9","#00a300","#228b22","#00cc00","#F0E442","#E69F00", "#D55E00")
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#F0E442", "#D55E00", "#CC79A7")

#correlations ####
if (volcano_switch == "On"){
  #correlations volcanofinder/Q95
  newtab<-data.frame(
    x=predictiontable$predictions[predictiontable$method=="Q_1_100_q95"],
    y=predictiontable$predictions[predictiontable$method=="vol_rel"],
    frequence=predictiontable$freq_benef[predictiontable$method=="Q_1_100_q95"])
  png(file.path(eval_dir,paste("corr_vol_Q95_",project,".png",sep="")),width = 800, height = 600)
  print(ggplot(newtab,aes(x=x,y=y,color=frequence))+
    geom_point()+
    labs(x="Q95",y="VolcanoFinder")+
    theme(panel.grid = element_line(color = "#b5dfed",
                                    linewidth = 0.5,
                                    linetype = 2))+
    theme(text=element_text(size=20))+
    scale_color_gradientn(colours =perfectpalette))
  dev.off()
}

if (maladapt_switch == "On"){
  #correlations MaLAdapt/Q95
  newtab<-data.frame(
    x=rep(predictiontable$predictions[predictiontable$method=="Q_1_100_q95"],length(mal_method_list)),
    y=predictiontable$predictions[predictiontable$method %in% mal_method_list],
    frequence=predictiontable$freq_benef[predictiontable$method %in% mal_method_list],
    features=predictiontable$method[predictiontable$method %in% mal_method_list])
  png(file.path(eval_dir,paste("corr_mal_Q95_",project,".png",sep="")),width = 900, height = 600)
  print(ggplot(newtab,aes(x=x,y=y))+
    geom_point(aes(shape=features,color=frequence))+
    geom_smooth(method = lm, aes(linetype=features,group=features))+
    labs(x="Q95",y="MaLAdapt")+
    theme(panel.grid = element_line(color = "#b5dfed",
                                    linewidth = 0.5,
                                    linetype = 2))+
    theme(text=element_text(size=20))+
    scale_color_gradientn(colours =perfectpalette))
  dev.off()
}

if (genomatnn_switch == "On"){
  #correlations genomatnn/Q95
  newtab<-data.frame()
  for (model in gnn_method_list){
    mod_table<-data.frame(
      x=c(predictiontable$predictions[predictiontable$method=="Q_1_100_q95"][seq(1, sum(predictiontable$method=="Q_1_100_q95"), 2)],predictiontable$predictions[predictiontable$method=="Q_1_100_q95"][seq(1, sum(predictiontable$method=="Q_1_100_q95"), 2)+1]),
      y=rep(predictiontable$predictions[predictiontable$method==model],2),
      id=rep(paste(predictiontable$sim_id[predictiontable$method==model],"_",predictiontable$start[predictiontable$method==model],"_",predictiontable$method[predictiontable$method==model]),2),
      frequence=c(predictiontable$freq_benef[predictiontable$method=="Q_1_100_q95"][seq(1, sum(predictiontable$method=="Q_1_100_q95"), 2)],predictiontable$freq_benef[predictiontable$method=="Q_1_100_q95"][seq(1, sum(predictiontable$method=="Q_1_100_q95"), 2)+1]),
      model=model)
    newtab<-rbind(newtab,mod_table)
    }
  png(file.path(eval_dir,paste("corr_gnn_Q95_",project,".png",sep="")),width = 800, height = 600)
  print(ggplot(newtab,aes(x=x,y=y,group=id))+
    geom_point(aes(color=frequence,shape=model))+
    geom_line(aes(linetype=model))+
    labs(x="Q95",y="genomatnn")+
    theme(panel.grid = element_line(color = "#b5dfed",
                                    linewidth = 0.5,
                                    linetype = 2))+
    theme(text=element_text(size=20))+
    scale_color_gradientn(colours =perfectpalette))
  dev.off()
}

#main curves ####

if(simul_type=="AI"){
  if(length(unique(predictiontable$classifier))==2  & length(unique(predictiontable$predictions))>1){
    if(length(method_list)<=5){
      #creates roc
      png(file.path(eval_dir,paste("roc_",project,".png",sep="")),width = 600, height = 600)
      print(ggplot(graphtable[graphtable$graph=="roc",],aes(x=x,y=y,group=method,color=method,shape=method))+
        geom_line()+
        geom_point(data = subset(graphtable[graphtable$graph=="roc",], ideal), fill="black")+
        geom_point(data = subset(graphtable[graphtable$graph=="roc",], focus), fill="white")+
        scale_shape_manual(values = 21:25 )+
        xlim(0,1)+
        ylim(0,1)+
        geom_abline(aes(intercept=0,slope=1), linetype="dotted")+
        labs(x="Taux de faux positifs",y="Taux de vrais positifs",title="")+
        scale_colour_manual(values=cbPalette)+
        theme(legend.position="none")+
        theme(panel.grid = element_line(color = "#b5dfed",
                                        linewidth = 0.5,
                                        linetype = 2))+
        theme(text=element_text(size=20))+
        theme(plot.title = element_text(hjust = 0, vjust = -4.5)))
      dev.off()

      png(file.path(eval_dir,paste("roc_leg_",project,".png",sep="")),width = 900, height = 600)
      print(ggplot(graphtable[graphtable$graph=="roc",],aes(x=x,y=y,group=method,color=method,shape=method))+
        geom_line()+
        geom_point(data = subset(graphtable[graphtable$graph=="roc",], ideal), fill="black")+
        geom_point(data = subset(graphtable[graphtable$graph=="roc",], focus), fill="white")+
        scale_shape_manual(values = 21:25 )+
        xlim(0,1)+
        ylim(0,1)+
        geom_abline(aes(intercept=0,slope=1), linetype="dotted")+
        labs(x="Taux de faux positifs",y="Taux de vrais positifs",title="")+
        scale_colour_manual(values=cbPalette)+
        theme(panel.grid = element_line(color = "#b5dfed",
                                        linewidth = 0.5,
                                        linetype = 2))+
        theme(text=element_text(size=20))+
        theme(plot.title = element_text(hjust = 0, vjust = -4.5)))
      dev.off()

      #creates precision/recall curve
      png(file.path(eval_dir,paste("precrec_",project,".png",sep="")),width = 600, height = 600)
      print(ggplot(graphtable[graphtable$graph=="precrec",],aes(x=x,y=y,group=method,color=method,shape=method))+
        geom_path()+
        geom_point(data = subset(graphtable[graphtable$graph=="precrec",], ideal), fill="black")+
        geom_point(data = subset(graphtable[graphtable$graph=="precrec",], focus), fill="white")+
        scale_shape_manual(values = 21:25 )+
        xlim(0,1)+
        labs(x="Taux de vrais positifs",y="Précision",title="")+
        scale_colour_manual(values=cbPalette)+
        theme(legend.position="none")+
        theme(panel.grid = element_line(color = "#b5dfed",
                                        linewidth = 0.5,
                                        linetype = 2))+
        theme(text=element_text(size=20))+
        theme(plot.title = element_text(hjust = 0, vjust = -4.5)))
      dev.off()

      png(file.path(eval_dir,paste("precrec_leg_",project,".png",sep="")),width = 900, height = 600)
      print(ggplot(graphtable[graphtable$graph=="precrec",],aes(x=x,y=y,group=method,color=method,shape=method))+
        geom_path()+
        geom_point(data = subset(graphtable[graphtable$graph=="precrec",], ideal), fill="black")+
        geom_point(data = subset(graphtable[graphtable$graph=="precrec",], focus), fill="white")+
        scale_shape_manual(values = 21:25 )+
        xlim(0,1)+
        labs(x="Taux de vrais positifs",y="Précision",title="")+
        scale_colour_manual(values=cbPalette)+
        theme(panel.grid = element_line(color = "#b5dfed",
                                        linewidth = 0.5,
                                        linetype = 2))+
        theme(text=element_text(size=20))+
        theme(plot.title = element_text(hjust = 0, vjust = -4.5)))
      dev.off()

      #creates mccf1
      png(file.path(eval_dir,paste("mccf1_",project,".png",sep="")),width = 600, height = 600)
      print(ggplot(graphtable[graphtable$graph=="mccf1",],aes(x=x,y=y,alpha=0.7,group=method,color=method,shape=method))+
        geom_path()+
        geom_point(data = subset(graphtable[graphtable$graph=="mccf1",], ideal), fill="black")+
        geom_point(data = subset(graphtable[graphtable$graph=="mccf1",], focus), fill="white")+
        scale_shape_manual(values = 21:25 )+
        xlim(0,1)+
        ylim(0,1)+
        geom_hline(aes(yintercept=0.5), linetype="dotted")+
        labs(x="F1",y="MCC normalisé",title="")+
        scale_colour_manual(values=cbPalette)+
        theme(legend.position="none")+ 
        theme(panel.grid = element_line(color = "#b5dfed",
                                        linewidth = 0.5,
                                        linetype = 2))+
        theme(text=element_text(size=20))+
        theme(plot.title = element_text(hjust = 0, vjust = -4.5)))
      dev.off()

      png(file.path(eval_dir,paste("mccf1_leg_",project,".png",sep="")),width = 900, height = 600)
      print(ggplot(graphtable[graphtable$graph=="mccf1",],aes(x=x,y=y,group=method,color=method,shape=method,alpha=0.7))+
        geom_path()+
        geom_point(data = subset(graphtable[graphtable$graph=="mccf1",], ideal), fill="black")+
        geom_point(data = subset(graphtable[graphtable$graph=="mccf1",], focus), fill="white")+
        scale_shape_manual(values = 21:25 )+
        geom_hline(aes(yintercept=0.5), linetype="dotted")+
        xlim(0,1)+
        ylim(0,1)+
        labs(x="F1",y="MCC normalisé",title="")+
        scale_colour_manual(values=cbPalette)+
        theme(panel.grid = element_line(color = "#b5dfed",
                                        linewidth = 0.5,
                                        linetype = 2))+
        theme(text=element_text(size=20))+
        theme(plot.title = element_text(hjust = 0, vjust = -4.5)))
      dev.off()


      #creates zoomed mccf1
      png(file.path(eval_dir,paste("zoomed_mccf1_",project,".png",sep="")),width = 600, height = 600)
      print(ggplot(graphtable[graphtable$graph=="mccf1",],aes(x=x,y=y,group=method,color=method,shape=method))+
        geom_path()+
        geom_point(data = subset(graphtable[graphtable$graph=="mccf1",], ideal), fill="black")+
        geom_point(data = subset(graphtable[graphtable$graph=="mccf1",], focus), fill="white")+
        scale_shape_manual(values = 21:25 )+
        geom_hline(aes(yintercept=0.5), linetype="dotted")+
        labs(x="F1",y="MCC normalisé",title="")+
        scale_colour_manual(values=cbPalette)+
        theme(legend.position="none")+ 
        theme(panel.grid = element_line(color = "#b5dfed",
                                        linewidth = 0.5,
                                        linetype = 2))+
        theme(text=element_text(size=20))+
        theme(plot.title = element_text(hjust = 0, vjust = -4.5)))
      dev.off()

      png(file.path(eval_dir,paste("zoomed_mccf1_leg_",project,".png",sep="")),width = 900, height = 600)
      print(ggplot(graphtable[graphtable$graph=="mccf1",],aes(x=x,y=y,group=method,color=method,shape=method))+
        geom_path()+
        geom_point(data = subset(graphtable[graphtable$graph=="mccf1",], ideal), fill="black")+
        geom_point(data = subset(graphtable[graphtable$graph=="mccf1",], focus), fill="white")+
        scale_shape_manual(values = 21:25 )+
        geom_hline(aes(yintercept=0.5), linetype="dotted")+
        labs(x="F1",y="MCC normalisé",title="")+
        scale_colour_manual(values=cbPalette)+
        theme(panel.grid = element_line(color = "#b5dfed",
                                        linewidth = 0.5,
                                        linetype = 2))+
        theme(text=element_text(size=20))+
        theme(plot.title = element_text(hjust = 0, vjust = -4.5)))
      dev.off()
    }else{
      #creates roc
      png(file.path(eval_dir,paste("roc_",project,".png",sep="")),width = 600, height = 600)
      print(ggplot(graphtable[graphtable$graph=="roc",],aes(x=x,y=y,group=method,shape=method))+
        geom_line(aes(color=method))+
        geom_point(data = subset(graphtable[graphtable$graph=="roc",], ideal), color="black")+
        geom_point(data = subset(graphtable[graphtable$graph=="roc",], focus), color="white")+
        xlim(0,1)+
        ylim(0,1)+
        geom_abline(aes(intercept=0,slope=1), linetype="dotted")+
        labs(x="Taux de faux positifs",y="Taux de vrais positifs",title="")+
        scale_colour_manual(values=cbPalette)+
        theme(legend.position="none")+
        theme(panel.grid = element_line(color = "#b5dfed",
                                        linewidth = 0.5,
                                        linetype = 2))+
        theme(text=element_text(size=20))+
        theme(plot.title = element_text(hjust = 0, vjust = -4.5)))
      dev.off()

      png(file.path(eval_dir,paste("roc_leg_",project,".png",sep="")),width = 900, height = 600)
      print(ggplot(graphtable[graphtable$graph=="roc",],aes(x=x,y=y,group=method,shape=method))+
        geom_line(aes(color=method))+
        geom_point(data = subset(graphtable[graphtable$graph=="roc",], ideal), color="black")+
        geom_point(data = subset(graphtable[graphtable$graph=="roc",], focus), color="white")+
        xlim(0,1)+
        ylim(0,1)+
        geom_abline(aes(intercept=0,slope=1), linetype="dotted")+
        labs(x="Taux de faux positifs",y="Taux de vrais positifs",title="")+
        scale_colour_manual(values=cbPalette)+
        theme(panel.grid = element_line(color = "#b5dfed",
                                        linewidth = 0.5,
                                        linetype = 2))+
        theme(text=element_text(size=20))+
        theme(plot.title = element_text(hjust = 0, vjust = -4.5)))
      dev.off()

      #creates precision/recall curve
      png(file.path(eval_dir,paste("precrec_",project,".png",sep="")),width = 600, height = 600)
      print(ggplot(graphtable[graphtable$graph=="precrec",],aes(x=x,y=y,group=method,shape=method))+
        geom_path(aes(color=method))+
        geom_point(data = subset(graphtable[graphtable$graph=="precrec",], ideal), color="black")+
        geom_point(data = subset(graphtable[graphtable$graph=="precrec",], focus), color="white")+
        xlim(0,1)+
        labs(x="Taux de vrais positifs",y="Précision",title="")+
        scale_colour_manual(values=cbPalette)+
        theme(legend.position="none")+
        theme(panel.grid = element_line(color = "#b5dfed",
                                        linewidth = 0.5,
                                        linetype = 2))+
        theme(text=element_text(size=20))+
        theme(plot.title = element_text(hjust = 0, vjust = -4.5)))
      dev.off()

      png(file.path(eval_dir,paste("precrec_leg_",project,".png",sep="")),width = 900, height = 600)
      print(ggplot(graphtable[graphtable$graph=="precrec",],aes(x=x,y=y,group=method,shape=method))+
        geom_path(aes(color=method))+
        geom_point(data = subset(graphtable[graphtable$graph=="precrec",], ideal), color="black")+
        geom_point(data = subset(graphtable[graphtable$graph=="precrec",], focus), color="white")+
        xlim(0,1)+
        labs(x="Taux de vrais positifs",y="Précision",title="")+
        scale_colour_manual(values=cbPalette)+
        theme(panel.grid = element_line(color = "#b5dfed",
                                        linewidth = 0.5,
                                        linetype = 2))+
        theme(text=element_text(size=20))+
        theme(plot.title = element_text(hjust = 0, vjust = -4.5)))
      dev.off()

      #creates mccf1
      png(file.path(eval_dir,paste("mccf1_",project,".png",sep="")),width = 600, height = 600)
      print(ggplot(graphtable[graphtable$graph=="mccf1",],aes(x=x,y=y,group=method,shape=method,alpha=0.7))+
        geom_path(aes(color=method))+
        geom_point(data = subset(graphtable[graphtable$graph=="mccf1",], ideal), color="black")+
        geom_point(data = subset(graphtable[graphtable$graph=="mccf1",], focus), color="white")+
        xlim(0,1)+
        ylim(0,1)+
        geom_hline(aes(yintercept=0.5), linetype="dotted")+
        labs(x="F1",y="MCC normalisé",title="")+
        scale_colour_manual(values=cbPalette)+
        theme(legend.position="none")+ 
        theme(panel.grid = element_line(color = "#b5dfed",
                                        linewidth = 0.5,
                                        linetype = 2))+
        theme(text=element_text(size=20))+
        theme(plot.title = element_text(hjust = 0, vjust = -4.5)))
      dev.off()

      png(file.path(eval_dir,paste("mccf1_leg_",project,".png",sep="")),width = 900, height = 600)
      print(ggplot(graphtable[graphtable$graph=="mccf1",],aes(x=x,y=y,group=method,shape=method,alpha=0.7))+
        geom_path(aes(color=method))+
        geom_point(data = subset(graphtable[graphtable$graph=="mccf1",], ideal), color="black")+
        geom_point(data = subset(graphtable[graphtable$graph=="mccf1",], focus), color="white")+
        geom_hline(aes(yintercept=0.5), linetype="dotted")+
        xlim(0,1)+
        ylim(0,1)+
        labs(x="F1",y="MCC normalisé",title="")+
        scale_colour_manual(values=cbPalette)+
        theme(panel.grid = element_line(color = "#b5dfed",
                                        linewidth = 0.5,
                                        linetype = 2))+
        theme(text=element_text(size=20))+
        theme(plot.title = element_text(hjust = 0, vjust = -4.5)))
      dev.off()


      #creates zoomed mccf1
      png(file.path(eval_dir,paste("zoomed_mccf1_",project,".png",sep="")),width = 600, height = 600)
      print(ggplot(graphtable[graphtable$graph=="mccf1",],aes(x=x,y=y,group=method,shape=method))+
        geom_path(aes(color=method))+
        geom_point(data = subset(graphtable[graphtable$graph=="mccf1",], ideal), color="black")+
        geom_point(data = subset(graphtable[graphtable$graph=="mccf1",], focus), color="white")+
        geom_hline(aes(yintercept=0.5), linetype="dotted")+
        labs(x="F1",y="MCC normalisé",title="")+
        scale_colour_manual(values=cbPalette)+
        theme(legend.position="none")+ 
        theme(panel.grid = element_line(color = "#b5dfed",
                                        linewidth = 0.5,
                                        linetype = 2))+
        theme(text=element_text(size=20))+
        theme(plot.title = element_text(hjust = 0, vjust = -4.5)))
      dev.off()

      png(file.path(eval_dir,paste("zoomed_mccf1_leg_",project,".png",sep="")),width = 900, height = 600)
      print(ggplot(graphtable[graphtable$graph=="mccf1",],aes(x=x,y=y,group=method,shape=method))+
        geom_path(aes(color=method))+
        geom_point(data = subset(graphtable[graphtable$graph=="mccf1",], ideal), color="black")+
        geom_point(data = subset(graphtable[graphtable$graph=="mccf1",], focus), color="white")+
        geom_hline(aes(yintercept=0.5), linetype="dotted")+
        labs(x="F1",y="MCC normalisé",title="")+
        scale_colour_manual(values=cbPalette)+
        theme(panel.grid = element_line(color = "#b5dfed",
                                        linewidth = 0.5,
                                        linetype = 2))+
        theme(text=element_text(size=20))+
        theme(plot.title = element_text(hjust = 0, vjust = -4.5)))
      dev.off()
    }
  }else{#makes a bunch of empty plots if there is AI everywhere/no AI when there should be
    #roc
    png(file.path(eval_dir,paste("roc_",project,".png",sep="")),width = 600, height = 600)
    plot(1)
    dev.off()
    png(file.path(eval_dir,paste("roc_leg_",project,".png",sep="")),width = 900, height = 600)
    plot(1)
    dev.off()
    #precision/recall curve
    png(file.path(eval_dir,paste("precrec_",project,".png",sep="")),width = 600, height = 600)
    plot(1)
    dev.off()
    png(file.path(eval_dir,paste("precrec_leg_",project,".png",sep="")),width = 900, height = 600)
    plot(1)
    dev.off()
    #mccf1
    png(file.path(eval_dir,paste("mccf1_",project,".png",sep="")),width = 600, height = 600)
    plot(1)
    dev.off()
    png(file.path(eval_dir,paste("mccf1_leg_",project,".png",sep="")),width = 900, height = 600)
    plot(1)
    dev.off()
    #zoomed mccf1
    png(file.path(eval_dir,paste("zoomed_mccf1_",project,".png",sep="")),width = 600, height = 600)
    plot(1)
    dev.off()
    png(file.path(eval_dir,paste("zoomed_mccf1_leg_",project,".png",sep="")),width = 900, height = 600)
    plot(1)
    dev.off()
  }
}

for(method in method_list){
  apriori_threshold<-summarytable$apriori_threshold[summarytable$method==method]
  best_threshold<-summarytable$best_threshold[summarytable$method==method]

  #stat distribution ####
  png(file.path(eval_dir,paste("dist_",method,"_",project,".png",sep="")),width = 600, height = 600)
  print(ggplot(predictiontable[predictiontable$method==method,],aes(x=predictions,group=as.factor(classifier),alpha=0.5,fill=as.factor(classifier),color=as.factor(classifier)))+
    geom_density()+
    geom_vline(aes(xintercept=as.numeric(best_threshold)),linetype="dashed")+
    geom_vline(aes(xintercept=as.numeric(apriori_threshold)))+
    xlim(0,1)+
    labs(x="critère de classification",y="densité")+
    theme(panel.grid = element_line(color = "#b5dfed",
                                    linewidth = 0.5,
                                    linetype = 2))+
    theme(text=element_text(size=20))+
    theme(legend.position="none"))
  dev.off()

  #explicitly threshold dependant plots ####
  etdp<-data.frame()
  for(cutoff in sort(unique(predictiontable$predictions[predictiontable$method==method]))){
    newcut<-data.frame(
      cutoffs=cutoff,
      fpr=sum(predictiontable$predictions[predictiontable$method==method]>=cutoff & !predictiontable$classifier[predictiontable$method==method])/sum(!predictiontable$classifier[predictiontable$method==method]),
      fdr=sum(predictiontable$predictions[predictiontable$method==method]>=cutoff & !predictiontable$classifier[predictiontable$method==method])/sum(predictiontable$predictions[predictiontable$method==method]>=cutoff))
    etdp<-rbind(etdp,newcut)
  }

  #fpr/threshold plot
  png(file.path(eval_dir,paste("fpr_",method,"_",project,".png",sep="")),width = 600, height = 600)
  print(ggplot(data=etdp,aes(x=cutoffs,y=fpr))+
    geom_area(aes(alpha=0.7,fill="red",color="red"))+
    geom_vline(aes(xintercept=as.numeric(best_threshold)),linetype="dashed")+
    geom_vline(aes(xintercept=as.numeric(apriori_threshold)))+
    xlim(0,1)+
    ylim(0,1)+
    labs(x="Seuil",y="Taux de faux positifs")+
    theme(panel.grid = element_line(color = "#b5dfed",
                                    linewidth = 0.5,
                                    linetype = 2))+
    theme(text=element_text(size=20))+
    theme(legend.position="none"))
  dev.off()

  #fdr/threshold plot
  png(file.path(eval_dir,paste("fdr_",method,"_",project,".png",sep="")),width = 600, height = 600)
  print(ggplot(data=etdp,aes(x=cutoffs,y=fdr))+
    geom_area(aes(alpha=0.7,fill="red",color="red"))+
    geom_vline(aes(xintercept=as.numeric(best_threshold)),linetype="dashed")+
    geom_vline(aes(xintercept=as.numeric(apriori_threshold)))+
    xlim(0,1)+
    ylim(0,1)+
    labs(x="Seuil",y="Taux de fausse découverte")+
    theme(panel.grid = element_line(color = "#b5dfed",
                                    linewidth = 0.5,
                                    linetype = 2))+
    theme(text=element_text(size=20))+
    theme(legend.position="none"))
  dev.off()
}

#latent variables plot ####
if(simul_type=="AI"){
  if(length(unique(predictiontable$classifier))==2  & length(unique(predictiontable$predictions))>1){
    latentgraph<-data.frame()
    for(sim in simul_vector){
      max_freq<-max(predictiontable$freq_benef[predictiontable$sim_id==sim])
      fixation_time<-unique(predictiontable$fixation_g_r_back[predictiontable$sim_id==sim])
      prop_intro<-unique(predictiontable$mean_p1[predictiontable$sim_id==sim])
      for(method in method_list){
        newpoints<-data.frame(method=method,sim=sim,max_freq=max_freq,fixation_time=fixation_time,prop_intro=prop_intro)

        fpr0.05tp<-sum(predictiontable$fpr0.05_bin_pred[predictiontable$sim_id==sim & predictiontable$method==method] & predictiontable$classifier[predictiontable$sim_id==sim & predictiontable$method==method])
        fpr0.05tn<-sum(!predictiontable$fpr0.05_bin_pred[predictiontable$sim_id==sim & predictiontable$method==method] & !predictiontable$classifier[predictiontable$sim_id==sim & predictiontable$method==method])
        fpr0.05fp<-sum(predictiontable$fpr0.05_bin_pred[predictiontable$sim_id==sim & predictiontable$method==method] & !predictiontable$classifier[predictiontable$sim_id==sim & predictiontable$method==method])
        fpr0.05fn<-sum(!predictiontable$fpr0.05_bin_pred[predictiontable$sim_id==sim & predictiontable$method==method] & predictiontable$classifier[predictiontable$sim_id==sim & predictiontable$method==method])

        newpoints$fpr0.05_sensitivity<-(fpr0.05tp)/(fpr0.05tp+fpr0.05fn)
        newpoints$fpr0.05_specificity<-(fpr0.05tn)/(fpr0.05tn+fpr0.05fp)
        
        latentgraph<-rbind(latentgraph,newpoints)
      }
    }

    if(length(unique(latentgraph$max_freq[latentgraph$max_freq<1 & latentgraph$max_freq>0]))>=3){
      for(method in method_list){
        newfit<-fitme( fpr0.05_sensitivity ~ 1+Matern(1|max_freq), data=latentgraph[latentgraph$max_freq<1 & latentgraph$max_freq>0 & latentgraph$method==method,], family=binomial())
        latentgraph$fitted_pow_freq[latentgraph$max_freq<1 & latentgraph$max_freq>0 & latentgraph$method==method]<-predict(newfit)
      }
    }else{
      latentgraph$fitted_pow_freq[latentgraph$max_freq<1 & latentgraph$max_freq>0]<-as.numeric(NA)
    }

    if(length(unique(latentgraph$fixation_time[latentgraph$max_freq==1]))>=3){
      for(method in method_list){
        newfit<-fitme( fpr0.05_sensitivity ~ 1+Matern(1|fixation_time), data=latentgraph[latentgraph$max_freq==1 & latentgraph$method==method,], family=binomial())
        latentgraph$fitted_pow_time[latentgraph$max_freq==1 & latentgraph$method==method]<-predict(newfit)
      }
    }else{
      latentgraph$fitted_pow_time[latentgraph$max_freq==1]<-as.numeric(NA)
    }

    if(length(unique(latentgraph$prop_intro))>=3){
      for(method in method_list){
        newfit<-fitme( fpr0.05_sensitivity ~ 1+Matern(1|prop_intro), data=latentgraph[latentgraph$method==method,], family=binomial())
        latentgraph$fitted_pow_prop[latentgraph$method==method]<-predict(newfit)
      }
    }else{
      latentgraph$fitted_pow_prop<-as.numeric(NA)
    }

    png(file.path(eval_dir,paste("freq_",project,".png",sep="")),width = 900, height = 600)
    if(sum(latentgraph$max_freq!=1)!=0){
      print(ggplot(latentgraph[latentgraph$max_freq!=1,],aes(x=max_freq,group=method,color=method,shape=method))+
        xlim(0,1)+
        ylim(0,1)+
        geom_point(aes(y=fpr0.05_sensitivity))+
        geom_line(aes(y=fitted_pow_freq))+
        labs(x="Fréquence de l'allèle bénéfique",y="Sensibilité",title="")+
        scale_colour_manual(values=cbPalette)+
        theme(panel.grid = element_line(color = "#b5dfed",
                                        linewidth = 0.5,
                                        linetype = 2))+
        theme(text=element_text(size=20))+
        theme(plot.title = element_text(hjust = 0, vjust = -4.5)))
    }else{
      print(ggplot())
    }
    dev.off()
    
    png(file.path(eval_dir,paste("time_",project,".png",sep="")),width = 900, height = 600)
    if(sum(latentgraph$max_freq==1)!=0){
      print(ggplot(latentgraph[latentgraph$max_freq==1,],aes(x=fixation_time,group=method,color=method,shape=method))+
        ylim(0,1)+
        geom_point(aes(y=fpr0.05_sensitivity))+
        geom_line(aes(y=fitted_pow_time))+
        labs(x="Temps de fixation",y="Sensibilité",title="")+
        scale_colour_manual(values=cbPalette)+
        theme(panel.grid = element_line(color = "#b5dfed",
                                        linewidth = 0.5,
                                        linetype = 2))+
        theme(text=element_text(size=20))+
        theme(plot.title = element_text(hjust = 0, vjust = -4.5)))
    }else{
      print(ggplot())
    }
    dev.off()

    png(file.path(eval_dir,paste("introprop_",project,".png",sep="")),width = 900, height = 600)
    print(ggplot(latentgraph,aes(x=prop_intro,group=method,color=method,shape=method))+
      ylim(0,1)+
      geom_point(aes(y=fpr0.05_sensitivity))+
      geom_line(aes(y=fitted_pow_prop))+
      labs(x="Proportion d'introgression",y="Sensibilité",title="")+
      scale_colour_manual(values=cbPalette)+
      theme(panel.grid = element_line(color = "#b5dfed",
                                      linewidth = 0.5,
                                      linetype = 2))+
      theme(text=element_text(size=20))+
      theme(plot.title = element_text(hjust = 0, vjust = -4.5)))
    dev.off()
  }else{
    png(file.path(eval_dir,paste("freq_",project,".png",sep="")),width = 600, height = 600)
    plot(1)
    dev.off()

    png(file.path(eval_dir,paste("time_",project,".png",sep="")),width = 600, height = 600)
    plot(1)
    dev.off()

    png(file.path(eval_dir,paste("introprop_",project,".png",sep="")),width = 600, height = 600)
    plot(1)
    dev.off()
  }
}

#tables ####
write.table(summarytable, file =paste(eval_dir, paste("performance_metrics_",project,".csv",sep=""),sep="/") , col.names = TRUE, row.names = FALSE)
write.table(predictiontable, file =paste(eval_dir, paste("predictions_",project,".csv",sep=""),sep="/"), col.names = TRUE, row.names = FALSE)
write.table(graphtable, file =paste(eval_dir, paste("justin_case_",project,".csv",sep=""),sep="/"), col.names = TRUE, row.names = FALSE)


compare_end <- Sys.time()
duration_compare <- as.numeric(compare_end - compare_start)
print(paste0("(getparams) Duration sims file creation : ", duration_compare, " s "))