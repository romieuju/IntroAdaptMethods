#####################Violin plot##############################################################################################################
#script that produces graphs showing the distribution of scores by method as a function of the distance from the window to the AI window 
#####################GOOD FOR 2L###############################################################################################################
library(ggplot2)
library(dplyr)
#path of comparison folder (folder created after end_evaluation rule in snakemake file)
eval_dir            <- "/home/romieuju/results/classification/human_2chro_s_inter/comparison"
setwd(eval_dir)
#name of the project
project <-"human_2chro_s_inter"
#name of the prediction file with all score value by windows by methods after end_evaluation rule execution
pred_data <- read.csv(file = "predictions_human_2chro_s_inter.csv", sep =" ")
#path to neutral (NI) project prediction file with all score value by windows by methods after end_evaluation rule execution
pred_data_IN_true <- read.csv(file = "/home/romieuju/results/human_2chro_neutral/comparison/predictions_human_2chro_neutral.csv", sep =" ")
#####to calculate volcano max########
max_vol <- max(pred_data[pred_data$method=="vol_rel",]$predictions_no_norm, pred_data_IN_true[pred_data_IN_true$method=="vol_rel",]$predictions_no_norm)
####auroc by non-Ai window#######
#auroc by adjacent windows class (from article_roc_curve_by_adjacent_window_exemple.R)
auroc_adj   <- read.csv(file = "auroc_AI_VS_ADJ_by_window_human_2chro_s_inter.csv", sep =" ")
names(auroc_adj)<-c("method", "50000",  "100000",  "150000", "200000",  "250000", "300000",  "350000", "400000",  "450000", "500000")
#auroc for second chromosome windows (from article_roc_curve_by_non_ai_window_type_exemple.R)
auroc_chro2 <- read.csv(file = "auroc_AI_VS_CHRO2_all_human_2chro_s_inter.csv", sep =" ")
names(auroc_chro2)[2] <- "auroc_chro2"
#auroc for independante neutral introgression windows (from article_roc_curve_by_non_ai_window_type_exemple.R)
auro_ni     <- read.csv(file = "auroc_AI_VS_ADJ_all_human_2chro_s_inter.csv", sep =" ")
names(auro_ni)[2] <- "auroc_ni"
auro<-merge(auroc_adj, auroc_chro2, by="method")
auro<-merge(auro, auro_ni, by="method")
pred_dat_IA <-pred_data[pred_data$start<1000000,] 
pred_data_IN <-pred_data[pred_data$start>=1000000,]
pred_data   <-pred_dat_IA
pred_data_IN$dist_to_mut <-100000000
pred_data_IN_true$dist_to_mut <-200000000
pred_data$dist_to_mut<-abs(pred_data$start-500000)
pred_data[pred_data$method=="gnn_1",]$start = pred_data[pred_data$method=="gnn_1",]$start-1
pred_data[pred_data$method=="gnn_2",]$start = pred_data[pred_data$method=="gnn_2",]$start-1
v_dist_name<-c("0","50","100","150","200","250","300","350","400","450","500","","IN")
colorp <- c("#0072B2","#00945C","#00945C","#00945C","#00945C","#00945C","#00945C","#00945C","#00945C","#00945C","#00945C","#00A300","#00CC00","#00CC00")
colorp_g <- c("#0072B2","#00945C","#00945C","#00945C","#00945C","#00945C","#00A300","#00CC00","#00CC00")
v_limit <- c(0,50000,100000,150000,200000, 250000, 300000, 350000, 400000, 450000, 500000, 550000, 100000000, 200000000)
v_dist_name_g<-c("0kb","100","200","300","400","500","","IN")
v_limit_g <- c(0,100000,200000, 300000, 400000,500000,600000, 100000000, 200000000)
r<-1.2e-8
v_dist_cM_g<-c()
for (dist in v_limit_g[-c(length(v_limit_g)-2, length(v_limit_g)-1,length(v_limit_g))]){
  cM<-as.character(round(50*log(1/(1-2*(dist*r))),2))
  v_dist_cM_g<-c(v_dist_cM_g, cM)
}
v_dist_cM<-c()
for (dist in v_limit[-c(length(v_limit)-2, length(v_limit)-1,length(v_limit))]){
  cM<-as.character(round(50*log(1/(1-2*(dist*r))),2))
  v_dist_cM<-c(v_dist_cM, cM)
}
v_dist_cM_g <- c(v_dist_cM_g, "", "chro 2", "NI")
v_dist_cM   <- c(v_dist_cM, "", "chro 2", "NI")
pred_data[pred_data$method=="gnn_1",]$dist_to_mut = abs(pred_data[pred_data$method=="gnn_1",]$start-500000)
pred_data[pred_data$method=="gnn_2",]$dist_to_mut = abs(pred_data[pred_data$method=="gnn_2",]$start-500000)
pred_data <- rbind(pred_data, pred_data_IN, pred_data_IN_true)
auroc_data<-data.frame()
pred_data_bis<-data.frame()
for (met in unique(pred_data$method)){
  if (met=="gnn_1" | met=="gnn_2"){
    x_names<-v_dist_cM_g
    limit_value <- v_limit_g
    color_p<-colorp_g
  }else{
    x_names<-v_dist_cM
    limit_value <- v_limit
    color_p<-colorp
  }
  if(met!="vol_rel"){
    pred_data_bis<-pred_data[, c("method", "dist_to_mut", "predictions_no_norm")]
    if(met =="gnn_1" | met=="gnn_2"){
      vec_auroc     <- as.numeric(t(auro[auro$method==met,][-1]))
      vec_auroc_adj <-vec_auroc[1:(length(vec_auroc)-2)] 
      vec_auroc_chro2_in <-vec_auroc[(length(vec_auroc)-1):length(vec_auroc)]
      n <- c(1:length(vec_auroc_adj))
      n_odd<-n%%2 != 0
      vec_auroc_adj_odd <- vec_auroc_adj[n_odd]
      vec_auroc_final<-c(vec_auroc_adj_odd,vec_auroc_chro2_in)
      auroc_data<-data.frame("dist_to_mut"= limit_value[-c(1,length(limit_value)-2)],"auroc"=vec_auroc_final)
    }else{
      auroc_data<-data.frame("dist_to_mut"= limit_value[-c(1,length(limit_value)-2)],"auroc"=as.numeric(t(auro[auro$method==met,][-1])))
    }
    print(met)
    print(limit_value[-c(1,length(limit_value)-2)])
    print(auroc_data)
    test_1 <- (ggplot(pred_data_bis[pred_data_bis$method==met,], aes(factor(dist_to_mut), as.numeric(predictions_no_norm)))+ 
                 geom_violin(alpha=0.3, na.rm = FALSE, scale = "width") +
                 geom_boxplot(width=0.05) +
                 scale_colour_manual(values =color_p, labels=x_names)+
                 geom_point(data = auroc_data, aes(x = factor(dist_to_mut), y = auroc), shape = 24 ,size = 5, colour = "black", fill="gray") +
                 ylim(-0.01,1.01)+
                 scale_x_discrete(limits=factor(limit_value), labels=x_names)+
                 labs(x="Distance from the window under AI", y="Score value")+
                 scale_fill_brewer(palette="Dark2") + theme_minimal() +theme(text=element_text(size=27, color="black"), legend.position = "none" ))
    print(test_1)
    ggsave(paste("Score_value_vs_distance_to_IA_wind_",met,"_",project,"_violin.pdf",sep=""), test_1, width = 11.69, height = 8.27, units = "in", dpi = 300, device = "pdf")
  }else{
    pred_data_bis<-pred_data[, c("method", "dist_to_mut", "predictions_no_norm")]
    auroc_data<-data.frame("dist_to_mut"= limit_value[-c(1,length(limit_value)-2)],"auroc"=as.numeric(t(auro[auro$method==met,][-1])))
    test_1 <- (ggplot(pred_data_bis[pred_data_bis$method==met,], aes(factor(dist_to_mut), as.numeric(predictions_no_norm/max_vol)))+ 
                 geom_violin(alpha=0.3, na.rm = FALSE, scale = "width") +
                 geom_boxplot(width=0.05) +
                 geom_point(data = auroc_data, aes(x = factor(dist_to_mut), y = auroc), shape = 24 ,size = 5, colour = "black", fill="gray") +
                 scale_colour_manual(values =color_p, labels=x_names)+
                 ylim(-0.01,1.01)+
                 scale_x_discrete(limits=factor(limit_value), labels=x_names)+
                 labs(x="Distance from the window under AI", y="Score value")+
                 scale_fill_brewer(palette="Dark2") + theme_minimal() +theme(text=element_text(size=27, color="black"), legend.position = "none" ))
    print(test_1)
    ggsave(paste("Score_value_vs_distance_to_IA_wind_",met,"_",project,"_violin.pdf",sep=""), test_1, width = 11.69, height = 8.27, units = "in", dpi = 300, device = "pdf")
  }
}
