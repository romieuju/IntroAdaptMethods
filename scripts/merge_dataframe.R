#    IntroAdapt: joint inference of introgression and selection
#    Copyright (C) 2024 Jules Romieu, ISEM/CBGP, CNRS/INRAE. 
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


##Script used to merge tables containing parameter values, latent variables and summary statistics and to obtain tables containing information of interest per window and for each simulated genomic data, named reftable_analysis_project and reftable_wind_analysis_project respectively. 
#Warning: the term reftable is not appropriate as these tables are not necessarily going to be used to train a method, but these names come from a script whose original purpose was to create reftables from simulated data and these file names have not been changed.  

# R package used for this script
library(extraDistr, quietly=TRUE)
library(ini, quietly = TRUE)
source("scripts/introadapt.R")

# read command line arguments
args = commandArgs(trailingOnly = TRUE)
if (length(args) == 0){
  stop("One positional argument necessary in the command (options file)")
  quit(save="no")
} else options_file = args[1]

create_reftable_start <- Sys.time()


# read options file
opts             <- get_project_options(options_file)
Settings         <- opts$Settings
Inference        <- opts$Inference
project          <- Settings$project
analysis         <- Settings$analysis
results_dir      <- Settings$results_dir
project_dir_bis  <- Settings$project_dir
wind_samples     <- Inference$wind_reftable_type
simu             <- Settings$num_of_sims

#Future tables merging parameter values, latent variables and summary statistics by genomic data or by window:
df_reftable      <- data.frame()
df_reftable_wind <- data.frame()

#For each genomic data :
for (sim in seq_len(simu)){
    #Folder
    simu_dir             <- paste(project,"sim" ,sim, sep="_")
    simu_dir             <- paste(results_dir, analysis, project, paste(analysis,simu_dir,sep="_"), sep="/")
    project_dir          <- paste(results_dir, analysis, project, sep="/")
   
    ##Futur "reftable" (dataframe with parameters, latent var and summary statistics values)
    #For all genomic data 
    reftable_file        <- paste(project_dir, paste("reftable_",analysis,"_",project,".csv",sep=""),sep="/")
    #For all genomic data and all window
    reftable_file_wind   <- paste(project_dir, paste("reftable_wind_",analysis,"_",project,".csv",sep=""),sep="/")
    #For all window of the interest genomic data (sim)
    reftable_by_sim_file <- paste(project_dir, paste("reftable_",analysis,"_",project,"_sim_",sim,".csv",sep=""),sep="/")
    
    # Read genomic data ini file :
    simu_file            <- paste(simu_dir, paste(analysis,"_",project,"_","sim_", sim ,".ini",sep=""), sep="/")
    opts_bis             <- get_sim_options(simu_file)
    Genome               <- opts_bis$Genome

    ##Files :
    #Parameters files :
    parameters_file    <- paste(simu_dir, paste("data_parameters_",analysis,"_",project,"_","sim_", sim,".csv",sep=""), sep="/")
    #Latent variable files:
    latente_var_1_file <- paste(simu_dir,paste("forwsim_latent_variable_",analysis,"_",project,"_","sim_", sim ,".txt",sep=""), sep="/")
    adv_mut_file       <- paste(simu_dir, paste("forwsim_freq_mut_don_in_rec_",analysis,"_",project,"_","sim_", sim ,".txt",sep=""), sep="/")
    #Summary stat file :
    statistics_file    <- paste(simu_dir, paste("sum_stat_",analysis,"_",project,"_","sim_", sim ,".csv",sep=""),sep="/")
    statistics_file_2  <- paste(simu_dir, paste("Sum_Stat_Mut_TreeSeq_",analysis,"_",project,"_","sim_", sim ,".csv",sep=""),sep="/")

    #read file : 
    parameters    <- read.csv(parameters_file)
    latente_var_1 <- read.csv(latente_var_1_file)
    adv_mut       <- read.csv(adv_mut_file)
    statistics    <- read.csv(statistics_file)
    statistics_2  <- read.csv(statistics_file_2)

    #create reftable dataframe and save the parameters :
    ref_table <- parameters

    #Get reftable column name :
    col_names <- colnames(ref_table)

    #remove subnumeraire column : 
    ref_table <- subset(ref_table, select = col_names[-1])
    #Merge refttable and ia_in_simu 
    ref_table <- cbind(ref_table, latente_var_1)
    ref_table <- cbind(ref_table, statistics)

    #Write reftable file for each genomic data :
    if (Inference$reftable_type == "estimation_by_sim" | Inference$reftable_type == "all"){
        write.table(ref_table, file = reftable_by_sim_file, col.names = TRUE, row.names = FALSE)
    }

    #Dataframe for all genomic data
    df_reftable <- rbind(df_reftable, ref_table)
    
    #Start creation of dataframe with all window summary stat value :
    #add on the .ini wind_samples (nbr of sample neutral window) and a parameters for creat this kind of reftable
    #dataframe with summary stat by window
    stat       <- statistics_2
    list_r     <- Genome$msprime_r_map_rates
    list_r_num <- as.numeric(strsplit(list_r," ")[[1]])
    keep_r     <- round(list_r_num,2)!=0.69
    list_r_num <- list_r_num[keep_r]

    list_r_pos <- Genome$msprime_r_map_positions
    list_int   <- as.integer(strsplit(list_r_pos," ")[[1]])
    list_int_1 <- list_int[-1]
    list_int_1 <- list_int_1[keep_r]
    list_int_0 <- as.integer(append(0,list_int_1))
    list_int_0 <- list_int_0[-c(length(list_int_0))]

    stat_wind <- stat[,c("start","end")]
    for (interv_recomb in 1:length(list_int_1)){
        stat_wind[which(stat_wind[,"start"]>=list_int_0[interv_recomb] & stat_wind[,"end"]<list_int_1[interv_recomb]),"r"]<-list_r_num[interv_recomb]
    }
    stat["r"]         <- stat_wind["r"]
    #data frame with advantageous mut introgressed from donor to recipient pop with N*s>1.0
    mut               <- adv_mut
    #advantageous mut position + sort
    vect_mut          <- mut$position
    #check if ia mutation exist :
    if (is.nan(mut$position[1])==FALSE){
        vect_mut_selec_coef  <- mut$selection_coef
        #vector for save futur window with IA index
        #vect_index_row_to_del <- c()
        vect_index_row_to_del <- list()
        #for each position
        for (adv_mut in vect_mut){
            window_of_interest    <- stat[stat$start<=adv_mut & stat$end>=adv_mut,]#keep columns where is the adv mut
            index                 <- rownames(window_of_interest) #save its index
            vect_index_row_to_del <- list(vect_index_row_to_del,as.integer(index)) #save the index in the index vector
        }
        #remove duplicate index
        vect_index_row_to_del      <- unlist(vect_index_row_to_del)
        vect_uniq_index_row_to_del <- unique(vect_index_row_to_del)
        #create a data frame with uniquely the neutral window
        params_stat_reftable           <- 0
        #vector with unique window 
        vect_index_agg                 <- as.vector(aggregate(vect_mut_selec_coef~vect_index_row_to_del,FUN=mean)["vect_index_row_to_del"][[1]])
        #selection coef mean by window
        vect_mean_s                     <- as.vector(aggregate(vect_mut_selec_coef~vect_index_row_to_del,FUN=mean)["vect_mut_selec_coef"][[1]])
        #selection coef max by window
        vect_max_s                      <- as.vector(aggregate(vect_mut_selec_coef~vect_index_row_to_del,FUN=max)["vect_mut_selec_coef"][[1]])
        #number of IA mut by window
        vect_count_mut                  <- as.vector(aggregate(vect_mut_selec_coef~vect_index_row_to_del,FUN=length)["vect_mut_selec_coef"][[1]])
        #Add this columns to reftable
        reftable_wind_adv               <- stat[vect_index_agg,]
        reftable_wind_adv["IA"]         <- 1
        reftable_wind_adv["mean_s"]     <- vect_mean_s
        reftable_wind_adv["max_s"]      <- vect_max_s
        reftable_wind_adv["ia_mut_nbr"] <- vect_count_mut
        #window without ia mut parameters 
        stat_without_mut               <- stat[-vect_uniq_index_row_to_del,]
        #Check if neutral window exist :
        if (nrow(stat_without_mut)!=0){
            stat_without_mut["IA"]         <- 0
            stat_without_mut["mean_s"]     <- 0.0
            stat_without_mut["max_s"]      <- 0.0
            stat_without_mut["ia_mut_nbr"] <- 0
            reftable_wind_neut             <- 0
            if (wind_samples == "model_train"){#if the user want a ratio 1:2 (introgressed window/window) 
                nbr_wind_adv <- length(reftable_wind_adv$start)
                nbr_wind_neu <- length(stat_without_mut$start)
                if (nbr_wind_adv*2<=nbr_wind_neu){
                    reftable_wind_neut <- sample_n(stat_without_mut,nbr_wind_adv*2)
                }else if (nbr_wind_adv<=nbr_wind_neu){
                    reftable_wind_neut <- sample_n(stat_without_mut,nbr_wind_adv)
                }else if (nbr_wind_adv>nbr_wind_neu){
                    reftable_wind_neut <- sample_n(stat_without_mut,nbr_wind_neu)
                }
            }else if (wind_samples == "classifier"){#if the user want a finit number of neutral window
                reftable_wind_neut <- stat_without_mut
            }
        } else {#if no neutral window in the simulation the reftable is only made of introgressed window 
            print("REFTABLE SCRIPT WARNING : Simulation without neutral window, only made of introgressed windows")
            reftable_wind_neut <- stat_without_mut #void neutral table
        }
        #merge neutral and advantageous dataframe
        reftable_all         <- rbind(reftable_wind_adv,reftable_wind_neut)
        #Sort reftable by start column
        reftable_all_order   <- reftable_all[order(reftable_all$start),]
        #merge reftable with parameters table
        latente_var_1_bis    <-latente_var_1[,c("mean_p1","prop_mut_do_in_rec")]
        parameters_bis       <- merge(parameters, latente_var_1_bis)
        params_stat_reftable <- merge(parameters_bis, reftable_all_order)
        if (is.na(latente_var_1["fixation_d_g"])){
            params_stat_reftable["fixation_g_d_back"] <- "NA"
        }else{
            params_stat_reftable["fixation_g_d_back"] <- (parameters["g_forward"]-latente_var_1["fixation_d_g"])*parameters["scaling_factor"]
        }   
        if (is.na(latente_var_1["fixation_r_g"])){
            params_stat_reftable["fixation_g_r_back"] <- "NA"
        }else{
            params_stat_reftable["fixation_g_r_back"] <- (parameters["g_forward"]-latente_var_1["fixation_r_g"])*parameters["scaling_factor"]
        }
    }else{#for simulation without IA :
        if (wind_samples == "model_train"){
            colname_ref                  <- c(colnames(parameters),"mean_p1", "prop_mut_do_in_rec", colnames(stat),"IA", "mean_s","max_s", "ia_mut_nbr")
            reftable_all_order           <- data.frame(matrix(nrow = 0, ncol = length(colname_ref)))
            colnames(reftable_all_order) <- colname_ref
            params_stat_reftable         <- reftable_all_order
            params_stat_reftable["fixation_g_d_back"] <- "NA"
            params_stat_reftable["fixation_g_r_back"] <- "NA"
        }else if (wind_samples == "classifier"){
            stat["IA"]                                <- 0
            stat["mean_s"]                            <- 0.0
            stat["max_s"]                             <- 0.0
            stat["ia_mut_nbr"]                        < -0
            latente_var_1_bis                         <- latente_var_1[,c("mean_p1","prop_mut_do_in_rec")]
            parameters_bis                            <- merge(parameters, latente_var_1_bis)
            params_stat_reftable                      <- merge(parameters_bis, stat)
            params_stat_reftable["fixation_g_d_back"] <- "NA"
            params_stat_reftable["fixation_g_r_back"] <- "NA"
        }
    }
    df_reftable_wind <- rbind(df_reftable_wind, params_stat_reftable)
}
if (Inference$reftable_type == "classification" | Inference$reftable_type == "all"){
    write.table(df_reftable_wind, file = reftable_file_wind, col.names = TRUE, row.names = FALSE)
}
if (Inference$reftable_type == "estimation" | Inference$reftable_type == "all"){
    write.table(df_reftable, file = reftable_file, col.names = TRUE, row.names = FALSE)
}



create_reftable_end <- Sys.time()
duration_create_reftable <- as.numeric(create_reftable_end - create_reftable_start)
print(paste0("(reftable) Duration reftable creation : ", duration_create_reftable, " s "))