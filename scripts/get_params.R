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

# This script is used to draws the evolutionary parameters from the prior, create simulation files for backward (msprime) and forward simulation (SLiM) and create parameters data frame. 

# R package used for this script
library(extraDistr, quietly=TRUE)
library(ini, quietly = TRUE)
source("scripts/introadapt.R")

# Variable with start time of getparams script
getparams_start <- Sys.time()

# read command line arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0){
  stop("Two positional arguments necessary in the command (RDS options file and batch number)")
  quit(save="no")
} else {
  options_file = args[1]
}

#Save the ini information on variables :
opts       = get_project_options(options_file)
Settings   = opts$Settings
Model      = opts$Model
Statistics = opts$Statistics
Sample     = opts$Sample
Genome     = opts$Genome

# set seed for random number generator and get seeds for msprime and SLiM
seed_coal = as.integer(round(runif(Settings$num_of_sims, 1, Settings$seed)))
seed_forw = as.integer(round(runif(Settings$num_of_sims, 1, Settings$seed)))
seed_mut  = as.integer(round(runif(Settings$num_of_sims, 1, Settings$seed)))
if (Settings$verbose >= 10) write(paste("Seed (msprime-coalescence)", seq_len(Settings$num_of_sims), ":", seed_coal), stdout())
if (Settings$verbose >= 10) write(paste("Seed (SLiM)", seq_len(Settings$num_of_sims), ":", seed_forw), stdout())
if (Settings$verbose >= 10) write(paste("Seed (msprime-mutation)", seq_len(Settings$num_of_sims), ":", seed_mut), stdout())



#For each genomic data to simulated do : 
for (sim in seq_len(Settings$num_of_sims)){
  # create simulation directory
  project     <-Settings$project
  project_dir <-Settings$project_dir
  analysis    <-Settings$analysis
  job_dir     <-paste(analysis,project,"sim", sim,sep="_")
  job_dir_path<-paste(project_dir, job_dir,sep="/")
  dir.create(job_dir_path, showWarnings = FALSE)
  if (Settings$verbose >= 10) write(paste("Job folder:", job_dir), stdout())

  #Define recombination map for msprime and SliM if chromosome number > 1
  #msprime recombination map:
  msprime_r_map_rates = paste(Genome$msprime_r_map_rates, collapse = ",")
  msprime_r_map_rates_num = as.numeric(unlist(strsplit(msprime_r_map_rates," ")))
  msprime_r_map_positions = paste(Genome$msprime_r_map_positions, collapse = ",")
  
  #SLiM recombination map:
  slim_r_map_rates = paste(Genome$slim_r_map_rates, collapse = ",")
  slim_r_map_rates_num = as.numeric(unlist(strsplit(slim_r_map_rates," ")))
  slim_r_map_positions    = paste(strsplit(Genome$slim_r_map_positions, split = " ")[[1]], collapse=",")

  #Random pick of demographic and evolutionnary parameters 
  if (Settings$verbose >= 10) write(paste("Simulation", sim), stdout())
  popAncO               = as.integer(eval(parse(text = Model$popAncO)))
  if ((is.vector(popAncO) == TRUE) && (length(popAncO) == Settings$num_of_sims)){
    popAncO             = popAncO[sim]
  }
  popO                  = as.integer(eval(parse(text = Model$popO)))
  if ((is.vector(popO) == TRUE) && (length(popO) == Settings$num_of_sims)){
    popO             = popO[sim]
  }
  popAncAB              = as.integer(eval(parse(text = Model$popAncAB)))
  if ((is.vector(popAncAB) == TRUE) && (length(popAncAB) == Settings$num_of_sims)){
    popAncAB             = popAncAB[sim]
  }
  popAncA               = as.integer(eval(parse(text = Model$popAncA)))
  if ((is.vector(popAncA) == TRUE) && (length(popAncA) == Settings$num_of_sims)){
    popAncA             = popAncA[sim]
  }
  popAncB               = as.integer(eval(parse(text = Model$popAncB)))
  if ((is.vector(popAncB) == TRUE) && (length(popAncB) == Settings$num_of_sims)){
    popAncB             = popAncB[sim]
  }
  popA1                 = as.integer(eval(parse(text = Model$popA1)))
  if ((is.vector(popA1) == TRUE) && (length(popA1) == Settings$num_of_sims)){
    popA1             = popA1[sim]
  }
  popB1                 = as.integer(eval(parse(text = Model$popB1)))
  if ((is.vector(popB1) == TRUE) && (length(popB1) == Settings$num_of_sims)){
    popB1             = popB1[sim]
  }
  popA2                 = as.integer(eval(parse(text = Model$popA2)))
  if ((is.vector(popA2) == TRUE) && (length(popA2) == Settings$num_of_sims)){
    popA2             = popA2[sim]
  }
  popB2                 = as.integer(eval(parse(text = Model$popB2)))
  if ((is.vector(popA2) == TRUE) && (length(popB2) == Settings$num_of_sims)){
    popB2             = popB2[sim]
  }
  migration_rate_r      = as.numeric(eval(parse(text = Model$migration_rate_r)))
  if ((is.vector(migration_rate_r) == TRUE) && (length(migration_rate_r) == Settings$num_of_sims)){
    migration_rate_r             = migration_rate_r[sim]
  }
  coef_selec_before_mig = as.numeric(eval(parse(text = Model$coef_selec_before_mig)))
  if ((is.vector(coef_selec_before_mig) == TRUE) && (length(coef_selec_before_mig) == Settings$num_of_sims)){
    coef_selec_before_mig             = coef_selec_before_mig[sim]
  }
  mu_advantageous       = as.numeric(eval(parse(text = Model$mu_advantageous)))
  if ((is.vector(mu_advantageous) == TRUE) && (length(mu_advantageous) == Settings$num_of_sims)){
    mu_advantageous             = mu_advantageous[sim]
  }
  mu_total              = as.numeric(eval(parse(text = Model$mu_total)))
  if ((is.vector(mu_total) == TRUE) && (length(mu_total) == Settings$num_of_sims)){
    mu_total             = mu_total[sim]
  }
  mu_neutral            = as.numeric(mu_total-mu_advantageous)
  if ((is.vector(mu_neutral) == TRUE) && (length(mu_neutral) == Settings$num_of_sims)){
    mu_neutral             = mu_neutral[sim]
  }
  #Define populations split (divergence time) :
  #first, draw randomly from the prior and the distribution the divergence time 
  generation_split_OA         = as.numeric(eval(parse(text=Model$generation_split_OA)))
  generation_split_AB         = as.numeric(eval(parse(text=Model$generation_split_AB)))
  generation_split_A          = as.numeric(eval(parse(text=Model$generation_split_A)))
  generations_migration_start = as.numeric(eval(parse(text=Model$generations_migration_start)))
    repeat{
        if(generation_split_AB>generation_split_A & generation_split_A>generations_migration_start){
        break
        }
        generation_split_AB         <- as.numeric(eval(parse(text=Model$generation_split_AB)))
        generation_split_A          <- as.numeric(eval(parse(text=Model$generation_split_A)))
        generations_migration_start <- as.numeric(eval(parse(text=Model$generations_migration_start)))
    }
    # Recipient divergence time have to be egal or inf to the donor divergence time
    generation_split_B = as.numeric(eval(parse(text=Model$generation_split_B))) 
    repeat{
        if(generation_split_A>=generation_split_B){
        break
        }
        generation_split_B = as.numeric(eval(parse(text=Model$generation_split_B)))
    }
    # the end of the Migration time have to be inf to the start of the migration
    if (as.numeric(eval(parse(text=Model$generations_migration_end))) == generations_migration_start) {
        generations_migration_end = generations_migration_start-1
    }else{
        generations_migration_end = as.numeric(eval(parse(text=Model$generations_migration_end)))
    }
    generations_forward = as.numeric(eval(parse(text=Model$generations_forward)))
    if (as.numeric(eval(parse(text=Model$sample_time_d))) == generations_forward-1) {
        sample_time_d = generations_forward-1
    }else{
        sample_time_d = as.numeric(eval(parse(text=Model$sample_time_d)))
    }
    if (as.numeric(eval(parse(text=Model$sample_time_d_s))) == generations_forward-1) {
        sample_time_d_s = generations_forward-1
    }else{
        sample_time_d_s = as.numeric(eval(parse(text=Model$sample_time_d_s)))
    }



  #Apply the scaling factor to the parameter values:
  scaling_factor                = Model$scaling_factor
  if (scaling_factor > 1){
    popAncO                     = round(popAncO/scaling_factor, 0)
    popO                        = round(popO/scaling_factor, 0)
    popAncAB                    = round(popAncAB/scaling_factor, 0)
    popAncA                     = round(popAncA/scaling_factor, 0)
    popAncB                     = round(popAncB/scaling_factor, 0)
    popA1                       = round(popA1/scaling_factor, 0)
    popB1                       = round(popB1/scaling_factor, 0)
    popA2                       = round(popA2/scaling_factor, 0)
    popB2                       = round(popB2/scaling_factor, 0)
    if (Model$add_adv_mutation_style == "manually"){
      coef_selec_before_mig       = coef_selec_before_mig*scaling_factor
      migration_rate_r            = migration_rate_r
    }else{
      coef_selec_before_mig       = coef_selec_before_mig*scaling_factor
      migration_rate_r            = migration_rate_r
    }
    mu_advantageous             = mu_advantageous*scaling_factor
    mu_total                    = mu_total*scaling_factor
    mu_neutral                  = as.numeric(mu_total-mu_advantageous)
    s                           = which(slim_r_map_rates_num!=0.5)
    slim_r_map_rates_num[s]     = (1/2)*(1-(1-2*slim_r_map_rates_num[s])**scaling_factor)
    m                           = which(round(msprime_r_map_rates_num,2)!=0.69)
    msprime_r_map_rates_num[m]  = (1/2)*(1-(1-2*msprime_r_map_rates_num[m])**scaling_factor)
    generation_split_OA         = round(generation_split_OA/scaling_factor, 0)
    generation_split_AB         = round(generation_split_AB/scaling_factor, 0)
    generation_split_A          = round(generation_split_A/scaling_factor, 0)
    generation_split_B          = round(generation_split_B/scaling_factor, 0)
    if (generations_migration_end == generations_migration_start-1){
      generations_migration_end   = round(generations_migration_start/scaling_factor, 0)-1
    }else{
      generations_migration_end   = round(generations_migration_end/scaling_factor, 0)
    }
    generations_migration_start = round(generations_migration_start/scaling_factor, 0)
    generations_forward         = round(generations_forward/scaling_factor, 0)
    sample_time_d = round(sample_time_d/scaling_factor, 0)
    if (sample_time_d == 0){
      sample_time_d = 1
    }
    sample_time_d_s = round(sample_time_d_s/scaling_factor, 0)
    if (sample_time_d_s == 0){
      sample_time_d_s = 1
    }
  }

    #Define forward generation time for forward simulation:
    generation_split_OA_back<-generation_split_OA
    generation_split_AB_back<-generation_split_AB
    generation_split_A_back<-generation_split_A
    generation_split_B_back<-generation_split_B
    generations_migration_start_back<-generations_migration_start
    generations_migration_end_back<-generations_migration_end
    sample_time_d_back<-sample_time_d
    sample_time_d_s_back<-sample_time_d_s
    generations_forward_back<-generations_forward

    temps_backward <- c(generation_split_OA, generation_split_AB, generation_split_A, generation_split_B, generations_migration_start, generations_migration_end, sample_time_d_s, sample_time_d,0)
    temps_backward_bis<-temps_backward[temps_backward<=generations_forward]
    temps_backward_keep<-temps_backward[temps_backward>generations_forward]
    temps_forward <- rev(cumsum(rev(diff(temps_backward_bis)))) + temps_backward_bis[1]
    temps_forward[1]<-generations_forward[1]
    temps_forward<-c(temps_backward_keep, temps_forward)

    generation_split_OA<-temps_forward[1]
    generation_split_AB<-temps_forward[2]    
    generation_split_A<-temps_forward[3]
    generation_split_B<-temps_forward[4]     
    generations_migration_start<-temps_forward[5]
    generations_migration_end<-temps_forward[6]
    sample_time_d_s<-temps_forward[7]
    sample_time_d<-temps_forward[8]

    generation_split_AB_back <-(generation_split_AB_back - generation_split_AB_back)+1
    generation_split_OA_back <-(generation_split_OA_back - generation_split_AB)
    
    

  #Define region under selection in donor pop in SliM :
  if (Model$add_adv_mutation_style == "mutation_rate"){
    L                          = as.integer(Genome$L)
    proportion_under_selection = Model$proportion_under_selection
    if (proportion_under_selection !=1.0){
      genomic_region_size        = Model$genomic_region_size
      max_region_proportion      = Model$max_region_proportion
      Nbr_Inter                  = as.integer(L/genomic_region_size)+1
      V_Inter                    = seq(from=genomic_region_size-1, to=L, by = genomic_region_size)
      V_binom                    = rbinom(Nbr_Inter, 1, as.numeric(proportion_under_selection[1]))==1
      v_selected_region          = as.integer(V_Inter[V_binom])
      v_selected_region          = na.omit(v_selected_region)
      selected_region_number     = length(v_selected_region)
      repeat{
        if (((selected_region_number*genomic_region_size)/L<max_region_proportion) & (selected_region_number>0)){
          break
        }
        V_binom                    = rbinom(Nbr_Inter, 1, as.numeric(proportion_under_selection[1]))==1
        v_selected_region          = as.integer(V_Inter[V_binom])
        selected_region_number     = length(v_selected_region)
      }
      v_selected_region          = paste(v_selected_region,collapse=',')
      vect_mut_pos               = paste(0,collapse=',')
      generation_mut             = paste(0,collapse=',')
    }else{
      genomic_region_size        = L
      proportion_under_selection = 1.0
      selected_region_number     = 1
      v_selected_region          = paste(L,collapse=',')
    }
    mix_vect_pos=0
  }else if (Model$add_adv_mutation_style == "manually"){
    if (Model$proportion_under_selection<=1){
      genomic_region_size        = as.integer(Genome$L)+1
      proportion_under_selection = Model$proportion_under_selection
      proportion_under_selection = as.numeric(proportion_under_selection[1])
      max_region_proportion      = as.numeric(Model$max_region_proportion)
      window_size                = Statistics$window_size
      selected_region_number     = 1
      v_selected_region          = paste(genomic_region_size,collapse=',')
      ##calcul the windows number with the window and genome length
      nbr_win<-genomic_region_size/window_size
      ##calcul the number of windows with futur advantageous mutation 
      nbr_mut<-as.numeric(nbr_win)*max_region_proportion
      ##calcul the position of the first mutation in the first window
      pos_start<-window_size*proportion_under_selection
      #dist_between_mut<-genomic_region_size/nbr_mut
      dist_between_mut<-window_size
      ##Create a vector with all the possible advantageous mutation positions 
      vect_mut_pos<-seq(pos_start,genomic_region_size,dist_between_mut)
    
      #sample the mutation position
      mix_vect_pos<-sample(vect_mut_pos,nbr_mut)
      mix_vect_pos<-paste(mix_vect_pos,collapse=',')
      mu_advantageous <- 0.0
      mu_neutral <- mu_total
      generation_mut<-as.numeric(eval(parse(text=Model$generations_forward)))-as.numeric(eval(parse(text=Model$generation_mutation)))
      generation_mut<-paste(generation_mut,collapse=',')
    }else{
      proportion_under_selection <- 0.0
      genomic_region_size <- as.integer(Genome$L)+1
      selected_region_number <- 1
      v_selected_region <- paste(Genome$L,collapse=',')
      mix_vect_pos <- paste(strsplit(Model$proportion_under_selection, split = " ")[[1]], collapse=",")
      window_size <- Statistics$window_size
      mu_advantageous <- 0.0
      mu_neutral <- mu_total
      generation_mut<-as.numeric(eval(parse(text=Model$generations_forward)))-as.numeric(eval(parse(text=Model$generation_mutation)))
      generation_mut<-paste(generation_mut,collapse=',')
    }
  }

  # write sim_n.ini file (backward simulation ini file)
  sim_ini <- list()
  sim_ini[["Simulation"]] = list(sim                          = sim) 
  sim_ini[["Demography"]] = list(popAncO                      = popAncO,
                                 popO                         = popO,
                                 popAncAB                     = popAncAB,
                                 popAncA                      = popAncA,
                                 popAncB                      = popAncB,
                                 popA1                        = popA1,
                                 popB1                        = popB1,
                                 popA2                        = popA2,
                                 popB2                        = popB2,
                                 generation_split_OA          = as.integer(generation_split_OA_back),
                                 generation_split_AB          = as.integer(generation_split_AB_back),
                                 generation_split_A           = as.integer(generation_split_A),
                                 generation_split_B           = as.integer(generation_split_B),
                                 generations_migration_start  = as.integer(generations_migration_start),
                                 generations_migration_end    = as.integer(generations_migration_end),
                                 sample_time_d                = as.integer(sample_time_d),
                                 sample_time_d_s              = as.integer(sample_time_d_s),
                                 generations_forward          = as.integer(generations_forward),
                                 coalescence_split            = Model$coalescence_split,
                                 slim_split                   = Model$slim_split)
  sim_ini[["Genome"]]     = list(mu_neutral                   = mu_neutral,
                                 mu_advantageous              = mu_advantageous,
                                 ttratio                      = 2.0,
                                 seq_error                    = 0.005,
                                 msprime_r_map_positions      = msprime_r_map_positions,
                                 msprime_r_map_rates          = paste(msprime_r_map_rates_num,collapse=" "))
  sim_ini[["Seeds"]]      = list(seed_coal                    = seed_coal[sim],
                                 seed_mut                     = seed_mut[sim])
  
  #Write simulation ini file:
  sim_ini_file <- paste0(job_dir_path, "/",analysis,"_",project,"_sim_", sim,".ini")
  write.ini(sim_ini, sim_ini_file)
  
  
  # write source file for SLiM (sim_n.eidos)
  source4slim <- paste0("setSeed(", seed_forw[sim], ");\n",
                        "defineConstant(\"popO\",", popO, ");\n",
                        "defineConstant(\"popAncAB\",",  popAncAB, ");\n",
                        "defineConstant(\"popAncA\",",  popAncA, ");\n",
                        "defineConstant(\"popAncB\",",  popAncB, ");\n",
                        "defineConstant(\"popA1\",",  popA1, ");\n",
                        "defineConstant(\"popB1\",",  popB1, ");\n",
                        "defineConstant(\"popA2\",",  popA2, ");\n",
                        "defineConstant(\"popB2\",",  popB2, ");\n",
                        "defineConstant(\"recipient\",",  Model$recipient, ");\n",
                        "defineConstant(\"donor\",",  Model$donor, ");\n",
                        "defineConstant(\"slim_split\",", Model$slim_split, ");\n",
                        "defineConstant(\"migration_rate_r\",",  migration_rate_r, ");\n",
                        "defineConstant(\"generation_split_AB\",", generation_split_AB, ");\n",
                        "defineConstant(\"generation_split_A\",", generation_split_A, ");\n",
                        "defineConstant(\"generation_split_B\",", generation_split_B, ");\n",
                        "defineConstant(\"generations_forward\",", generations_forward, ");\n",
                        "defineConstant(\"sample_time_d\",", sample_time_d, ");\n",
                        "defineConstant(\"sample_time_d_s\",", sample_time_d_s, ");\n",
                        "defineConstant(\"migration_start\",", generations_migration_start, ");\n",
                        "defineConstant(\"migration_end\",", generations_migration_end, ");\n",
                        "defineConstant(\"s_before_m\",", coef_selec_before_mig, ");\n",
                        "defineConstant(\"i_allele_state_d\",", Model$i_allele_state_d, ");\n",
                        "defineConstant(\"i\",", sim, ");\n",
                        "defineConstant(\"project\",\"", Settings$project, "\");\n",
                        "defineConstant(\"analysis\",\"", Settings$analysis, "\");\n",
                        "defineConstant(\"results_dir\",\"", Settings$results_dir, "\");\n",
                        "defineConstant(\"L\",", Genome$L, ");\n",
                        "defineConstant(\"genomic_region_size\",", genomic_region_size, ");\n",
                        "defineConstant(\"proportion_under_selection\",", proportion_under_selection, ");\n",
                        "defineConstant(\"v_selected_region\",c(",v_selected_region, "));\n",
                        "defineConstant(\"v_mut_position\",c(",mix_vect_pos, "));\n",
                        "defineConstant(\"v_mut_generation\",c(",generation_mut, "));\n",
                        "defineConstant(\"mu_neutral\",",  mu_neutral, ");\n",
                        "defineConstant(\"mu_advantageous\",",  mu_advantageous, ");\n",
                        "defineConstant(\"ends\",c(",slim_r_map_positions, "));\n",
                        "defineConstant(\"rates\",c(",paste(slim_r_map_rates_num,collapse=","), "));\n",
                        "defineConstant(\"n_samples\",c(",paste(strsplit(Sample$n_samples, split = " ")[[1]], collapse=","), "));\n")
  #Write simulation eidos file
  write(source4slim, file = paste0(job_dir_path, "/",analysis,"_",project,"_sim_", sim, ".eidos"))
  
  #write Dataframe with all the demographique and the evolutionary parameters values (one row by genetic data, warning: values with scaling factor)
  data_simu = data.frame(num_sim                    = sim,
                         add_mutation_type          = Model$add_adv_mutation_style,
                         scaling_factor             = Model$scaling_factor,
                         N_AncO                     = popAncO,
                         N_O                        = popO,
                         N_AncAB                    = popAncAB,
                         N_AncA                     = popAncA,
                         N_AncB                     = popAncB,
                         N_A1                       = popA1,
                         N_B1                       = popB1,
                         N_A2                       = popA2,
                         N_B2                       = popB2,
                         split_AncO                 = generation_split_OA_back,
                         split_AB                   = generation_split_AB_back,
                         split_A                    = generation_split_A_back,
                         split_B                    = generation_split_B_back,
                         migration_start            = generations_migration_start_back,
                         migration_end              = generations_migration_end_back,
                         g_sample_d                 = sample_time_d_back,
                         g_sample_d_s               = sample_time_d_s_back,
                         g_forward                  = generations_forward_back,
                         genome_l                   = as.integer(Genome$L)+1,
                         migration_rate             = migration_rate_r,
                         select_coef                = coef_selec_before_mig,
                         mu_neutral                 = mu_neutral,
                         mu_advantageous            = mu_advantageous,
                         mu_total                   = mu_total,
                         recomb_rate                = mean(slim_r_map_rates_num[slim_r_map_rates_num!=0.5]),
                         genomic_region_size        = genomic_region_size, 
                         proportion_under_selection = proportion_under_selection,
                         position_mut_av            = Model$proportion_under_selection,
                         selected_region_number     = selected_region_number
                         )
  #Write parameters values dataframe:
  csv_name = paste0(job_dir_path,"/data_parameters_",analysis,"_",project,"_sim_",sim, ".csv")
  write.csv(data_simu,file=csv_name, row.names = TRUE)

}

#Duration time of the script:
getparams_end <- Sys.time()
getparams_duration <- as.numeric(getparams_end - getparams_start)
print(paste0("(getparams) Duration sims file creation : ", getparams_duration, " s "))