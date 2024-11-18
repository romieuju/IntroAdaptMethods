#    IntroAdapt: joint inference of introgression and selection
#    Copyright (C) 2022 Jules Romieu, ISEM/CBGP, CNRS/INRAE. 
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

### This script is used to write the informations save in the genome file with information about recombination rate, genome lenght and chromosome number in the new ini file in the folder of the project.

library(ini, quietly = TRUE)
source("scripts/introadapt.R")

# Variable with start time of setproject
setproject_start <- Sys.time()

# read command line arguments
args = commandArgs(trailingOnly = TRUE)
if (length(args) == 0){
  stop("One positional argument necessary in the command (options file)")
  quit(save="no")
} else options_file = args[1]

# read options file
opts     = get_project_options(options_file)
Settings = opts$Settings
Model    = opts$Model
Stat     = opts$Statistics
Infer    = opts$Inference
Tools    =  opts$Tools

# print script info to screen
print_info("setproject.R", Settings$verbose, project = Settings$project)

# set seed for random number generator
set.seed(Settings$seed)

# Name of analysis and project directory (result_dir/analysis_dir/project_dir)
analysis_dir <- paste(Settings$results_dir, Settings$analysis, sep="/")
project_dir  <- paste(analysis_dir, Settings$project, sep="/")
if (Settings$verbose >= 10) write(paste("Project folder:", project_dir), stdout())

#Recovers the sample sizes for each population (n_sample) and the size of the total sample (size)
n_samples    <- Model$n_samples
n_sample_vec <- as.vector(as.numeric(unlist(strsplit(n_samples," "))))
size         <- sum(n_sample_vec )
Sample <-   list(size      = size,
                 n_samples = n_samples)

#Recovers genome information from genome file (recombination rate, genome lenght, chromosome number) :
#If user doesn't use a genome file and have juste one chromosome else r_rates = file
if (Model$r_map_rate!="file"){
  r_rates           <- as.numeric(eval(parse(text =Model$r_map_rate)))
  nchr              <- 1
  L                 <- Model$genome_l
  L_msprime         <- as.integer(as.integer(L)+1)
  chr_ends          <- L
  msprime_positions <- paste("0",L_msprime,sep=" ")
  slim_positions    <- as.integer(L)
  Genome <-   list(nchr        = nchr,
                 chr_ends      = chr_ends,
                 L             = L,
                 msprime_r_map = list(rates     = r_rates,
                                      positions = msprime_positions),
                 slim_r_map    = list(rates     = r_rates,
                                      positions = slim_positions))
}else{
  Genome <- read_genome_info(Settings$genome_file)
}
if (Settings$verbose >= 10) write(paste("Genome length:", Genome$L+1), stdout()) 


# Create project options file with information from genome file. 
project_ini <- list()
project_ini[["Settings"]] = list(seed           = Settings$seed,
                                 simul_type     = Settings$simul_type,
                                 project        = Settings$project,
                                 config_file    = Settings$config_file,
                                 analysis       = Settings$analysis,
                                 results_dir    = Settings$results_dir,
                                 verbose        = Settings$verbose,
                                 num_of_sims    = Settings$num_of_sims,
                                 project_dir    = project_dir)
project_ini[["Model"]]    = list(populations                 = Model$populations,
                                 popAncO                     = Model$popAncO,
                                 popO                        = Model$popO,
                                 popAncAB                    = Model$popAncAB,
                                 popA1                       = Model$popA1,
                                 popB1                       = Model$popB1,
                                 popAncA                     = Model$popAncA,
                                 popAncB                     = Model$popAncB,
                                 popA2                       = Model$popA2,
                                 popB2                       = Model$popB2,
                                 recipient                     = Model$recipient,
                                 donor                      = Model$donor,
                                 coalescence_split           = Model$coalescence_split,
                                 slim_split                  = Model$slim_split,
                                 slim_conditional            = Model$slim_conditional,
                                 add_adv_mutation_style      = Model$add_adv_mutation_style,
                                 migration_rate_r            = Model$migration_rate_r,
                                 generation_split_OA         = Model$generation_split_OA,
                                 generation_split_AB         = Model$generation_split_AB,
                                 generation_split_A          = Model$generation_split_A,
                                 generation_split_B          = Model$generation_split_B,
                                 generation_mutation         = Model$generation_mutation,
                                 generations_migration_start = Model$generations_migration_start,
                                 generations_migration_end   = Model$generations_migration_end,
                                 generations_forward         = Model$generations_forward,
                                 sample_time_d_s             = Model$sample_time_d_s,
                                 sample_time_d               = Model$sample_time_d,
                                 coef_selec_before_mig       = Model$coef_selec_before_mig,
                                 mu_advantageous             = Model$mu_advantageous,
                                 mu_total                    = Model$mu_total,
                                 i_allele_state_d            = Model$i_allele_state_d,
                                 scaling_factor              = Model$scaling_factor,
                                 genomic_region_size         = Model$genomic_region_size,
                                 proportion_under_selection  = paste(Model$proportion_under_selection, collapse=" "),
                                 max_region_proportion       = Model$max_region_proportion)
project_ini[["Statistics"]] = list(zns                = Stat$zns,
                                   donor_pop_stat    = Stat$donor_pop_stat,
                                   prop_i_by_r_interv = Stat$prop_i_by_r_interv,
                                   isphased           = Stat$isphased,
                                   window_size   = Stat$window_size,
                                   window_start  = Stat$window_start,
                                   window_end    = Stat$window_end,
                                   window_step   = Stat$window_step
                                   )
project_ini[["Sample"]]   = list(size         = Sample$size,
                                 n_samples = paste(Sample$n_samples, collapse=" "))
project_ini[["Genome"]]   = list(nchr                    = Genome$nchr,
                                 L                       = Genome$L,
                                 chr_ends                = paste(Genome$chr_ends, collapse=" "),
                                 msprime_r_map_positions = paste(Genome$msprime_r_map$positions, collapse=" "),
                                 msprime_r_map_rates     = paste(Genome$msprime_r_map$rates, collapse=" "),
                                 slim_r_map_positions    = paste(Genome$slim_r_map$positions, collapse=" "),
                                 slim_r_map_rates        = paste(Genome$slim_r_map$rates, collapse=" "))
project_ini[["Inference"]]= list(volcano_threshold       = Infer$volcano_threshold,
                                 genomatnn_threshold     = Infer$genomatnn_threshold,
                                 maladapt_threshold      = Infer$maladapt_threshold,
                                 q95_threshold           = Infer$q95_threshold,
                                 volcano_switch          = Infer$volcano_switch,
                                 maladapt_switch         = Infer$maladapt_switch,
                                 maladapt_trained_path   = Infer$maladapt_trained_path,
                                 maladapt_features       = Infer$maladapt_features,
                                 genomatnn_switch        = Infer$genomatnn_switch,
                                 genomatnn_trained_mod   = Infer$genomatnn_trained_mod,
                                 wind_reftable_type      = Infer$wind_reftable_type,
                                 reftable_type           = Infer$reftable_type)

#Save the new ini file with sample and genome information :
write.ini(project_ini, file = paste0(project_dir,"/project_options.ini"))

#Check that the number of chromosome is the same everywhere:
if(Stat$nb_chr !=Genome$nchr & Infer$volcano_switch=="On"){
  stop("Inconsistent chromosome number between the two .ini files")
  }

setproject_end <- Sys.time()
setproject_duration <- as.numeric(setproject_end - setproject_start)
print(paste0("(setproject) Duration project file creation : ", setproject_duration, " s "))