#    TimeAdapt: joint inference of demography and selection
#    Copyright (C) 2024  Jules Romieu, CBGP/UM, CNRS/INRAE
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
#

###Script containing the functions used in the R scripts in the pipeline

library(dplyr)

### PRINT INFO ######################################################################################
print_info = function(script_name, verbose, project = NA, batch = NA, sim = NA)
{
  if (verbose >= 1) write("#########################################",stdout())
  if (verbose >= 0)
  {
    info_message = paste0("TimeAdapt - ", script_name)
    if (!is.na(project)) info_message = paste0(info_message, " - project ", project)
    if (!is.na(batch)) info_message = paste0(info_message, " - batch ", batch)
    if (!is.na(sim)) info_message = paste0(info_message, " - simulation ", sim)
    write(info_message, stdout())
  }
  if (verbose >= 1)
  {
    write("by Romieu Jules", stdout())
    write("INRAE & Montpllier Universitet", stdout())
    write("jules.romieu@umontpllier.fr", stdout())
    write("#########################################", stdout())
  } 
}

### GET PROJECT OPTIONS ###############################################################################
get_project_options = function(options_file){
  options = read.ini(options_file)
  # SETTINGS
  options$Settings$project        = options$Settings$project
  options$Settings$config_file    = options$Settings$config_file
  options$Settings$analysis       = options$Settings$analysis
  options$Settings$results_dir    = options$Settings$results_dir
  options$Settings$genome_file    = options$Settings$genome_file
  options$Settings$verbose        = as.integer(options$Settings$verbose)
  options$Settings$num_of_sims    = as.integer(options$Settings$num_of_sims)
  options$Settings$seed           = as.integer(options$Settings$seed)
  options$Settings$simul_type     = options$Settings$simul_type

  # MODEL
  options$Model$populations                 = as.integer(options$Model$populations)
  options$Model$popAncO                     = options$Model$popAncO
  options$Model$popO                        = options$Model$popO
  options$Model$popAncAB                    = options$Model$popAncAB
  options$Model$popA1                       = options$Model$popA1
  options$Model$popB1                       = options$Model$popB1
  options$Model$popAncA                     = options$Model$popAncA
  options$Model$popAncB                     = options$Model$popAncB
  options$Model$popA2                       = options$Model$popA2
  options$Model$popB2                       = options$Model$popB2
  options$Model$recipient                   = options$Model$recipient
  options$Model$donor                       = options$Model$donor
  options$Model$coalescence_split           = options$Model$coalescence_split
  options$Model$slim_split                  = as.integer(options$Model$slim_split)
  options$Model$slim_conditional            = options$Model$slim_conditional
  options$Models$add_adv_mutation_style     = options$Models$add_adv_mutation_style
  options$Model$migration_rate_r            = options$Model$migration_rate_r
  options$Model$generation_split_OA         = options$Model$generation_split_OA
  options$Model$generation_split_AB         = options$Model$generation_split_AB 
  options$Model$generation_split_A          = options$Model$generation_split_A
  options$Model$generation_split_B          = options$Model$generation_split_B
  options$Model$generations_migration_start = options$Model$generations_migration_start
  options$Model$generations_migration_end   = options$Model$generations_migration_end
  options$Model$generations_forward         = options$Model$generations_forward
  options$Model$sample_time_d_s             = options$Model$sample_time_d_s
  options$Model$sample_time_d               = options$Model$sample_time_d
  options$Model$generation_mutation         = options$Model$generation_mutation
  options$Model$coef_selec_before_mig       = options$Model$coef_selec_before_mig
  options$Model$mu_advantageous             = options$Model$mu_advantageous
  options$Model$mu_total                    = options$Model$mu_total
  options$Model$i_allele_state_d            = options$Model$i_allele_state_d
  options$Model$genome_l                    = options$Model$genome_l
  options$Model$r_map_rate                  = options$Model$r_map_rate
  options$Model$n_samples                   = options$Model$n_samples
  options$Model$scaling_factor              = as.integer(options$Model$scaling_factor)
  options$Model$genomic_region_size         = as.integer(options$Model$genomic_region_size)
  options$Model$proportion_under_selection  = options$Model$proportion_under_selection
  options$Model$max_region_proportion       = as.numeric(options$Model$max_region_proportion)
  
  # INFERENCE METHODS

  #Volcanofinder
  options$Inference$volcano_switch         = options$Inference$volcano_switch
  options$Inference$volcano_threshold      = as.numeric(options$Inference$volcano_threshold)
  #MaLAdapt
  options$Inference$maladapt_switch        = options$Inference$maladapt_switch
  options$Inference$maladapt_trained_path  = options$Inference$maladapt_trained_path
  options$Inference$maladapt_features      = options$Inference$maladapt_features
  options$Inference$maladapt_threshold     = as.numeric(options$Inference$maladapt_threshold)
  #Genomatnn
  options$Inference$genomatnn_switch       = options$Inference$genomatnn_switch
  options$Inference$genomatnn_trained_mod  = options$Inference$genomatnn_trained_mod
  options$Inference$genomatnn_threshold    = as.numeric(options$Inference$genomatnn_threshold)
  #MaLAdapt table
  options$Inference$wind_reftable_type     = options$Inference$wind_reftable_type
  options$Inference$reftable_type          = options$Inference$reftable_type

  #Q95
  options$Inference$q95_threshold          = as.numeric(options$Inference$q95_threshold)

  
  # STATISTICS
  options$Statistics$zns                    = options$Statistics$zns
  options$Statistics$donor_pop_stat         = options$Statistics$donor_pop_stat
  options$Statistics$prop_i_by_r_interv     = options$Statistics$prop_i_by_r_interv
  options$Statistics$isphased               = options$Statistics$isphased
  options$Statistics$window_size            = as.integer(options$Statistics$window_size)
  options$Statistics$window_start           = as.integer(options$Statistics$window_start)
  options$Statistics$window_end             = as.integer(options$Statistics$window_end)
  options$Statistics$window_step            = as.integer(options$Statistics$window_step)

  return(options)
}

### GET SIMULATION OPTIONS ###############################################################################
get_sim_options = function(options_file){
  option = read.ini(options_file)
  # SIMULATION
  option$Simulation$sim               = as.integer(option$Simulation$sim)
  option$Simulation$generation        = as.integer(option$Simulation$generation)

  # DEMOGRAPHY
  option$Demography$popAncO                     = as.integer(option$Demography$popAncO)
  option$Demography$popO                        = as.integer(option$Demography$popO)
  option$Demography$popAncAB                    = as.integer(option$Demography$popAncAB)
  option$Demography$popA1                       = as.integer(option$Demography$popA1)
  option$Demography$popB1                       = as.integer(option$Demography$popB1)
  option$Demography$popAncA                     = as.integer(option$Demography$popAncA)
  option$Demography$popAncB                     = as.integer(option$Demography$popAncB)
  option$Demography$popA2                       = as.integer(option$Demography$popA2)
  option$Demography$popB2                       = as.integer(option$Demography$popB2)
  option$Demography$generation_split_OA         = as.integer(option$Demography$generation_split_OA)
  option$Demography$generation_split_AB         = as.integer(option$Demography$generation_split_AB)
  option$Demography$generation_split_A          = as.integer(option$Demography$generation_split_A)
  option$Demography$generation_split_B          = as.integer(option$Demography$generation_split_B)
  option$Demography$generations_forward         = as.integer(option$Demography$generations_forward)
  option$Demography$generations_migration_start = as.integer(option$Demography$generations_migration_start)
  option$Demography$generations_migration_end   = as.integer(option$Demography$generations_migration_end)
  option$Demography$sample_time_d               = as.integer(option$Model$sample_time_d)
  option$Demography$coalescence_split           = option$Demography$coalescence_split
  option$Demography$slim_split                  = as.integer(option$Demography$slim_split)

  # GENOME
  option$Genome$mu_neutral              = option$Genome$mu_neutral
  option$Genome$mu_advantageous         = option$Genome$mu_advantageous
  option$Genome$ttratio                 = option$Genome$ttratio
  option$Genome$seq_error               = option$Genome$seq_error
  option$Genome$msprime_r_map_positions = option$Genome$msprime_r_map_positions
  option$Genome$msprime_r_map_rates     = option$Genome$msprime_r_map_rates

  # SEEDS
  option$Seeds$seed_coal              = as.integer(option$Seeds$seed_coal)
  option$Seeds$seed_mut               = as.integer(option$Seeds$seed_mut)

  return(option)
}


### CHECK FILE HEADER ######################################################################################
check_file_header = function(expected_header, file_header){
  missing =! is.element(expected_header, file_header)
  if (any(missing)) return(list(ok = FALSE, missing = expected_header[missing]))
  else              return(list(ok = TRUE,  missing = NA)) 
}

### READ GENOME INFO ######################################################################################
read_genome_info = function(file){
  info = read.table(file,header=T)
  expected_header = c("Chromosome","Position","Recombination_rate")
  header_check = check_file_header(expected_header, file_header = colnames(info))
  if(header_check$ok){
    rates            = info$Recombination_rate
    positions        = info$Position
    nchr             = nlevels(as.factor(info$Chromosome))
    chr_ends_index   = c(which(diff(info$Chromosome)!=0), length(info$Chromosome))
    rescaling_values = c(0,cumsum(info$Position[chr_ends_index]))  
    chromo = 1
    for (i in seq_along(info$Position)){
      positions[i] = info$Position[i] + rescaling_values[chromo]
      if (any(i == chr_ends_index)) chromo = chromo + 1
    }
    chr_ends = as.integer(positions[chr_ends_index])

    slim_positions =  numeric()
    slim_rates =  numeric()
    msprime_positions =  0
    msprime_rates =  numeric()
    for (chromo in seq_len(nchr)){
      slim_positions    = c(slim_positions,    positions[which(info$Chromosome==chromo)]-1)
      slim_rates        = c(slim_rates,        rates[which(info$Chromosome==chromo)])
      msprime_positions = c(msprime_positions, positions[which(info$Chromosome==chromo)])
      msprime_rates     = c(msprime_rates,     rates[which(info$Chromosome==chromo)])
      if (chromo!=max(info$Chromosome)){
        slim_positions    = c(slim_positions,    positions[chr_ends_index[chromo]])
        slim_rates        = c(slim_rates,        0.5)
        msprime_positions = c(msprime_positions, positions[chr_ends_index[chromo]]+1)
        msprime_rates     = c(msprime_rates,     log(2))
      }
    } 
    slim_positions = as.integer(slim_positions)
    msprime_positions = as.integer(msprime_positions)
    L = as.integer(slim_positions[length(slim_positions)]) # total genome length
    return( list(nchr     = nchr,
                 chr_ends = chr_ends,
                 L        = L,
                 msprime_r_map = list(rates     = msprime_rates,
                                      positions = msprime_positions),
                 slim_r_map    = list(rates     = slim_rates,
                                      positions = slim_positions) ))
  }else{
    stop(paste("Missing columns in input file:",header_check$missing))
  }
}
