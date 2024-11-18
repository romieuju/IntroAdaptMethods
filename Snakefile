#    IntroAdapt: Performance test of Adaptive Introgression methods
#    Copyright (C) 2024  Jules Romieu and Ghislain Camarata CNRS and UM
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

import configparser
import os 

# check if user specified a config file in command line (normal behaviour except for tests)
try:
  config["options_file"]
except:
  config["options_file"] = 'tests/config_project.ini'
options_file = config["options_file"]

# check if user specified a python path in command line
try:
  config["python_path"]
except:
  config["python_path"] = 'python'
python_path  = config["python_path"]

# READ OPTIONS FILE
options      = configparser.ConfigParser()
options.read(options_file)

# SETTINGS parameters
project          = options.get('Settings','project')
analysis         = options.get('Settings','analysis')
results_dir      = options.get('Settings', 'results_dir')
simul_type       = options.get('Settings', 'simul_type')
volcano_path     = options.get('Inference','volcano_path')
volcano_div      = options.get('Inference','volcano_div')
volcano_nblocks  = options.getint('Inference','volcano_nblocks')
blocks           = range(1,volcano_nblocks+1)
num_of_sims      = options.getint('Settings','num_of_sims')
sims             = range(1,num_of_sims+1)
nb_chr           = options.getint('Statistics','nb_chr')
chromosomes      = range(1,nb_chr+1)
genome_file      = options.get('Settings','genome_file')
seed             = options.getint('Settings','seed')
keep_files       = True
dir_path         = results_dir


#Contol to the presence of config file
if not os.path.exists(results_dir) :
    sys.exit("The results folder with path : "+results_dir+" does not exist")

#Way to add advantageous mutation in forward simulation
add_adv_mutation_style = str(options.get('Model','add_adv_mutation_style'))
slim_split             = options.get('Model','slim_split')

if add_adv_mutation_style == "manually" :
    script_slim = 'scripts/forwsim_manually.slim'
elif add_adv_mutation_style == "mutation_rate" :
    script_slim = 'scripts/forwsim_mutation_rate.slim'

#AI tools parameters :
wind_reftable_type = options.get('Inference','wind_reftable_type')
reftable_type      = options.get('Inference','reftable_type')

#AI classification methods parameters :
volcano_switch        = options.get('Inference','volcano_switch')
maladapt_switch       = options.get('Inference','maladapt_switch')
maladapt_trained_path = options.get('Inference','maladapt_trained_path')
maladapt_features     = [i for i in options.get('Inference','maladapt_features').split()]
genomatnn_switch      = options.get('Inference','genomatnn_switch')
gnn_t_dir             = options.get('Inference','genomatnn_trained_dir')
gnn_models            = [i for i in options.get('Inference','genomatnn_trained_mod').split()]

#Sample number by pop names (used for VolcanoFinder)
rec_pop   = options.get('Model','recipient')
n_samples = [int(i) for i in options.get('Model',"n_samples").split()]
if rec_pop == '"A1"' :
    n_rec=n_samples[0]
elif rec_pop == '"A2"' :
    n_rec=n_samples[1]
elif rec_pop == '"B1"' :
    n_rec=n_samples[2]
elif rec_pop == '"B2"' :
    n_rec=n_samples[3]
n_haplo_vol=2*n_rec


# evaluate our inference methods
rule end_evaluation :
    input :
        expand('{d}/{a}/{p}/comparison/corr_vol_Q95_{p}.png', d = dir_path, a = analysis, p = project)+
        expand('{d}/{a}/{p}/comparison/dist_vol_rel_{p}.png', d = dir_path, a = analysis, p = project)+
        expand('{d}/{a}/{p}/comparison/fpr_vol_rel_{p}.png', d = dir_path, a = analysis, p = project)+
        expand('{d}/{a}/{p}/comparison/fdr_vol_rel_{p}.png', d = dir_path, a = analysis, p = project) if volcano_switch == "On" else [],
        expand('{d}/{a}/{p}/comparison/corr_mal_Q95_{p}.png', d = dir_path, a = analysis, p = project)+
        expand('{d}/{a}/{p}/comparison/dist_mal_{f}_{p}.png', d = dir_path, a = analysis, p = project, f=maladapt_features)+
        expand('{d}/{a}/{p}/comparison/fpr_mal_{f}_{p}.png', d = dir_path, a = analysis, p = project, f=maladapt_features)+
        expand('{d}/{a}/{p}/comparison/fdr_mal_{f}_{p}.png', d = dir_path, a = analysis, p = project, f=maladapt_features) if maladapt_switch == "On" else [],
        expand('{d}/{a}/{p}/comparison/corr_gnn_Q95_{p}.png', d = dir_path, a = analysis, p = project)+
        expand('{d}/{a}/{p}/comparison/dist_gnn_{id}_{p}.png', d = dir_path, a = analysis, p = project, id=range(1,len(gnn_models)+1))+
        expand('{d}/{a}/{p}/comparison/fpr_gnn_{id}_{p}.png', d = dir_path, a = analysis, p = project, id=range(1,len(gnn_models)+1))+
        expand('{d}/{a}/{p}/comparison/fdr_gnn_{id}_{p}.png', d = dir_path, a = analysis, p = project, id=range(1,len(gnn_models)+1)) if genomatnn_switch == "On" else [],
        expand('{d}/{a}/{p}/comparison/dist_Q_1_100_q95_{p}.png', d = dir_path, a = analysis, p = project),
        expand('{d}/{a}/{p}/comparison/fdr_Q_1_100_q95_{p}.png', d = dir_path, a = analysis, p = project),
        expand('{d}/{a}/{p}/comparison/fpr_Q_1_100_q95_{p}.png', d = dir_path, a = analysis, p = project),
        expand('{d}/{a}/{p}/comparison/roc_{p}.png', d = dir_path, a = analysis, p = project)+
        expand('{d}/{a}/{p}/comparison/roc_leg_{p}.png', d = dir_path, a = analysis, p = project)+
        expand('{d}/{a}/{p}/comparison/precrec_{p}.png', d = dir_path, a = analysis, p = project)+
        expand('{d}/{a}/{p}/comparison/precrec_leg_{p}.png', d = dir_path, a = analysis, p = project)+
        expand('{d}/{a}/{p}/comparison/mccf1_{p}.png', d = dir_path, a = analysis, p = project)+
        expand('{d}/{a}/{p}/comparison/mccf1_leg_{p}.png', d = dir_path, a = analysis, p = project)+
        expand('{d}/{a}/{p}/comparison/zoomed_mccf1_{p}.png', d = dir_path, a = analysis, p = project)+
        expand('{d}/{a}/{p}/comparison/zoomed_mccf1_leg_{p}.png', d = dir_path, a = analysis, p = project)+
        expand('{d}/{a}/{p}/comparison/freq_{p}.png', d = dir_path, a = analysis, p = project)+
        expand('{d}/{a}/{p}/comparison/introprop_{p}.png', d = dir_path, a = analysis, p = project)+
        expand('{d}/{a}/{p}/comparison/time_{p}.png', d = dir_path, a = analysis, p = project) if simul_type == "AI" else [],
        expand('{d}/{a}/{p}/comparison/justin_case_{p}.csv', d = dir_path, a = analysis, p = project),
        expand('{d}/{a}/{p}/comparison/performance_metrics_{p}.csv', d = dir_path, a = analysis, p = project),
        expand('{d}/{a}/{p}/comparison/predictions_{p}.csv', d = dir_path, a = analysis, p = project),

rule evaluation :
    input :
        '{d}/{a}/{p}/maladapt_out_{p}.csv' if maladapt_switch == "On" else '{d}/{a}/{p}/reftable_wind_{a}_{p}.csv',
        expand('{{d}}/{{a}}/{{p}}/{{a}}_{{p}}_sim_{s}/VolcanoFinder/volcanotest_{s}_chr_{c}.out', s = sims, c=chromosomes) if volcano_switch == "On" else [],
        expand('{{d}}/{{a}}/{{p}}/{{a}}_{{p}}_sim_{s}/genomatnn/{mod}/predictions.txt', s = sims, mod=gnn_models) if genomatnn_switch == "On" else [],
        script      = 'scripts/ai_methods_eval.R',
        project_ini = '{d}/{a}/{p}/project_options.ini',
    output :
        ['{d}/{a}/{p}/comparison/corr_vol_Q95_{p}.png',
        '{d}/{a}/{p}/comparison/fpr_vol_rel_{p}.png',
        '{d}/{a}/{p}/comparison/fdr_vol_rel_{p}.png',
        '{d}/{a}/{p}/comparison/dist_vol_rel_{p}.png'] if volcano_switch == "On" else [],
        ['{d}/{a}/{p}/comparison/corr_mal_Q95_{p}.png']+
        [['{d}/{a}/{p}/comparison/dist_mal_'+ mal_feat +'_{p}.png','{d}/{a}/{p}/comparison/fpr_mal_'+ mal_feat +'_{p}.png','{d}/{a}/{p}/comparison/fdr_mal_'+ mal_feat +'_{p}.png'] for mal_feat in maladapt_features] if maladapt_switch == "On" else [],
        ['{d}/{a}/{p}/comparison/corr_gnn_Q95_{p}.png']+
        [['{d}/{a}/{p}/comparison/dist_gnn_'+ str(mod_id) +'_{p}.png','{d}/{a}/{p}/comparison/fpr_gnn_'+ str(mod_id) +'_{p}.png','{d}/{a}/{p}/comparison/fdr_gnn_'+ str(mod_id) +'_{p}.png'] for mod_id in range(1,len(gnn_models)+1)] if genomatnn_switch == "On" else [],
        '{d}/{a}/{p}/comparison/dist_Q_1_100_q95_{p}.png',
        '{d}/{a}/{p}/comparison/fpr_Q_1_100_q95_{p}.png',
        '{d}/{a}/{p}/comparison/fdr_Q_1_100_q95_{p}.png',
        ['{d}/{a}/{p}/comparison/roc_{p}.png',
        '{d}/{a}/{p}/comparison/roc_leg_{p}.png',
        '{d}/{a}/{p}/comparison/precrec_{p}.png',
        '{d}/{a}/{p}/comparison/precrec_leg_{p}.png',
        '{d}/{a}/{p}/comparison/mccf1_{p}.png',
        '{d}/{a}/{p}/comparison/mccf1_leg_{p}.png',
        '{d}/{a}/{p}/comparison/zoomed_mccf1_{p}.png',
        '{d}/{a}/{p}/comparison/zoomed_mccf1_leg_{p}.png',
        '{d}/{a}/{p}/comparison/freq_{p}.png',
        '{d}/{a}/{p}/comparison/introprop_{p}.png',
        '{d}/{a}/{p}/comparison/time_{p}.png'] if simul_type == "AI" else [],
        '{d}/{a}/{p}/comparison/justin_case_{p}.csv',
        '{d}/{a}/{p}/comparison/performance_metrics_{p}.csv',
        '{d}/{a}/{p}/comparison/predictions_{p}.csv',
    shell :
        'Rscript {input.script} {input.project_ini}'


# Run MaLAdapt
rule early_end_maladapt :
    input : 
        expand('{d}/{a}/{p}/maladapt_out_{p}.csv', d = dir_path, a = analysis, p = project),
rule maladapt :
    input:
        script      = 'scripts/run_maladapt.py',
        project_ini = '{d}/{a}/{p}/project_options.ini',
        maladapt_ref='{d}/{a}/{p}/maladapt_tab_{p}.csv',
        trained_model=[maladapt_trained_path+"/"+feat+"_model.sav" for feat in maladapt_features],
        feature_lists=["config_file/maladapt_feature_lists/"+feat+".txt" for feat in maladapt_features],
    output:
        '{d}/{a}/{p}/maladapt_out_{p}.csv',
    resources:
        runtime_min = 30
    conda:
        "maladapt"
    shell:
        '{python_path} {input.script} {input.project_ini}'

rule early_end_maladapt_input:
    input : 
        expand('{d}/{a}/{p}/maladapt_tab_{p}.csv', d = dir_path, a = analysis, p = project),
rule maladapt_input:
    input :
        script      = 'scripts/maladapt_input.py',
        project_ini = '{d}/{a}/{p}/project_options.ini',
        windowed_table='{d}/{a}/{p}/reftable_wind_{a}_{p}.csv',
    output :
        temp('{d}/{a}/{p}/maladapt_tab_{p}.csv')
    resources:
        runtime_min = 30
    shell:
        'time {python_path} {input.script} {input.project_ini}'



# Run VolcanoFinder
rule early_end_merge_outputs :
    input : 
        expand('{d}/{a}/{p}/{a}_{p}_sim_{s}/VolcanoFinder/volcanotest_{s}_chr_{c}.out', d = dir_path, a = analysis, p = project, s = sims, c=chromosomes),
rule merge_outputs :
    input : 
        expand('{{d}}/{{a}}/{{p}}/{{a}}_{{p}}_sim_{{s}}/VolcanoFinder/volcanotest_{{s}}_chr_{{c}}.out_{b}_'+str(volcano_nblocks), b=blocks),
    output :
        '{d}/{a}/{p}/{a}_{p}_sim_{s}/VolcanoFinder/volcanotest_{s}_chr_{c}.out',
    shell :
        '{volcano_path} -m {output} {volcano_nblocks}'

rule early_end_volcanofinder :
    input : 
        expand('{d}/{a}/{p}/{a}_{p}_sim_{s}/VolcanoFinder/volcanotest_{s}_chr_{c}.out_{b}_'+str(volcano_nblocks), d = dir_path, a = analysis, p = project, s = sims, c=chromosomes, b=blocks),
rule volcanofinder :
    input:
        AFF='{d}/{a}/{p}/{a}_{p}_sim_{s}/VolcanoFinder/AFF_volcano_{s}_chr_{c}.txt',
        SFS='{d}/{a}/{p}/{a}_{p}_sim_{s}/VolcanoFinder/SFS_volcano_{s}_chr_{c}.txt',
        dgrid='{d}/{a}/{p}/{a}_{p}_sim_{s}/VolcanoFinder/volcanolook_{s}_chr_{c}_dvalues',
        lookup='{d}/{a}/{p}/{a}_{p}_sim_{s}/VolcanoFinder/volcanolook_{s}_chr_{c}_lookuptable',
    output:
        temp('{d}/{a}/{p}/{a}_{p}_sim_{s}/VolcanoFinder/volcanotest_{s}_chr_{c}.out_{b}_'+str(volcano_nblocks)),
    params :
        prefix='{d}/{a}/{p}/{a}_{p}_sim_{s}/VolcanoFinder/volcanotest_{s}_chr_{c}.out',
        LookupPrefix='{d}/{a}/{p}/{a}_{p}_sim_{s}/VolcanoFinder/volcanolook_{s}_chr_{c}',
    resources:
        runtime_min = 30
    shell:#./VolcanoFinder -pbig g FreqFile SpectFile LookupPrefix OutFile BLOCK NBLOCK
        'time {volcano_path} -pbig 1000 {input.AFF} {input.SFS} {params.LookupPrefix} {params.prefix} {wildcards.b} {volcano_nblocks}'#pas imposés de 1000 entre sites testés, avec divergence donnée par fichier de configuration (mettre valeur négative pour grille par défaut), polarisation et modèle 1

rule early_end_pb_lookup_table :
    input :
        expand('{d}/{a}/{p}/{a}_{p}_sim_{s}/VolcanoFinder/volcanolook_{s}_chr_{c}_dvalues', d = dir_path, a = analysis, p = project, s = sims, c=chromosomes),
        expand('{d}/{a}/{p}/{a}_{p}_sim_{s}/VolcanoFinder/volcanolook_{s}_chr_{c}_lookuptable', d = dir_path, a = analysis, p = project, s = sims, c=chromosomes),
rule pb_lookup_table :
    input:
        SFS='{d}/{a}/{p}/{a}_{p}_sim_{s}/VolcanoFinder/SFS_volcano_{s}_chr_{c}.txt',
    output:
        '{d}/{a}/{p}/{a}_{p}_sim_{s}/VolcanoFinder/volcanolook_{s}_chr_{c}_dvalues',
        '{d}/{a}/{p}/{a}_{p}_sim_{s}/VolcanoFinder/volcanolook_{s}_chr_{c}_lookuptable',
    params :
        prefix='{d}/{a}/{p}/{a}_{p}_sim_{s}/VolcanoFinder/volcanolook_{s}_chr_{c}',
    resources:
        runtime_min = 30
    shell:#./VolcanoFinder –p SpectFile D P MODEL nmin nmax xmin xmax LookupPrefix
        'time {volcano_path} -p {input.SFS} {volcano_div} 1 1 {n_haplo_vol} {n_haplo_vol} 1 {n_haplo_vol} {params.prefix}'#meant to speedup evaluation

# Run genomatnn
rule early_end_genomatnn :
    input : 
        expand('{d}/{a}/{p}/{a}_{p}_sim_{s}/genomatnn/{mod}/predictions.pdf', d = dir_path, a = analysis, p = project, s = sims, mod=gnn_models),
        expand('{d}/{a}/{p}/{a}_{p}_sim_{s}/genomatnn/{mod}/predictions.txt', d = dir_path, a = analysis, p = project, s = sims, mod=gnn_models),
        expand('{d}/{a}/{p}/{a}_{p}_sim_{s}/genomatnn/{mod}.hdf5', d = dir_path, a = analysis, p = project, s = sims, mod=gnn_models),
        expand('{d}/{a}/{p}/{a}_{p}_sim_{s}/genomatnn/don-rec-sis_{s}.vcf.gz', d = dir_path, a = analysis, p = project, s = sims),
        expand('{d}/{a}/{p}/{a}_{p}_sim_{s}/genomatnn/don-rec-sis_{s}.vcf.gz.tbi', d = dir_path, a = analysis, p = project, s = sims),
rule genomatnn :
    input:
        don_id      ='{d}/{a}/{p}/{a}_{p}_sim_{s}/genomatnn/donnor_{s}.indlist',
        sis_id      ='{d}/{a}/{p}/{a}_{p}_sim_{s}/genomatnn/sister_{s}.indlist',
        rec_id      ='{d}/{a}/{p}/{a}_{p}_sim_{s}/genomatnn/recipient_{s}.indlist',
        mod_ori     = gnn_t_dir + '{mod}.hdf5',
        gnn_config  ='{d}/{a}/{p}/{a}_{p}_sim_{s}/genomatnn/config_{s}.toml',
        zipped_vcf  ='{d}/{a}/{p}/{a}_{p}_sim_{s}/genomatnn/don-rec-sis_{s}.vcf.gz',
        vcf_index   ='{d}/{a}/{p}/{a}_{p}_sim_{s}/genomatnn/don-rec-sis_{s}.vcf.gz.tbi',
    output:
        mod_copy    =temp('{d}/{a}/{p}/{a}_{p}_sim_{s}/genomatnn/{mod}.hdf5'),
        gnn_graph   ='{d}/{a}/{p}/{a}_{p}_sim_{s}/genomatnn/{mod}/predictions.pdf',
        gnn_score   ='{d}/{a}/{p}/{a}_{p}_sim_{s}/genomatnn/{mod}/predictions.txt',
    resources:
        runtime_min = 30
    conda:
        "genomatnn"
    shell:
        """
        cp {input.mod_ori} {output.mod_copy}
        time genomatnn apply {input.gnn_config} {output.mod_copy}
        """

#vcf pre genomatnn adjustment
rule early_end_genomatnn_vcf :
    input :
        expand('{d}/{a}/{p}/{a}_{p}_sim_{s}/genomatnn/don-rec-sis_{s}.vcf.gz', d = dir_path, a = analysis, p = project, s = sims),
        expand('{d}/{a}/{p}/{a}_{p}_sim_{s}/genomatnn/don-rec-sis_{s}.vcf.gz.tbi', d = dir_path, a = analysis, p = project, s = sims),
rule genomatnn_vcf :
    input:
        gnn_vcf     ='{d}/{a}/{p}/{a}_{p}_sim_{s}/genomatnn/don-rec-sis_{s}.vcf',
    output:
        zipped_vcf  =temp('{d}/{a}/{p}/{a}_{p}_sim_{s}/genomatnn/don-rec-sis_{s}.vcf.gz'),
        vcf_index   =temp('{d}/{a}/{p}/{a}_{p}_sim_{s}/genomatnn/don-rec-sis_{s}.vcf.gz.tbi'),
    resources:
        runtime_min = 30
    conda:
        "genomatnn"
    shell:
        """ 
        bgzip {input.gnn_vcf}
        tabix {output.zipped_vcf}
        """

# creation of input files for genomatnn and volcanofinder
rule early_end_input_files :
    input : 
        [expand('{d}/{a}/{p}/{a}_{p}_sim_{s}/VolcanoFinder/SFS_volcano_{s}_chr_{c}.txt', d = dir_path, a = analysis, p = project, s = sims, c=chromosomes),
        expand('{d}/{a}/{p}/{a}_{p}_sim_{s}/VolcanoFinder/AFF_volcano_{s}_chr_{c}.txt', d = dir_path, a = analysis, p = project, s = sims, c=chromosomes)] if volcano_switch == "On" else [],
        [expand('{d}/{a}/{p}/{a}_{p}_sim_{s}/genomatnn/donnor_{s}.indlist', d = dir_path, a = analysis, p = project, s = sims, c=chromosomes),
        expand('{d}/{a}/{p}/{a}_{p}_sim_{s}/genomatnn/sister_{s}.indlist', d = dir_path, a = analysis, p = project, s = sims, c=chromosomes),
        expand('{d}/{a}/{p}/{a}_{p}_sim_{s}/genomatnn/recipient_{s}.indlist', d = dir_path, a = analysis, p = project, s = sims, c=chromosomes),
        expand('{d}/{a}/{p}/{a}_{p}_sim_{s}/genomatnn/config_{s}.toml', d = dir_path, a = analysis, p = project, s = sims, c=chromosomes),
        expand('{d}/{a}/{p}/{a}_{p}_sim_{s}/genomatnn/don-rec-sis_{s}.vcf', d = dir_path, a = analysis, p = project, s = sims, c=chromosomes)] if genomatnn_switch == "On" else [],

rule input_files :
    input:
        script       = 'scripts/tree_seq_to_input.py',
        sim_ini      = '{d}/{a}/{p}/{a}_{p}_sim_{s}/{a}_{p}_sim_{s}.ini',
        project_ini  = '{d}/{a}/{p}/project_options.ini',
        mutsim_trees = '{d}/{a}/{p}/{a}_{p}_sim_{s}/Mut_TreeSeq_{a}_{p}_sim_{s}.trees',
        awk_script   = 'scripts/parse_gnn.awk',
        gnn_template = 'config_file/gnn_config.toml',
    output:
        ['{d}/{a}/{p}/{a}_{p}_sim_{s}/VolcanoFinder/SFS_volcano_{s}_chr_'+str(ch_num)+'.txt' for ch_num in chromosomes]+
        ['{d}/{a}/{p}/{a}_{p}_sim_{s}/VolcanoFinder/AFF_volcano_{s}_chr_'+str(ch_num)+'.txt' for ch_num in chromosomes] if volcano_switch == "On" else [],
        ['{d}/{a}/{p}/{a}_{p}_sim_{s}/genomatnn/donnor_{s}.indlist',
        '{d}/{a}/{p}/{a}_{p}_sim_{s}/genomatnn/sister_{s}.indlist',
        '{d}/{a}/{p}/{a}_{p}_sim_{s}/genomatnn/recipient_{s}.indlist',
        '{d}/{a}/{p}/{a}_{p}_sim_{s}/genomatnn/config_{s}.toml',
        '{d}/{a}/{p}/{a}_{p}_sim_{s}/genomatnn/don-rec-sis_{s}.vcf'] if genomatnn_switch == "On" else [],
    resources:
        runtime_min = 30
    conda:
        "introadapt_end"
    shell:
        'time {python_path} {input.script} {input.project_ini} {input.sim_ini}'


# Create merge dataframe (parameters, latent var and sum stat) by windows and for each simulated data :
rule early_merge_dataframe :
    input : 
        reftable_file        = [expand('{d}/{a}/{p}/reftable_wind_{a}_{p}.csv', d = dir_path, a = analysis, p = project) if reftable_type == "classification" or reftable_type == "all" else [],
                                expand('{d}/{a}/{p}/reftable_{a}_{p}.csv', d = dir_path, a = analysis, p = project) if reftable_type == "estimation" or reftable_type == "all" else []],
        reftable_file_by_sim = [expand('{d}/{a}/{p}/reftable_{a}_{p}_sim_'+str(sim_number)+'.csv', d = dir_path, a = analysis, p = project) for sim_number in sims] if reftable_type == "estimation_by_sim" or reftable_type == "all" else [],

rule merge_dataframe:
    input:
        script                = 'scripts/merge_dataframe.R',
        options_file          = '{d}/{a}/{p}/project_options.ini',
        sim_summary_stat_wind = expand('{d}/{a}/{p}/{a}_{p}_sim_{s}/Sum_Stat_Mut_TreeSeq_{a}_{p}_sim_{s}.csv', d = dir_path, a = analysis, p = project, s = sims),
        sim_summary_stat_all  = expand('{d}/{a}/{p}/{a}_{p}_sim_{s}/sum_stat_{a}_{p}_sim_{s}.csv', d = dir_path, a = analysis, p = project, s = sims),
        ai_freq_rec           = expand('{d}/{a}/{p}/{a}_{p}_sim_{s}/forwsim_freq_mut_don_in_rec_{a}_{p}_sim_{s}.txt', d = dir_path ,a = analysis ,p = project ,s = sims),
        latent_var            = expand('{d}/{a}/{p}/{a}_{p}_sim_{s}/forwsim_latent_variable_{a}_{p}_sim_{s}.txt', d = dir_path ,a = analysis ,p = project ,s = sims),
        data_params           = expand('{{d}}/{{a}}/{{p}}/{{a}}_{{p}}_sim_{s}/data_parameters_{{a}}_{{p}}_sim_{s}.csv',d = dir_path, a = analysis, p = project, s = sims)
    output:
       reftable_file        = [['{d}/{a}/{p}/reftable_wind_{a}_{p}.csv'] if reftable_type == "classification" or reftable_type == "all" else [],
                                 ['{d}/{a}/{p}/reftable_{a}_{p}.csv'] if reftable_type == "estimation" or reftable_type == "all" else []],
       reftable_file_by_sim = [['{d}/{a}/{p}/reftable_{a}_{p}_sim_'+str(sim_number)+'.csv' for sim_number in sims] if reftable_type == "estimation_by_sim" or reftable_type == "all" else []],
    resources:
        runtime_min = 10
    conda:
        "introadapt_end"
    shell:
        'Rscript {input.script} {input.options_file}'

# create summary statistics table
rule early_end_summary_stats :
    input : 
        sim_summary_stat_wind = expand('{d}/{a}/{p}/{a}_{p}_sim_{s}/Sum_Stat_Mut_TreeSeq_{a}_{p}_sim_{s}.csv', d = dir_path, a = analysis, p = project, s = sims),
        sim_summary_stat_all  = expand('{d}/{a}/{p}/{a}_{p}_sim_{s}/sum_stat_{a}_{p}_sim_{s}.csv', d = dir_path, a = analysis, p = project, s = sims)
rule summary_stats:
    input:
        script        = 'scripts/sum_stat_computation.py',
        sim_ini       = '{d}/{a}/{p}/{a}_{p}_sim_{s}/{a}_{p}_sim_{s}.ini',
        project_ini   = '{d}/{a}/{p}/project_options.ini',
        mutsim_trees  = '{d}/{a}/{p}/{a}_{p}_sim_{s}/Mut_TreeSeq_{a}_{p}_sim_{s}.trees'
    output:
        sim_summary_stat_wind = '{d}/{a}/{p}/{a}_{p}_sim_{s}/Sum_Stat_Mut_TreeSeq_{a}_{p}_sim_{s}.csv',
        sim_summary_stat_all  = '{d}/{a}/{p}/{a}_{p}_sim_{s}/sum_stat_{a}_{p}_sim_{s}.csv'
    resources:
        runtime_min = 30
    conda:
        "introadapt_end"
    shell:
        '{python_path} {input.script} {input.project_ini} {input.sim_ini}'


# simulation of mutations with msprime (add neutral simulation)
rule early_end_mutation_sim :
    input : 
        mut_treeseq     = expand('{d}/{a}/{p}/{a}_{p}_sim_{s}/Mut_TreeSeq_{a}_{p}_sim_{s}.trees', d = dir_path, a = analysis, p = project, s = sims),
        i_interv_recomb = expand('{d}/{a}/{p}/{a}_{p}_sim_{s}/Recomb_interval_intro_{a}_{p}_sim_{s}.csv', d = dir_path, a = analysis, p = project, s = sims)

rule mutation_sim:
    input:
        script        = 'scripts/mutsim.py',
        sim_ini       = '{d}/{a}/{p}/{a}_{p}_sim_{s}/{a}_{p}_sim_{s}.ini',
        project_ini   = '{d}/{a}/{p}/project_options.ini',
        forwsim_trees = '{d}/{a}/{p}/{a}_{p}_sim_{s}/forwsim_final_{a}_{p}_sim_{s}.trees'
    output:
        mut_treeseq     = '{d}/{a}/{p}/{a}_{p}_sim_{s}/Mut_TreeSeq_{a}_{p}_sim_{s}.trees',
        i_interv_recomb = '{d}/{a}/{p}/{a}_{p}_sim_{s}/Recomb_interval_intro_{a}_{p}_sim_{s}.csv'
    resources:
        runtime_min = 30
    conda:
        "introadapt_end"
    shell:
        '{python_path} {input.script} {input.project_ini} {input.sim_ini}'

# simulation with SLiM (forward simulation + variable latent files)
rule early_end_forward_sim:
    input:
        forwsim_trees      = expand('{d}/{a}/{p}/{a}_{p}_sim_{s}/forwsim_final_{a}_{p}_sim_{s}.trees', d = dir_path ,a = analysis ,p = project ,s = sims),
        ai_freq_rec        = expand('{d}/{a}/{p}/{a}_{p}_sim_{s}/forwsim_freq_mut_don_in_rec_{a}_{p}_sim_{s}.txt', d = dir_path ,a = analysis ,p = project ,s = sims),
        ai_freq_rec_samp   = expand('{d}/{a}/{p}/{a}_{p}_sim_{s}/forwsim_freq_mut_don_in_rec_samp_{a}_{p}_sim_{s}.txt', d = dir_path ,a = analysis ,p = project ,s = sims),
        latent_var         = expand('{d}/{a}/{p}/{a}_{p}_sim_{s}/forwsim_latent_variable_{a}_{p}_sim_{s}.txt', d = dir_path ,a = analysis ,p = project ,s = sims)

rule forward_sim:
    input:
        script        = expand(script_slim),
        slim_options  = '{d}/{a}/{p}/{a}_{p}_sim_{s}/{a}_{p}_sim_{s}.eidos',
        coalsim_trees = '{d}/{a}/{p}/{a}_{p}_sim_{s}/coalsim_{a}_{p}_sim_{s}.trees'
    output:
        forwsim_trees      = '{d}/{a}/{p}/{a}_{p}_sim_{s}/forwsim_final_{a}_{p}_sim_{s}.trees',
        ai_freq_rec        = '{d}/{a}/{p}/{a}_{p}_sim_{s}/forwsim_freq_mut_don_in_rec_{a}_{p}_sim_{s}.txt',
        ai_freq_rec_samp   = '{d}/{a}/{p}/{a}_{p}_sim_{s}/forwsim_freq_mut_don_in_rec_samp_{a}_{p}_sim_{s}.txt',
        latent_var         = '{d}/{a}/{p}/{a}_{p}_sim_{s}/forwsim_latent_variable_{a}_{p}_sim_{s}.txt'
    resources:
        runtime_min = 120
    conda:
        "introadapt_end"
    shell:
        'time slim -l 0 -d "option_file=\'{input.slim_options}\'" {input.script}'

# Backward simulation (msprime)
rule early_end_backward_sim:
    input:
        coalsim_trees = expand('{d}/{a}/{p}/{a}_{p}_sim_{s}/coalsim_{a}_{p}_sim_{s}.trees',d = dir_path ,a = analysis ,p = project, s = sims),

rule backward_sim:
    input:
        script      = 'scripts/coalsim.py',
        sim_ini     = '{d}/{a}/{p}/{a}_{p}_sim_{s}/{a}_{p}_sim_{s}.ini',
        project_ini = '{d}/{a}/{p}/project_options.ini'
    output:
        coalsim_trees = '{d}/{a}/{p}/{a}_{p}_sim_{s}/coalsim_{a}_{p}_sim_{s}.trees'
    resources:
        runtime_min = 120
    conda:
        "introadapt_end"
    shell:
        '{python_path} {input.script} {input.project_ini} {input.sim_ini} graph'


# read parameters and sample from priors (create forward and backward config file and parameters values files)
rule early_end_parameters:
    input:
        slim_options     = expand('{d}/{a}/{p}/{a}_{p}_sim_{s}/{a}_{p}_sim_{s}.eidos', d = dir_path, a = analysis, p = project, s = sims),
        sim_ini          = expand('{d}/{a}/{p}/{a}_{p}_sim_{s}/{a}_{p}_sim_{s}.ini', d = dir_path, a = analysis, p = project, s = sims),
        data_params      = expand('{{d}}/{{a}}/{{p}}/{{a}}_{{p}}_sim_{s}/data_parameters_{{a}}_{{p}}_sim_{s}.csv',d = dir_path, a = analysis, p = project, s = sims)
    resources:
        runtime_min = 10
rule draw_parameters:
    input:
        script      = 'scripts/get_params.R',
        options_ini = '{d}/{a}/{p}/project_options.ini'
    output:
        slim_options     = expand('{{d}}/{{a}}/{{p}}/{{a}}_{{p}}_sim_{s}/{{a}}_{{p}}_sim_{s}.eidos', s = sims),
        sim_ini          = expand('{{d}}/{{a}}/{{p}}/{{a}}_{{p}}_sim_{s}/{{a}}_{{p}}_sim_{s}.ini', s = sims),
        data_params      = expand('{{d}}/{{a}}/{{p}}/{{a}}_{{p}}_sim_{s}/data_parameters_{{a}}_{{p}}_sim_{s}.csv', s = sims),
    resources:
        runtime_min = 10
    conda:
        "introadapt_end"
    shell:
        'Rscript {input.script} {input.options_ini}'


# read project parameters (config file to project file)
rule setup_project:
    input:
        script       = 'scripts/set_project.R',
        options_file = options_file
    output:
       expand('{d}/{a}/{p}/project_options.ini', d = dir_path, a = analysis, p = project),
    resources:
        runtime_min   = 10
    conda:
        "introadapt_end"
    shell:
        'Rscript {input.script} {input.options_file}'