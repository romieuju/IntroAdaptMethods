## IntroAdapt

IntroAdaptMethod is a pipeline thta use simulated genomic data to test the perfromances of 4 adaptive introgression (AI) classification methods: A genomic scan method based on the Q95(1%,100%) summary statistic (Racimo et al., 2017 : https://doi.org/10.1093/molbev/msw216), VolcanoFinder (Setter et al., 2020 : https://doi.org/10.1371/journal.pgen.1008867
), genomatnn (Gower et al., 2021 :  https://doi.org/10.7554/eLife.64669) and MaLAdapt (Zhang et al., 2023 : https://doi.org/10.1093/molbev/msad001). 
If you have any questions about the use and limitations of this pipeline, please contact Jules Romieu (jules.romieu@umontpellier.fr). 

Warning : The pipeline has only been tested under a Linux operating system.

## I/ Installing Conda :
The pipeline requires a Conda environment to be installed, so you will need to install Conda (https://conda.io/projects/conda/en/latest/index.html) to run it.

## II/ Installing the environment :
The introadaptmethods.yml file can be used to install the environment. In the folder containing the yml file, open a command terminal and type :

```shell
conda env create -f introadaptmethods.yml
```

The process of creating the environment can take a few minutes. 
To check that the environment is correctly installed and to activate the Conda environment, type :

```shell
conda activate introadaptmethods
```

## III/Install AI inference methods:

a) Volcanofinder (Setter et al., 2020)
To install Volcanofinder, go to this site:
https://degiorgiogroup.fau.edu/vf.html

b) genomatnn (Gower et al., 2021)
Go to the :
https://github.com/grahamgower/genomatnn
and follow the installation instructions.

c) MaLAdapt (Zhang et al., 2023)
The original scripts can be found on the git:
https://github.com/xzhang-popgen/maladapt
We use a conda environment dedicated to this tool. To install it, open your terminal:

```shell
conda env create -f maladapt.yml -n maladapt
conda activate maladapt
```

The scripts we use to infer AI from simulated genetic data with MaLadapt can be found in the scripts folder: ‘run_maladapt.py’ for a standard way to standardized test data and in the MaLAdapt git for their way of standardized in the original article. The maladapt_input.py script contains the names of the features and summary statistics used. 

Warning 1: The pipeline does not include the ability to train a new ETC model. For this we used the 5trainMaLAdapt_adapted.py script (in the scripts/standelone_scripts folder), a modified version of the script available on their git. This script requires a training table in the MaLAdapt format and a file containing the features to be kept (in the maladapt/feature/ folder).

Warning 2: The pipeline does not produce the graphs shown in the article (ROC curves and violin plot). Examples of the R scripts used to obtain these graphs can be found in the scripts/standelone_scripts folder. They cannot be run directly after the last rule in the snakemake pipeline (end_evaluation) and the paths to the files must be entered manually. 

## IV/Pipeline architecture:

All the information needed for the pipeline to function is contained in an .ini file, which must be given to the terminal to launch the pipeline. The pipeline architecture is constructed as follows:

A/Data simulation part:

1) From the configuration file, create the project + file (setproject.R)

2) Picks the parameter values then creates the configuration files for the backward (msprime) and forward (SLiM) simulators (getparams.R)

3) Simulates the genealogy of the ancestral populations to the populations of interest (donor and recipient) to the ancestral population to all the populations (outgroup + ingroup) in backward (coalsim.py)

4) Using the tree sequence obtained with the backward simulator, simulate the genealogy from the ancestral population to the populations of interest up to the present day. Add the advantageous mutation and the migration event. (.eidos file)

5) Add the neutral mutations to the tree sequence obtained (mutsim.py)

B/Analysis of AI inference methods:

6) Compute summary statistics for each simulated data and each window of each simulated data (For MaLAdapt + Q95)

7) Select the features and summary statistics required for MaLAdapt and create the input files for genomatnn and VolcanoFinder

8) Run method analysis (MaLAdapt, genomatnn and VolcanoFinder)

9) Evaluate methods (for all genome windows, AI, Adjacent and other chromosomes)
We advise you to run part A independently of part B. 

## V/Use of the pipeline:
1) Fill in the .ini configuration file with the parameter values of interest (an example is An example is available in /config_file
/ : config_project_example.ini and genome_example.txt all the ini file use for the manuscrit are save in /config_file
/config_file_article/)

2) Open a terminal in the IntroAdaptMethods folder and type :

```shell
conda activate introadaptmethods
snakemake -s Snakefile rule_to_reach -C options_file=‘config_file_path/config_file_name’ --cores core_number --use-conda
```

Replace :
- rule_to_reach: with the rule at which the pipeline should end (visible in the Snakemake file)
- config_file_path with the path of the folder containing the configuration file
- config_file_name by the name of the configuration file
- core_number by the number of cores dedicated to the pipeline

## VI/The snakemake rules:

- setup_project (script: set_project.R): Adds the information included in the genome file to the project configuration file.
- early_end_parameters (script : get_params.R):
Picks parameter values and creates the initiation files for the Backward and Forward simulator configuration files.
- early_end_backward_sim (script : coalsim.py)
From the configuration file, runs the backward simulations. 
- early_end_forward_sim (script : scripts with .slim)
From the tree sequence obtained in backward, runs the forward simulations 
- early_end_mutation_sim (scipts : mutsim.py)
From the tree sequences obtained in forward add the neutral mutations in forward
- early_end_summary_stats (script : sum_stat_computation.py)
Recover the genotype matrix from the tree sequence and calculate the summary statistics per window. 
- early_end_project_stats (scripts : merge_dataframe.R)
Create tables containing parameter values, latent variables and summary statistics for each simulated data item (reftable) and for each window (reftable_wind).
- early_end_input_files (scripts : tree_seq_to_input.py)
Creates input files for VolcanoFinder (SFS and AFF) and genomatnn (vcf, list of individuals) from genotype matrix
- early_end_genomatnn_vcf:
Compresses vcf files ( bgzip) + indexing (tabix) 
- early_end_genomatnn : 
execute genomatnn 
- early_end_pb_lookup_table :
Creates lookup table for volcanofinder
- early_end_volcanofinder:
Execute volcanofinder (analyses 1 site every 1000 bp)
-early_end_merge_outputs:
Merges volcanofinder result blocks
- early_end_maladapt_input (script : maladapt_imput.py)
Transforms the summary statistics table per window into a table in MaLAdapt format
- early_end_maladapt (script : run_maladapt.py)
Execute MaLAdapt
- end_evaluation (script : ai_method_eval.R)
Compares the performances of the methods (caution: for all windows in the simulated datasets) 

## VII/Config option file features :

| Parameter name | type | description |
|---|---|---------------|
|**[Settings]**|||
| seed        | integer | Seed use for the project and simulation |
| config_file | String | Path and name of this config file |
| result_dir  | String | Path of the resultat directory (structure result_dir/analysis/project) |
| analysis    | String | Name of the analysis (structure result_dir/analysis/project) |
| project     | String | Name of the project (structure result_dir/analysis/project) |
| genome_file | String | Path and name of the genome config file (with informartion about length, chromosome and recombination rate) |
| num_of_sim  | Integer | Number of genetic data to be simulated (integer) |
| simulation_type | String | Simulation type (AI : adaptive introgression or NI : Neutral introgression) |
|**[Model]**||| (In our simulation: O: outgroup, A1 : donor, A2 : donor sister, B1 : recipient sister, B2 : recipient)
| populations | Integer | Population number in present time (5) warning : don't change this |
| popAncO     | Integer or Function | Ancestral population of outgroup and ingroup size |
| popO        | Integer or Function | Outgroup population size |
| popAncAB    | Integer or Function | Ancestral population of donor and recipient population size |
| popA1       | Integer or Function | Donor sister (in our simulation) population size |
| popB1       | Integer or Function | Recipient sister (in our simulation) population size |
| popAncA     | Integer or Function | Ancestral donor and its sister population size |
| popAncB     | Integer or Function | Ancestral recipient and its sister population size |
| popA2       | Integer or Function | Donor (in our simulation) population size |
| popB2       | Integer or Function | Recipient (in our simulation) population size |
| recipient   | String              | Recipient population ("B2" in our simulation) |
| donor       | String              | donor population ("A2" in our simulation) |
| coalescence_split | String | Way to make divergence between populations in backward simulations (split for add_population_split and mass_migration for add_mass_migration) |
| slim_split  | Integer | If divergence between population occur in forward simulation (1) or in backward simulation (0) |
| add_adv_mutation_style | String | To choose how to add advantageous mutations , manually to choose the location of the mutation and mutation_rate so that advantageous mutations appear according to an advantageous mutation rate |
| migration_rate_r       | Float or Function | Migration rate from donor to recipient |
| generation_split_OA | Integer or Function | Generation (in backward) split between outgroup and ancestral ingroup population |
| generation_split_AB | Integer or Function | Generation (in backward) split between ancestral population of the donor and recipient populations |
| generation_split_A  | Integer or Function | Generation (in backward) split between ancestral population of donor and its sister population |
| generation_split_B  | Integer or Function | Generation (in backward) split between ancestral population of recipient and its sister population |
| generations_migration_start | Integer or Function | Generation (in backward) of migration start |
| generations_migration_end | Integer or Function | Generation (in backward) of migration end |
| generations_forward | Integer or Function | Generation  (in backward) of forward (slim) simulation start |
| generation_mutation | Integer or Function | Generation (in backward) of advantageous mutation occurence in ancestral or donor pop |
| sample_time_d_s     | Integer or Function | Sampling generation (in backward) of the donor sister population |
| sample_time_d       | Integer or Function | Sampling generation  (in backward) of the donor population |
| coef_selec_before_mig | Float or Function | Selection coefficient of the advantageous mutation |
| mu_advantageous | Float | Advantageous mutation rate |
| mu_total | Float | Total mutation rate (Neutral mutation rate = total mutation rate - advantageous mutation rate) |
| genome_l | Integer | Length of the genome -1 in base pair |
| r_map_rate | Float or string | Recombination rate or file if the recombination rate is specify in genome.txt file |
| n_samples | List of integer | Number of individual sample by population (donor_sister donor recipient_sister recipient outgroup) |
| scaling_factor | Integer | Value of the scaling factor |
| genomic_region_size | Integer | Size of genome regions (bp) able to host advantageous mutations (if add_adv_mutation_style == mutation_rate) |
| proportion_under_selection | Float or Integer | Proportion of the genome able to host advantageous mutations (if add_adv_mutation_style == manually, this parameters = position of the advantageous
 mutation) | 
| max_region_proportion | Float | Maximum proportion of the genome that can host advantageous mutations (used if max_region_proportion == mutation_rate) |
|**[Statistics]**|||
| zns | String | If the Kelly ZnS summary statistic is to be used in MaLAdapt inferences (Yes or No) |
| donor_pop_stat | String | From which donor population the summary statistics should be used (from the donor population, "A2" or its sister population, "A1") |
| prop_i_by_r_interv | String | If the proportion of introgression and adaptive introgression should be calculated for each recombination interval (Yes or No, warning : don't change this parameters) |
| isphased | String | If genetic data are phased or not and some phased summary stat can be calculated (Yes or No) |
| window_size | Integer | Windows size |
| window_start | Integer | Position in the genome from which the definition of the windows must be defined (position in base pairs) |
| window_end | Integer | Position in the genome where the windows must end (position in base pairs) |
| window_step | Integer | For sliding windows, step between the windows |
| nb_chr | Integer | Chromosome number in the genome |
|**[Inference]**|||
| volcano_switch | String | If VolcanoFinder have to be use in the analysis (On/off) |
| volcano_path | String | Path of VolcanoFinder folder |
| volcano_div | Integer or file | If divergence between donor and recipient have to be estimate (-1) or is given by the user (file with different values of divergence) |
| volcano_nblocks | Integer | Number of blocks into which the genome must be divided for VolcanoFinder analysis |
| maladapt_switch | String | If MaLadapt have to be use in the analysis (On/off) |
| maladapt_trained_path | String | Path to MaLadapt's trained models |
| maladapt_features | String or list of string | Name of MaLAdapt's trained models to be used |
| genomatnn_switch | String | If genomatnn have to be use in the analysis (On/off) |
| genomatnn_trained_dir | String | Path to genomatnn's trained models |
| genomatnn_trained_mod | String or list of string | Name of trained models to be used for genomatnn |
| volcano_threshold | Float | Score threshold to be used to classify windows as AI or non-AI with VolcanoFinder |
| genomatnn_threshold | Float | Score threshold to be used to classify windows as AI or non-AI with genomatnn |
| maladapt_threshold | Float | Score threshold to be used to classify windows as AI or non-AI with MaLadapt |
| q95_threshold | Float | Score threshold to be used to classify windows as AI or non-AI with Q95(1,100) |
| wind_reftable_type | String | If the merge table (parameters, latent var and sum stat values) is a reftable (model_train) or a test data (classifier) (warning : may be not work with the actual version of the pipeline)
| reftable_type | String | Type of dataframe to obtain after the merging of parameters, latent var and sum stat dataframe, windows dataframe for all genetic data (classification), genome dataframe for all genetic data (estimation) or both (all), for performance ai classification method use all or classification |

## VIII/Files after pipeline execution :
The pipeline will create an analysis folder in the result_dir folder, containing a project folder (result_dir/analysis/project). The project folder will contain all the files and folders for the different replicates of simulated genomic data. 

In results/analysis/project/ folder :
analysis_project_sim_n/  : folder with files for each simulated data replicat
- project_options.ini : config file of the project
- reftable_wind_analysis_project.csv : Merge dataframe with parameters, latent varibale and summary statistics for all windows of all genetic data replicat.

In analysis_project_sim_n folder (sim_n = simulated genetic data number n)
- analysis_project_sim_n.eidos                                : Forward (SLiM) config file 
- analysis_project_sim_n.ini                                  : Backward (msprime) config file 
- data_parameters_analysis_project_sim_n.csv                  : Parameters values with scaling factor (from get_params.R script)
- forwsim_final_analysis_project_sim_n.trees                  : Tree sequence at the end of forward (SLiM) simulation (from slim script)
- forwsim_freq_mut_don_in_rec_analysis_project_sim_n.txt      : AI mutation informations in all recipient population (from slim script)
- forwsim_freq_mut_don_in_rec_samp_analysis_project_sim_n.txt : AI mutation informations in sampled recipient individuals (from slim script)
- forwsim_latent_variable_analysis_project_sim_n.txt          : AI Latente variables from slim script and mutsim.py script
- forwsim_start_analysis_project_sim_n.trees                  : Tree sequence at the first generation of slim script (use when AI mutation is lost and forward simulation restart).
- Mut_TreeSeq_analysis_project_sim_n.trees                    : Tree sequence after add neutral mutation and simplifying with genealogy of sampled individual (from mutsim.py)
- Recomb_interval_intro_analysis_project_sim_n.csv            : Information about introgression for each recombination interval in recipient population (from mutsim.py)   
- sum_stat_analysis_project_sim_n.csv                         : Summary statistics at genomic scale of all windows summary statistics 
- Sum_Stat_Mut_TreeSeq_analysis_project_sim_n.csv             : Windows summary statistics 

In analysis_project_sim_n/genomatnn folder :
- config_n.toml       : genomatnn config file (see genomatnn git)
- donor_n.indlist     : donnor samples names list (in vcf)
- recipient_n.indlist : recipient samples names (in vcf)
- sister_n.indlist    : recipient sister samples names (in vcf)
- don-rec_sis_n.vcf   : vcf with variant information for donor, recipient and recipient sister samples (tempory file)

In analysis_project_sim_n/genomatnn/trained_CNN_used_name
- predictions.pdf : AI probability by windows graph
- prediction.txt  : AI probability by window values

analysis_project_sim_n/VolconaFinder
- AFF_volcano_n_chr_chronumber.txt : AFF file for replicat number n and chromosome specifiying by chronumber in volcanofinder file
- SFS_volcano_n_chr_chronumber.txt : SFS file of all the chromosome or the genome. 
- volcanotest_n_chr_chronumber.out : Volcanofinder result for the replicat number n and the chromosome number chronumber (test site : position, likelihood ratio, alpha and D)
- volcavotest_n_chr_chronumber_dvalues : D values test to estimates LR value 

In comparison folder : Folder containing the output files of method performance comparisons for all the test genetic data available in the project. 
- performance_metric_project.csv                               : Classification metrics and statistics for a priori thresold (define in .ini with the method_threshold parameters), mccf1 thresold  and FPR< or = 0.05 thresold. 
- Prediction_project.csv                                       : Method score value (prediction column) by genetic data and by window (define by sim, start and end column), with window true class type (AI=1/non-AI=0), predicted class type (for a priori, mccf and fpr<=0.05 thresold) and some latent variable (AI mut freq in rec, AI mut fixation time in rec and MaLAdapt introgression proportion)
- Classification curves : FDR by score by method, FPR by score by method, density score by class (AI = red and non-AI = blue), mccf1 (
https://doi.org/10.48550/arXiv.2006.11278), ROC, Precision-recall.



Warning : In the pipeline, performance tests are carried out using all the windows in the folder. For example, if a project contains 200 simulations with AI, the genome is made up of 1 chromosome of 1Mb with a mutation under AI and the non-overlapping windows are 50kb long, then the performance tests will be carried out on the 4000 windows, including 200 windows under AI and 3800 non-AI. If the genome is made up of 2 chromosomes, the performance tests will be carried out by taking into account the windows of the first chromosome and the second without differentiating between them. If the user wishes to calculate classification metrics for a test dataset containing a certain type of non-AI window (Adjacente or neutral chromosome for example). They can use the method score values stored in the Prediction_project.csv file and keep the sim, start, end, method, classifier and prediction columns for the windows they are interested in, and then calculate their own classification metrics. 