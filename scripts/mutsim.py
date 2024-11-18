#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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

#This script will take tree sequence from slim simulation, add neutral mutation with msprime and create a dataframe with recombinaison intervall and the propotion of introgression of each interval. 

#Package use in this python script
import sys
import msprime
import numpy as np
import pandas as pd
import introadapt
import tskit
import time
import statistics

def main():

  mutsim_start = time.time()

  # Ignoring some warnings when calculating summary stats with missing data,
  np.seterr(invalid='ignore',divide='ignore')

  # Get options for project and simulation:
  options = introadapt.get_options(proj_options_file = sys.argv[1], sim_options_file = sys.argv[2])

  #Save the results_dir path :
  project_dir = options["results_dir"]+"/"+options["analysis"]+"/"+options["project"]
  job_dir = project_dir+"/"+options["analysis"]+"_"+options["project"]+"_sim_"+options["sim"]

  #Define parameters :
  donor = options["donor"]
  recipient =  options["recipient"]
  g_end_sim = options["generations_forward"]
  g_samp = g_end_sim - (g_end_sim-1)
  sample_time_do = options["sample_time_d"]
  sample_time_do_s = options["sample_time_d_s"]
  g_mig = options["generations_migration_start"]
  g_between_end_and_mig = g_end_sim - g_mig
  g_samp_don = g_end_sim - sample_time_do
  g_samp_don_s = g_end_sim - sample_time_do_s
  g_mig = (g_end_sim) - options["generations_migration_start"]
  geno_len = options["L"]+1
  prop_i_by_r_interv = options["prop_i_by_r_interv"]

  n_popA1=options["n_samples"][0]
  n_popA2=options["n_samples"][1]
  n_popB1=options["n_samples"][2]
  n_popB2=options["n_samples"][3]
  n_popO=options["n_samples"][4]


  #Define SLiM recipient pop index : 
  if recipient == '"A1"' :
    slim_pop_recipient = 5
    slim_pop_recipient_sister = 6
    n_r = options["n_samples"][0]
  elif recipient == '"A2"' :
    slim_pop_recipient = 6
    slim_pop_recipient_sister = 5
    n_r = options["n_samples"][1]
  elif recipient == '"B1"' :
    slim_pop_recipient = 7
    slim_pop_recipient_sister = 8
    n_r = options["n_samples"][2]
  elif recipient == '"B2"' :
    slim_pop_recipient = 8
    slim_pop_recipient_sister = 7
    n_r = options["n_samples"][3]
  elif recipient == '"O"' :
    slim_pop_recipient = 1
    n_r = options["n_samples"][4]
  else :
    print("wrong recipient population name")

  #Define SLiM donor pop index : 
  if donor == '"A1"' :
    slim_pop_donor = 5
    slim_pop_donor_sister = 6
  elif donor == '"A2"' :
    slim_pop_donor = 6
    slim_pop_donor_sister = 5
  elif donor == '"B1"' :
    slim_pop_donor = 7
    slim_pop_donor_sister = 8
  elif donor == '"B2"' :
    slim_pop_donor = 8
    slim_pop_donor_sister = 7
  elif donor == '"O"' :
    slim_pop_donor = 1
  else :
    print("wrong donor population name")

  #Define outgroup SLiM index
  outgroup_pop = 1

  # print program name
  introadapt.print_info(sys.argv[0],options["verbose"],sim=options["sim"])

  # set random seed:
  np.random.seed(options["seed_mut"])

  #IA mut file 
  mut_ia_list_file = pd.read_csv(job_dir+"/forwsim_freq_mut_don_in_rec_"+options["analysis"]+"_"+options["project"]+"_sim_"+options["sim"]+".txt")

  # read tree sequence from SLiM output file:
  treesq_all = tskit.load(job_dir+"/forwsim_final_"+options["analysis"]+"_"+options["project"]+"_sim_"+options["sim"]+".trees")
  
  #calcul mean ancestry % of donor pop in recipiant pop (adapt from MaLAdapt (Zhang et al.; 2023 : https://doi.org/10.1093/sysbio/syad033)
  mean_p1_all = introadapt.calc_p1ancestry(treesq_all, slim_pop_donor, slim_pop_recipient, n_r, g_between_end_and_mig)

  #Add neutral mutation :
  # Add neutral mutation in tree sequence
  add_mutation_start = time.time()
  msprime_seed = np.random.randint(1, 2**32-1)
  if options["verbose"]>=10 : print("SEED: " + str(msprime_seed) )
  mut_treesq = msprime.sim_mutations(treesq_all,
                                    rate = options["mu_neutral"],
                                    random_seed = msprime_seed,
                                    keep=True, 
                                    model=msprime.SLiMMutationModel(type=0))
  if options["verbose"]>=10 : print("Number of mutations " + str(mut_treesq.num_mutations))
  if options["verbose"]>=10 : print("Number of sites " + str(mut_treesq.num_sites))
  add_mutation_end = time.time()
  duration_add_mutation = add_mutation_end - add_mutation_start
  print("(mutsim) mutation duration : "+str(duration_add_mutation)+" s ")


  # Outgroup nodes samples
  sample_nodes_out = [x.id for x in mut_treesq.nodes() if ((x.time == 1) and (x.population == outgroup_pop))]
  # Donor sister nodes samples
  sample_nodes_non_don = [x.id for x in mut_treesq.nodes() if ((x.time == g_samp_don_s) and (x.population == slim_pop_donor_sister))]
  # Donor nodes samples
  sample_nodes_don = [x.id for x in mut_treesq.nodes() if ((x.time == g_samp_don) and (x.population == slim_pop_donor))]
  # non-recipient nodes samples (recipient sister, non introgressed)
  sample_nodes_non_rec = [x.id for x in mut_treesq.nodes() if ((x.time == 1) and (x.population == slim_pop_recipient_sister))]
  # recipient nodes samples (introgressed)
  sample_nodes_rec = [x.id for x in mut_treesq.nodes() if ((x.time == 1) and (x.population == slim_pop_recipient))]
  # All node samples
  sample_nodes = sample_nodes_out + sample_nodes_non_don + sample_nodes_don +sample_nodes_non_rec + sample_nodes_rec

  simplify_start = time.time()
  #simplify tree with sampled samples uniquely 
  mut_treesq_simp = mut_treesq.simplify(sample_nodes, keep_input_roots=True, filter_individuals = True, filter_populations = False)
  simplify_end = time.time()
  duration_simplify = simplify_end - simplify_start
  print("(mutsim) simplify duration : "+str(duration_simplify)+" s ")

  #Save tree sequences :
  mut_treesq_simp.dump(job_dir+"/Mut_TreeSeq_"+options["analysis"]+"_"+options["project"]+"_sim_"+options["sim"]+".trees")

  ####################
  ###Track ancestry###
  ####################

  ##Part of the script used to retrieve latent variables from the SLiM tree sequence for each recombination interval of the recipient population and to recover the proportion of introgression (in terms of nucleotides or recombination intervals).
  
  track_ancestry_start = time.time()
  #sampled recipient id node :
  node_recipient_sampled = [x.id for x in mut_treesq.nodes() if ((x.time == 1.0) and (x.population == slim_pop_recipient))]
  #Simplify the tree seq with uniquly the sampled recipient ind
  mut_treesq_rec_samp = mut_treesq.simplify(node_recipient_sampled, keep_input_roots=True, filter_individuals = True, filter_populations = False)
  #Number of mutation in the sampled recipient population :
  nbr_all_mut_in_samp_rec = mut_treesq_rec_samp.num_mutations

  sample_size = len(node_recipient_sampled)
  if prop_i_by_r_interv != "No" :
    #donor node id before the migration : 
    node_donor_mig = [x.id for x in mut_treesq.nodes() if ((x.time == g_mig) and (x.population == slim_pop_donor))]

    #Track recipient sampling nodes ancestry :
    edges_pop_donor_mig = mut_treesq.tables.link_ancestors(node_recipient_sampled, node_donor_mig)
    #Transforme ancestry table in dico :
    dict_edges_donor = edges_pop_donor_mig.asdict()
    #Remove metadata informations :
    list_del = ("metadata","metadata_offset","metadata_schema")
    list(map(dict_edges_donor.__delitem__, filter(dict_edges_donor.__contains__,list_del)))
    #Transforme dico in dataframe :
    data_edges_donor = pd.DataFrame({k:list(v) for k,v in dict_edges_donor.items()})
    #Sort by left column :
    data_edges_donor = data_edges_donor.sort_values(by=["left"])

    #################################################################################################################################
    ####### ia and i recombination interval and nucleotide lenght by ind 
    #################################################################################################################################
    node_recipient_mig = [x.id for x in mut_treesq.nodes() if ((x.time == g_mig) and (x.population == slim_pop_recipient))]
    #parent recipient pop edge (before mig)
    edges_pop_recipient_mig = mut_treesq.tables.link_ancestors(node_recipient_sampled, node_recipient_mig)
    #Transforme ancestry table in dico :
    dict_edges_recipient    = edges_pop_recipient_mig.asdict()
    #Remove metadata informations :
    list_del = ("metadata","metadata_offset","metadata_schema")
    list(map(dict_edges_recipient.__delitem__, filter(dict_edges_recipient.__contains__,list_del)))
    #Transforme dico in dataframe :
    data_edges_recipient     = pd.DataFrame({k:list(v) for k,v in dict_edges_recipient.items()})
    #Sort by left column :
    data_edges_recipient     = data_edges_recipient.sort_values(by=["left"])
    data_edge_all            = pd.concat([data_edges_recipient, data_edges_donor],axis=0)
    #control if introgression exist in tree seq
    len_data_edge_donor = len(data_edges_donor)
    if len_data_edge_donor != 0 :
        #Genome lenght 
        genome_size                              = options["L"]+1
        #recombination interval lenght (in nucleotids)
        data_edges_donor["recomb_int_n_lenght"] = data_edges_donor["right"]-data_edges_donor["left"]
        #take pair recipient node 
        recipient_pair_node = []
        for node in node_recipient_sampled :
            if node % 2 ==0 :
                recipient_pair_node.append(node)
        #check precense of ia mutations
        is_ai      = pd.isna(mut_ia_list_file["position"][0])
        if is_ai != True :
            #take their positions
            ia_mut_pos = mut_ia_list_file["position"]
            #create a empty dataframe for recomb int with ai mutation 
            concat_ia_recomb_int                     = pd.DataFrame(columns = ['left', 'right', 'parent', 'child'])
            #for each ia mutation get recom int 
            for pos in ia_mut_pos :
                ia_recomb_int         = data_edges_donor[((data_edges_donor['left']<=pos) & (data_edges_donor['right']>= pos))]
                concat_ia_recomb_int  = pd.concat([concat_ia_recomb_int, ia_recomb_int], axis=0)
            #remove duplicate (because some ia mut can be on the same recomb int)
            concat_ia_recomb_int_uniq = concat_ia_recomb_int.drop_duplicates()
            #biallelic genome size (in nucleotid)
            ind_genome_lenght_n      = genome_size*2
            list_node1               = []
            list_node2               = []
            list_n_prop_ia           = []
            list_n_prop_i            = []
            list_i_prop_recomb_int   = []
            liste_prop_recomb_int    = []
            for node in recipient_pair_node :
                #Get recombination interval by ind
                recomb_int_ind                   = data_edge_all.loc[((data_edge_all['child']==node) | (data_edge_all['child'] == node+1))]
                #get row size of the recombination interval
                recomb_int_len                   = len(recomb_int_ind)
                #Get IA recombination interval by ind
                concat_ia_recomb_int_uniq_by_ind = concat_ia_recomb_int_uniq.loc[((concat_ia_recomb_int_uniq['child']==node) | (concat_ia_recomb_int_uniq['child'] == node+1))]
                #i
                concat_i_recomb_int_uniq_by_ind  = data_edges_donor.loc[((data_edges_donor['child']==node) | (data_edges_donor['child'] == node+1))]
                #get IA row size of the recombination interval
                recom_int_ia_len                 = len(concat_ia_recomb_int_uniq_by_ind)
                #i
                recom_int_i_len                  = len(concat_i_recomb_int_uniq_by_ind)
                ##append all list
                #Get IA proportion (in recombination interval : nbr_of_ia_recomb_int/total_nbr_of_recomb_int)
                liste_prop_recomb_int.append(recom_int_ia_len/recomb_int_len)
                list_i_prop_recomb_int.append(recom_int_i_len/recomb_int_len)
                list_node1.append(node)
                list_node2.append(node+1)
                #Get ia proportion (in nucleotide : sum(ia_n_lenght)/total_n_genome_lenght)
                list_n_prop_ia.append(sum(concat_ia_recomb_int_uniq_by_ind["recomb_int_n_lenght"])/ind_genome_lenght_n)
                list_n_prop_i.append(sum(concat_i_recomb_int_uniq_by_ind["recomb_int_n_lenght"])/ind_genome_lenght_n)
            #save result in a dataframe
            prop_n_ia_recomb_int_by_ind = pd.DataFrame(list(zip(list_node1, list_node2, list_n_prop_ia, list_n_prop_i, liste_prop_recomb_int, list_i_prop_recomb_int )), columns =['node1', 'node2', 'n_prop_ia','n_prop_i' ,'prop_recomb_int', 'prop_recomb_int_i'])
            # Proportion of recombination interval with ia (mean, variance, max and min)
            mean_prop_ia_recom_in       = statistics.mean(prop_n_ia_recomb_int_by_ind["prop_recomb_int"])
            var_prop_ia_recom_in        = statistics.variance(prop_n_ia_recomb_int_by_ind["prop_recomb_int"])
            max_prop_ia_recom_in        = max(prop_n_ia_recomb_int_by_ind["prop_recomb_int"])
            min_prop_ia_recom_in        = min(prop_n_ia_recomb_int_by_ind["prop_recomb_int"])
            # Proportion of nucleotide with ia (mean, variance, max and min)
            mean_prop_ia_n              = statistics.mean(prop_n_ia_recomb_int_by_ind["n_prop_ia"])
            var_prop_ia_n               = statistics.variance(prop_n_ia_recomb_int_by_ind["n_prop_ia"])
            max_prop_ia_n               = max(prop_n_ia_recomb_int_by_ind["n_prop_ia"])
            min_prop_ia_n               = min(prop_n_ia_recomb_int_by_ind["n_prop_ia"])
            # Proportion of recombination interval with i (mean, variance, max and min)
            mean_prop_i_recom_in        = statistics.mean(prop_n_ia_recomb_int_by_ind["prop_recomb_int_i"])
            var_prop_i_recom_in         = statistics.variance(prop_n_ia_recomb_int_by_ind["prop_recomb_int_i"])
            max_prop_i_recom_in         = max(prop_n_ia_recomb_int_by_ind["prop_recomb_int_i"])
            min_prop_i_recom_in         = min(prop_n_ia_recomb_int_by_ind["prop_recomb_int_i"])
            # Proportion of nucleotide with i (mean, variance, max and min)
            mean_prop_i_n               = statistics.mean(prop_n_ia_recomb_int_by_ind["n_prop_i"])
            var_prop_i_n                = statistics.variance(prop_n_ia_recomb_int_by_ind["n_prop_i"])
            max_prop_i_n                = max(prop_n_ia_recomb_int_by_ind["n_prop_i"])
            min_prop_i_n                = min(prop_n_ia_recomb_int_by_ind["n_prop_i"])
        #if ia mut file is empty (no IA)
        else :
            #biallelic genome size (in nucleotid)
            ind_genome_lenght_n      = genome_size*2
            list_node1               = []
            list_node2               = []
            list_n_prop_i            = []
            list_i_prop_recomb_int   = []
            for node in recipient_pair_node :
                #Get recombination interval by ind
                recomb_int_ind                   = data_edge_all.loc[((data_edge_all['child']==node) | (data_edge_all['child'] == node+1))]
                #get nbr of recombination interval by ind
                recomb_int_len                   = len(recomb_int_ind)
                #get introgressed recombination interval by ind
                concat_i_recomb_int_uniq_by_ind = data_edges_donor.loc[((data_edges_donor['child']==node) | (data_edges_donor['child'] == node+1))]
                #get nbr of introgressed recombination interval by ind
                recom_int_i_len                 = len(concat_i_recomb_int_uniq_by_ind)
                ##append all list
                list_i_prop_recomb_int.append(recom_int_i_len/recomb_int_len)
                list_node1.append(node)
                list_node2.append(node+1)
                list_n_prop_i.append(sum(concat_i_recomb_int_uniq_by_ind["recomb_int_n_lenght"])/ind_genome_lenght_n)
            #save result in a dataframe
            prop_i_recomb_int_by_ind = pd.DataFrame(list(zip(list_node1, list_node2, list_n_prop_i, list_i_prop_recomb_int )), columns =['node1', 'node2','n_prop_i', 'prop_recomb_int_i'])
            mean_prop_i_recom_in     = statistics.mean(prop_i_recomb_int_by_ind["prop_recomb_int_i"])
            var_prop_i_recom_in      = statistics.variance(prop_i_recomb_int_by_ind["prop_recomb_int_i"])
            max_prop_i_recom_in      = max(prop_i_recomb_int_by_ind["prop_recomb_int_i"])
            min_prop_i_recom_in      = min(prop_i_recomb_int_by_ind["prop_recomb_int_i"])
            mean_prop_i_n            = statistics.mean(prop_i_recomb_int_by_ind["n_prop_i"])
            var_prop_i_n             = statistics.variance(prop_i_recomb_int_by_ind["n_prop_i"])
            max_prop_i_n             = max(prop_i_recomb_int_by_ind["n_prop_i"])
            min_prop_i_n             = min(prop_i_recomb_int_by_ind["n_prop_i"])
            mean_prop_ia_recom_in    = 0.0
            var_prop_ia_recom_in     = 0.0
            max_prop_ia_recom_in     = 0.0
            min_prop_ia_recom_in     = 0.0
            mean_prop_ia_n           = 0.0
            var_prop_ia_n            = 0.0
            max_prop_ia_n            = 0.0
            min_prop_ia_n            = 0.0
    #if data_edges_donor is empty (no I and no IA)
    else :
        mean_prop_ia_recom_in = 0.0
        var_prop_ia_recom_in  = 0.0
        max_prop_ia_recom_in  = 0.0
        min_prop_ia_recom_in  = 0.0
        mean_prop_ia_n        = 0.0
        var_prop_ia_n         = 0.0
        max_prop_ia_n         = 0.0
        min_prop_ia_n         = 0.0
        mean_prop_i_recom_in  = 0.0
        var_prop_i_recom_in   = 0.0
        max_prop_i_recom_in   = 0.0
        min_prop_i_recom_in   = 0.0
        mean_prop_i_n         = 0.0
        var_prop_i_n          = 0.0
        max_prop_i_n          = 0.0
        min_prop_i_n          = 0.0

    #################################################################################################################################
    #######Unique introgressed recomb interval for the pop 
    #################################################################################################################################

    ###Compt donor nodes by unique interval :
    #Array with unique recombination point : 
    array_unique_int_do   = np.sort(pd.concat([data_edges_donor['left'],data_edges_donor['right']]).unique())
    #create unique intervalle dataframe :
    data_interv_unique_do = pd.DataFrame()
    data_interv_unique_do["left"]=array_unique_int_do[:-1]
    data_interv_unique_do["right"]=array_unique_int_do[1:]

    ###Obtain list with for all unique intervals parent donor nodes comptage :####VERY LONG####
    #Empty list for sum of node with introgression for each recomb intervalle :
    list_node_compt_by_uniq_interv =[]
    #Arrays with each recomb intervalles 
    interv_do                         = pd.IntervalIndex.from_arrays(data_edges_donor.left, data_edges_donor.right, closed='both')
    interv_do_uniq                    = pd.IntervalIndex.from_arrays(data_interv_unique_do.left, data_interv_unique_do.right, closed='neither')
    #Sum for each unique recomb intervalle of donor intervalle that overlap it :
    list_node_compt_by_uniq_interv    =[np.sum(interv_do.overlaps(interv_do_uniq[x])==True) for x in range(len(interv_do_uniq))]
    #Create a column with number of introgressed node by unique intervalles.
    data_interv_unique_do["compt_do"] = list_node_compt_by_uniq_interv

    #merge consecutive duplicate :
    cols             = ['compt_do']
    merge_cons_dupli = data_interv_unique_do[(data_interv_unique_do[cols].shift() != data_interv_unique_do[cols]).any(axis=1)].reset_index()
    list_real_right  = []
    if len(merge_cons_dupli)> 0 :
      for ri in range(len(merge_cons_dupli)-1) : 
            list_real_right.append(merge_cons_dupli.left[ri+1])
      list_real_right.append(data_interv_unique_do.right[len(data_interv_unique_do)-1])
      merge_cons_dupli["right"]            = list_real_right
      merge_cons_dupli["proportion"]       = merge_cons_dupli["compt_do"]/sample_size
      merge_cons_dupli                     = merge_cons_dupli.drop(columns=['index'])
      #remove interval with 0 intro 
      merge_cons_dupli                     = merge_cons_dupli[merge_cons_dupli["compt_do"]!=0]
      #span intro calcul 
      merge_cons_dupli["intro_span"]       = merge_cons_dupli["right"]-merge_cons_dupli["left"]
      #proportion intro calcul 
      merge_cons_dupli["proportion_intro"] = merge_cons_dupli["proportion"]*merge_cons_dupli["intro_span"]
      merge_cons_dupli["len_intro"]        = merge_cons_dupli["compt_do"]*merge_cons_dupli["intro_span"]
    else : 
      merge_cons_dupli                     = pd.DataFrame({'left': [0.0] , 'right': [geno_len], 'compt_do' : [0], 'proportion' : [0], 'intro_span' : [0],  'proportion_intro' : [0], 'len_intro' : [0]})
    merge_cons_dupli.to_csv(job_dir+"/Recomb_interval_intro_"+options["analysis"]+"_"+options["project"]+"_sim_"+options["sim"]+".csv")
  else :
    merge_cons_dupli = pd.DataFrame({'left': [0.0] , 'right': [geno_len], 'compt_do' : [0], 'proportion' : [0], 'intro_span' : [0],  'proportion_intro' : [0], 'len_intro' : [0]})

  ###Add sampled recipient ind mutation number in csv :
  csv_mut_rec_samp                          = pd.read_csv(job_dir+"/forwsim_latent_variable_"+options["analysis"]+"_"+options["project"]+"_sim_"+options["sim"]+".txt")
  #All mutation (sweep + neutral + introgressed mutation) in sampled recipient population :
  csv_mut_rec_samp['all_mut_in_rec_s']      = nbr_all_mut_in_samp_rec
  #Ratio between advantageous introgressed mutations and all mutation in sampled recipient population :
  csv_mut_rec_samp['prop_mut_do_in_rec']    = csv_mut_rec_samp['uniq_mut_do_in_rec_s']/nbr_all_mut_in_samp_rec
  csv_mut_rec_samp['mean_p1']               = mean_p1_all
  csv_mut_rec_samp['mean_prop_ia_recom_in'] = mean_prop_ia_recom_in
  csv_mut_rec_samp['var_prop_ia_recom_in']  = var_prop_ia_recom_in
  csv_mut_rec_samp['max_prop_ia_recom_in']  = max_prop_ia_recom_in
  csv_mut_rec_samp['min_prop_ia_recom_in']  = min_prop_ia_recom_in
  csv_mut_rec_samp['mean_prop_i_recom_in']  = mean_prop_i_recom_in
  csv_mut_rec_samp['var_prop_i_recom_in']   = var_prop_i_recom_in
  csv_mut_rec_samp['max_prop_i_recom_in']   = max_prop_i_recom_in
  csv_mut_rec_samp['min_prop_i_recom_in']   = min_prop_i_recom_in
  csv_mut_rec_samp['mean_prop_ia_n']        = mean_prop_ia_n
  csv_mut_rec_samp['var_prop_ia_n']         = var_prop_ia_n
  csv_mut_rec_samp['max_prop_ia_n']         = max_prop_ia_n
  csv_mut_rec_samp['min_prop_ia_n']         = min_prop_ia_n
  csv_mut_rec_samp['mean_prop_i_n']         = mean_prop_i_n
  csv_mut_rec_samp['var_prop_i_n']          = var_prop_i_n
  csv_mut_rec_samp['max_prop_i_n']          = max_prop_i_n
  csv_mut_rec_samp['min_prop_i_n']          = min_prop_i_n
  if len(merge_cons_dupli)> 0 :
    csv_mut_rec_samp["prop_intro_all"]      = sum(merge_cons_dupli["proportion_intro"])/geno_len
    #calcul the span intro for all the sample :
    csv_mut_rec_samp["len_intro_all"]       = sum(merge_cons_dupli["len_intro"])
  else :
    csv_mut_rec_samp["prop_intro_all"]      = [0.0]
    csv_mut_rec_samp["len_intro_all"]       = [0]
  csv_mut_rec_samp.to_csv(job_dir+"/forwsim_latent_variable_"+options["analysis"]+"_"+options["project"]+"_sim_"+options["sim"]+".txt", index=False)
  
  track_ancestry_end = time.time()
  duration_track_ancesty = track_ancestry_end - track_ancestry_start
  print("(mutsim) track ancestry duration : "+str(duration_track_ancesty)+" s ")

  mutsim_end = time.time()
  duration_mutsim = mutsim_end - mutsim_start
  print("(mutsim) mutsim duration : "+str(duration_mutsim)+" s ")
############################################################################################################
if __name__ == "__main__":
    main()


