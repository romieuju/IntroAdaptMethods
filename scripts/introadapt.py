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


###Script containing the functions used in the Python scripts in the pipeline 

import configparser
import pandas as pd
import numpy as np
import scipy.stats as st
import allel
import pytest

### PRINT INFO ######################################################################################
def print_info(script_name,verbose,project=None,batch=None,sim=None):
  if verbose >=1 :
    print("#########################################")
  if verbose >=0 :
    if batch is not None:
      project = None
      if sim is None :
        print("TimeAdapt - %s - batch %s" %(script_name,batch) )
      else : 
        print("TimeAdapt - %s - batch %s - simulation %s" %(script_name,batch,sim) )
    elif project is not None:
      print("TimeAdapt - %s - project %s" %(script_name,project) )
    else :
      print("TimeAdapt - %s" %script_name )
  if verbose >=1 :
    print("by Jules Romieu")
    print("INRAE & Montpellier universitet")
    print("jules.romieu@umontpellier.fr")
    print("#########################################")

### GET PROJECT OPTIONS ###############################################################################
def get_project_options(proj_options_file):
  proj_options = configparser.ConfigParser()
  proj_options.read(proj_options_file)

  # Settings
  assert 'Settings' in proj_options,"Missing [Settings] section in file"
  assert 'project' in proj_options['Settings'],"Missing 'project' parameter in file"
  simul_type         = proj_options.get('Settings','simul_type')
  project            = proj_options.get('Settings','project')
  config_file        = proj_options.get('Settings','config_file')
  analysis           = proj_options.get('Settings','analysis')
  results_dir        = proj_options.get('Settings','results_dir')
  try:
    verbose      = proj_options.getint('Settings','verbose')
  except:
    verbose      = int(proj_options.getfloat('Settings','verbose'))

  # Model
  assert 'Model' in proj_options,"Missing [Model] section in project options file"
  recipient             = proj_options.get('Model', 'recipient')
  donor                 = proj_options.get('Model', 'donor')
  slim_split            = proj_options.getint('Model', 'slim_split')

 # Sample
  assert 'Sample' in proj_options,"Missing [Sample] section in project options file"
  sample_size           = proj_options.getint('Sample','size')
  n_samples             = [int(i) for i in proj_options.get('Sample',"n_samples").split()]
  
  # Genome
  assert 'Genome' in proj_options,"Missing [Genome] section in project options file"
  nchr                  = proj_options.getint('Genome','nchr')
  chr_ends              = [int(i) for i in proj_options.get("Genome","chr_ends").split()]
  L                     = proj_options.getint("Genome",'L')

  # Inference methods
  assert 'Inference' in proj_options,"Missing [Inference] section in project options file"
  volcano_switch        = proj_options.get('Inference', 'volcano_switch')
  maladapt_switch       = proj_options.get('Inference', 'maladapt_switch')
  maladapt_trained_path = proj_options.get('Inference', 'maladapt_trained_path')
  maladapt_features     = [i for i in proj_options.get('Inference', 'maladapt_features').split()]
  genomatnn_switch      = proj_options.get('Inference', 'genomatnn_switch')
  genomatnn_trained_mod = [i for i in proj_options.get('Inference', 'genomatnn_trained_mod').split()]
  wind_reftable_type    = proj_options.get('Inference', 'wind_reftable_type')
  volcano_threshold     = proj_options.get('Inference', 'volcano_threshold')
  genomatnn_threshold   = proj_options.get('Inference', 'genomatnn_threshold')
  maladapt_threshold    = proj_options.get('Inference', 'maladapt_threshold')
  q95_threshold         = proj_options.get('Inference', 'q95_threshold')
  
  #  Statistics 
  assert 'Statistics' in proj_options,"Missing [Statistics] section in project options file"
  zns                      = proj_options.get('Statistics', 'zns')
  donor_pop_stat           = proj_options.get('Statistics', 'donor_pop_stat')
  prop_i_by_r_interv       = proj_options.get('Statistics', 'prop_i_by_r_interv')
  Is_phased                = proj_options.get('Statistics','isphased')
  window_size              = proj_options.getint('Statistics', 'window_size')
  window_start             = proj_options.getint('Statistics', 'window_start')
  window_end               = proj_options.getint('Statistics', 'window_end')
  window_step              = proj_options.getint('Statistics', 'window_step')

  # TODO
  return {"simul_type" : simul_type,
          "project" : project,
          "analysis" : analysis,
          "config_file" : config_file,
          "results_dir" : results_dir,
          "verbose" : verbose,
          "recipient": recipient,
          "donor": donor,
          "slim_split": slim_split,
          "nchr" : nchr,
          "chr_ends" : chr_ends,
          "L" : L,
          "sample_size": sample_size,
          "n_samples" : n_samples,
          "window_size" : window_size,
          "window_start" : window_start,
          "window_end" : window_end,
          "window_step" : window_step, 
          "volcano_switch" : volcano_switch,
          "maladapt_switch": maladapt_switch,
          "maladapt_trained_path": maladapt_trained_path,
          "maladapt_features": maladapt_features,
          "genomatnn_switch":genomatnn_switch,
          "genomatnn_trained_mod":genomatnn_trained_mod,
          "wind_reftable_type":wind_reftable_type,
          "volcano_threshold":volcano_threshold,
          "genomatnn_threshold":genomatnn_threshold,
          "maladapt_threshold":maladapt_threshold,
          "q95_threshold":q95_threshold,
          "zns": zns,
          "donor_pop_stat":donor_pop_stat,
          "prop_i_by_r_interv":prop_i_by_r_interv,
          "is_phased" : Is_phased}

### GET SIM OPTIONS ######################################################################################
def get_sim_options(sim_options_file):
  sim_options = configparser.ConfigParser()
  sim_options.read(sim_options_file)
  sim                         = sim_options.get('Simulation', 'sim')
  coalescence_split           = sim_options.get('Demography','coalescence_split')
  slim_split                  = sim_options.getint('Demography','slim_split')
  popAncO                     = sim_options.getint('Demography','popAncO')
  popO                        = sim_options.getint('Demography','popO')
  popAncAB                    = sim_options.getint('Demography','popAncAB')
  popA1                       = sim_options.getint('Demography','popA1')
  popB1                       = sim_options.getint('Demography','popB1')
  popAncA                     = sim_options.getint('Demography','popAncA')
  popA2                       = sim_options.getint('Demography','popA2')
  popAncB                     = sim_options.getint('Demography','popAncB')
  popB2                       = sim_options.getint('Demography','popB2')
  generation_split_OA         = sim_options.getint('Demography','generation_split_OA')
  generation_split_AB         = sim_options.getint('Demography','generation_split_AB')
  generation_split_A          = sim_options.getint('Demography','generation_split_A')
  generation_split_B          = sim_options.getint('Demography','generation_split_B')
  generations_forward         = sim_options.getint('Demography', 'generations_forward')
  sample_time_d               = sim_options.getint('Demography', 'sample_time_d')
  sample_time_d_s             = sim_options.getint('Demography', 'sample_time_d_s')
  generations_migration_start = sim_options.getint('Demography','generations_migration_start')
  generations_migration_end   = sim_options.getint('Demography','generations_migration_end')
  coalescence_split           = sim_options.get('Demography', "coalescence_split")
  mu_neutral                  = sim_options.getfloat('Genome','mu_neutral')
  assert mu_neutral>=0, "Verify neutral mutation rate, it can be only positive values and zero"
  mu_advantageous             = sim_options.getfloat('Genome','mu_advantageous')
  assert mu_advantageous>=0, "Verify advantageous mutation rate, it can be only positive values and zero"
  ttratio                     = sim_options.getfloat('Genome','ttratio')
  seq_error                   = sim_options.getfloat('Genome','seq_error')
  msprime_r_map_positions     = [int(i) for i in sim_options.get("Genome","msprime_r_map_positions").split()]
  msprime_r_map_rates         = [float(i) for i in sim_options.get("Genome","msprime_r_map_rates").split()]
  seed_coal                   = sim_options.getint('Seeds','seed_coal')
  seed_mut                    = sim_options.getint('Seeds','seed_mut')
  return {"sim":sim,
          "coalescence_split":coalescence_split,
          "slim_split":slim_split,
          "popAncO":popAncO,
          "popO":popO,
          "popAncAB":popAncAB,
          "popA1":popA1,
          "popB1":popB1,
          "popAncA":popAncA,
          "popA2":popA2,
          "popAncB":popAncB,
          "popB2":popB2,
          "generation_split_OA":generation_split_OA,
          "generation_split_AB":generation_split_AB,
          "generation_split_A":generation_split_A,
          "generation_split_B":generation_split_B,
          "coalescence_split":coalescence_split,
          "generations_forward":generations_forward,
          "sample_time_d":sample_time_d,
          "sample_time_d_s":sample_time_d_s,
          "generations_migration_start":generations_migration_start,
          "generations_migration_end":generations_migration_end,
          "mu_neutral":mu_neutral,
          "mu_advantageous":mu_advantageous,
          "ttratio":ttratio, 
          "seq_error":seq_error,
          "msprime_r_map" : {"positions" : msprime_r_map_positions,
                             "rates" : msprime_r_map_rates},
          "seed_coal":seed_coal, 
          "seed_mut":seed_mut}

### GET OPTIONS ######################################################################################
def get_options(proj_options_file,sim_options_file):
  options = get_project_options(proj_options_file)
  options.update(get_sim_options(sim_options_file))
  return options



### Define window for calculated summary statistics
def windows_funct(len_gen, wind, step) :
    start_w = 0
    windows = []
    start_pos = []
    i=0
    count=0
    while i <= len_gen:
        if (count%wind==1)& (start_w == 0):
            start_w +=1
            start_pos.append(i)
            i+=1
            count+=1
        if (count%wind==0) & (start_w > 0):
            end_pos = i
            if start_w > 0:
                windows.append([start_pos[0],end_pos])
            start_w = 0
            i=start_pos[0]+step
            start_pos = []
            count = 0
        else:
            i+=1
            count+=1
    return windows

#Methode use to compte parent node by uniq interval, need unique interval and interval by parent node (use in mutsim.py)
def compt_node_by_uniq_interv(uniq_interv_df, ancestry_df) :
    list_comp =[]
    #for each interval (x) in unique :
    for x in range(len(uniq_interv_df)) :
        #creat comptor for parent donnor node :
        compt = 0
        #Create an interval between left and right :
        x_int = pd.Interval(uniq_interv_df.left[x], uniq_interv_df.right[x])
        #for each interval in donnor ancestry table (i)
        for i in range(len(ancestry_df)) :
            #Create an interval between left and right :
            i_interv = pd.Interval(ancestry_df.left[i], ancestry_df.right[i])
            #If the unique interval overlaps the ancestry interval, compt += 1  :
            if i_interv.overlaps(x_int) == True :
                compt +=1
        #append on the unique interval index the comptage 
        list_comp.append(compt)
    return list_comp

#Proportion of introgression from MaLAdapt (meanp1 in original script)
def calc_p1ancestry (tree_seq, id_pop_d, id_pop_r, n_pop_r_size, t_sinceadm):
    ts = tree_seq
    any_ancestry = ancestry_p_varies(ts, id_pop_d, id_pop_r, n_pop_r_size, t_sinceadm)
    meanp1 = sum(any_ancestry)/len(any_ancestry)
    return meanp1

def ancestry_p_varies(ts, idpopd, idpopr, npoprsize, duration): #pop=source pop
    n=npoprsize*2
    mixtime=duration
    p = [x.id for x in ts.nodes() if ((x.population == int(idpopd)) and (x.time == mixtime))] #source pop
    today = [x.id for x in ts.nodes() if ((x.population == idpopr) and (x.time == 1))] #assuming p3 is recipient
    tree_p = [sum([t.num_tracked_samples(u) for u in p])/n for t in ts.trees(tracked_samples=today, sample_counts=True)]
    return tree_p


























### MAKE EMPTY GENOTYPE ARRAY ····························································
def empty_genotype_array(n_loci, n_samples, ploidy=2, allele=-1):
  """
  Creates a genotype array with all values as missing (-1) for a given number
  of samples, loci and ploidy

  :return: empty_ga
  """
  empty_ga = allel.GenotypeArray(np.full((n_loci, n_samples, ploidy), allele, dtype='i1'), dtype='i1')
  return empty_ga
def test_empty_genotype_array():
  test_ga = empty_genotype_array(3, 4, ploidy=2, allele=-1)
  assert type(test_ga) is allel.model.ndarray.GenotypeArray
  assert len(test_ga) == 3
  assert all(test_ga[0,0] == [-1,-1])
  assert all(test_ga[2,3] == [-1,-1])
  test_ga = empty_genotype_array(3, 4, ploidy=2, allele=0)
  assert all(test_ga[1,1] == [0,0])
### end MAKE EMPTY GENOTYPE ARRAY ····························································

### SNP CALLING FROM SIMULATED READS (WITH SEQUENCING ERROR)  ····························
def snp_calling(true_genotype, f_num_reads, error_rate=0.005, reads_th=8, score_th=5, ratio_th=10, damage=False, transversion=True):
    """
    snp_calling function takes perfect simulated data from one locus of one 
    diploid individual and adds missing data and error according to the number 
    of reads of the site, error rate of the sequencing technology and, for 
    ancient DNA not sequenced from damage repair (dr) libraries, creates 
    missing data for transition SNPs (since they cannot be distinguished from
    aDNA damage)
    :param true_genotype:
    :param f_num_reads:
    :param error_rate:
    :param reads_th:
    :param score_th:
    :param ratio_th:
    :param dr:
    :param transversion:
    :return:
    """
    if damage is True and transversion is False:
        genotype_call = [-1, -1]
    elif f_num_reads >= reads_th:
        derived_count = sum(true_genotype)
        p_derived = float(derived_count) / 2. * (1 - float(error_rate)) + (1 - float(derived_count) / 2.) * float(error_rate)
        derived_reads = st.binom.rvs(f_num_reads, p_derived)
        ancestral_reads = f_num_reads - derived_reads
        if f_num_reads >= (score_th*2):
            if derived_reads == 0:
                genotype_call = [0, 0]
            elif ancestral_reads == 0:
                genotype_call = [1, 1]
            else:
                if (derived_reads >= score_th) & (ancestral_reads < score_th):
                  genotype_call = [1, 1]
                elif (derived_reads < score_th) & (ancestral_reads >= score_th):
                  genotype_call = [0, 0]
                elif (derived_reads >= score_th) & (ancestral_reads >= score_th):
                  ratio_of_scores = derived_reads / ancestral_reads
                  if (ratio_of_scores >= 1 / ratio_th) & (ratio_of_scores <= ratio_th):
                    if (derived_count == 1):
                      genotype_call = true_genotype
                    elif (st.binom.rvs(1, 0.5) == 1):
                      genotype_call = [0, 1]
                    else:
                      genotype_call = [1, 0]
                  elif derived_reads > ancestral_reads:
                    genotype_call = [1, 1]
                  else:
                    genotype_call = [0, 0]
    else:
        genotype_call = [-1, -1]
    return genotype_call
def test_snp_calling():
  np.random.seed(1234)
  genotype_call = snp_calling( [0, 1], 100, error_rate=0.005, reads_th=1,
                score_th=10, ratio_th=3, damage=False, transversion=True)
  assert genotype_call == [0,1]
  genotype_call = snp_calling( [0, 1], 1, error_rate=0.005, reads_th=10,
                score_th=10, ratio_th=3, damage=False, transversion=True)
  assert genotype_call == [-1,-1]
### end SNP CALLING FROM SIMULATED READS (WITH SEQUENCING ERROR)  ····························

### SIMULATE SEQUENCING  ····························
def sequencing(ts, ssize, ttr, seq_error, damage, cov):
  if len(cov) != ssize:
    msg = "Number of coverage values (length=" + str(len(cov)) + \
          ") and number of samples (ssize=" + str(ssize) + \
          ") do not match"
    raise ValueError(msg)

  geno_data = empty_genotype_array(n_loci=ts.num_sites,
                                   n_samples=ssize,
                                   ploidy=2)
  positions = []
  locus = 0
  for variant in ts.variants():
    positions.append(round(variant.position))
    #print(variant.position)
    var_genotypes = variant.genotypes
    # print("--------------------------------")
    # print(var_genotypes)
    num_reads = np.random.poisson(lam=cov, size=ssize)
    transversion_snp = True
    if np.random.random() < ttr / (ttr + 1):
      transversion_snp = False
    # print(num_reads)
    for i in range(0, 2 * ssize, 2):
      if len(variant.alleles)==2:
        genotype_call = snp_calling(true_genotype=var_genotypes[i:(i + 2)],
                                                  f_num_reads=num_reads[int(i / 2)],
                                                  error_rate=seq_error,
                                                  damage=damage[int(i / 2)],
                                                  transversion=transversion_snp)
      else:
        genotype_call = [-1, -1] # this removes all SNP with more than two alleles
      geno_data[locus, int(i / 2)] = genotype_call
    locus = locus + 1
    #print(locus)
  return geno_data, positions
#def test_sequencing():
  # np.random.seed(1234)
  # TODO create a ts and some test from it
### end SIMULATE SEQUENCING  ····························

### TEST DATA FOR SUMMARY STATISTISCS
#                                  ind 0   ind 1   ind 2   ind 3  
test_ga_A = allel.GenotypeArray([[[ 0, 0],[ 0, 0],[ 0, 0],[ 0, 0]], # locus 0
                                 [[ 0, 1],[ 0, 1],[ 0, 1],[ 1, 1]], # locus 1
                                 [[-1,-1],[-1,-1],[-1,-1],[-1,-1]], # locus 2
                                 [[ 1, 1],[ 1, 1],[ 0, 1],[ 1, 1]], # locus 3
                                 [[ 0, 1],[ 0, 0],[ 0, 1],[ 0,-1]]],# locus 4
                                 dtype='i1')
#              0   1   2   3   4
test_pos_A = (10,123,234,299,340)
#                                  ind 0   ind 1   ind 2   ind 3  
test_ga_B = allel.GenotypeArray([[[ 0, 1],[ 0, 1],[ 1, 1],[ 0, 0]], # locus 0
                                 [[ 0, 1],[ 0, 1],[ 0, 1],[ 1, 1]], # locus 1
                                 [[ 0, 1],[ 0, 1],[-1,-1],[ 1, 1]], # locus 2
                                 [[ 0, 1],[ 0, 0],[ 0, 1],[ 1, 1]], # locus 3
                                 [[ 0, 1],[ 0, 1],[ 0, 0],[ 1, 1]], # locus 4
                                 [[ 0, 1],[ 0, 1],[ 0, 0],[ 1, 1]], # locus 5
                                 [[ 1, 1],[ 0, 1],[ 0, 1],[ 1, 1]], # locus 6
                                 [[ 0, 1],[ 0, 1],[ 0, 0],[ 1, 1]], # locus 7
                                 [[ 1, 1],[ 0, 1],[ 0, 0],[ 0, 1]], # locus 8
                                 [[ 0, 1],[ 0, 1],[ 0, 0],[ 1, 1]], # locus 9
                                 [[ 1, 1],[ 0, 1],[ 0, 1],[ 1, 1]], # locus 10
                                 [[ 0, 1],[ 0, 1],[ 0, 0],[ 1, 1]], # locus 11
                                 [[ 0, 0],[ 0, 1],[ 0, 1],[ 0, 1]], # locus 12
                                 [[-1,-1],[-1,-1],[-1,-1],[ 1, 1]], # locus 13
                                 [[ 1, 1],[ 1, 1],[ 0, 0],[ 1, 1]], # locus 14
                                 [[ 1, 1],[ 1, 1],[ 0, 0],[ 1, 1]], # locus 15
                                 [[ 1, 1],[ 0, 1],[ 0, 0],[ 0, 1]], # locus 16
                                 [[ 1, 1],[ 1, 1],[ 1, 0],[ 1, 1]], # locus 17
                                 [[ 1, 1],[ 1, 1],[ 0, 0],[ 1, 1]], # locus 18
                                 [[ 0, 1],[ 0, 0],[ 0, 1],[-1,-1]]],# locus 19
                                 dtype='i1')
#             0  1  2  3  4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19
test_pos_B = (4,10,50,77,99,123,134,150,178,201,209,234,256,270,299,311,315,340,358,378)

#                                 ga,        pos, nchr, chr_ends, w_size, expected_S
testdata_S = [pytest.param(test_ga_A, test_pos_A,    1,    [400],     50,          3, id="A"),
              pytest.param(test_ga_B, test_pos_B,    1,    [400],     50,         19, id="B")]
#                                 ga,        pos, nchr, chr_ends, w_size, expected_Pi
testdata_Pi = [pytest.param(test_ga_A, test_pos_A,    1,    [400],     50, 0.00315476, id="A"),
               pytest.param(test_ga_B, test_pos_B,    1,    [400],     50, 0.02418452, id="B")]
#                                 ga,        pos, nchr, chr_ends, w_size,  expected_WT
testdata_WT = [pytest.param(test_ga_A, test_pos_A,    1,    [400],     50, 0.00289256, id="A"),
               pytest.param(test_ga_B, test_pos_B,    1,    [400],     50, 0.01831956, id="B")]







### SAVE RESULTS from st.describe() INTO DICT ····························
def save_moments_2_dict(moments,sumstats,sample_name,sep,stat_name):
  sumstats[sample_name+sep+"min"+stat_name] = float(np.ma.getdata(moments[1][0]))
  sumstats[sample_name+sep+"max"+stat_name] = float(np.ma.getdata(moments[1][1]))
  sumstats[sample_name+sep+"m"+stat_name] = float(np.ma.getdata(moments[2]))
  sumstats[sample_name+sep+"v"+stat_name] = float(np.ma.getdata(moments[3]))
  sumstats[sample_name+sep+"s"+stat_name] = float(np.ma.getdata(moments[4]))
  sumstats[sample_name+sep+"k"+stat_name] = float(np.ma.getdata(moments[5]))
### end SAVE RESULTS from st.describe() INTO DICT  ····························

### CALCULATE SUMMARY STATISTICS : SINGLE SAMPLE  ····························
def single_sample_sumstats(ga,pos,nchr,chr_ends,w_size,sumstats,name="",sep="_",quiet=True):
  ac = ga.count_alleles()
  sfs = allel.sfs_folded(ac)
  # total number of segregating sites (excluding monomorphic and all missing data)
  segsites = sum(sfs)-sfs[0]
  sumstats[name+sep+"S"]=segsites
  if quiet is False: print("Segregating sites: "+ str(segsites))
  # Observed heterozygosity and Inbreeding coefficient
  ho = st.describe(allel.heterozygosity_observed(ga), nan_policy='omit')
  save_moments_2_dict(ho,sumstats,name,sep,"Ho")
  fis = st.describe(allel.inbreeding_coefficient(ga), nan_policy='omit')
  save_moments_2_dict(fis,sumstats,name,sep,"Fis")
  if quiet is False: print("Ho: " + str(ho) + "; Fis: " + str(fis) )
  # pairwise genetic diversity
  total_pi = allel.sequence_diversity(pos, ac, start=1, stop=chr_ends[nchr-1])
  sumstats[name+sep+"Pi"]=total_pi
  w_pi, _, _, _ = allel.windowed_diversity(pos, ac, size=w_size, start= 1, stop=chr_ends[0])
  for chromo in range(1,nchr):
    temp, _, _, _ = allel.windowed_diversity(pos, ac, size=w_size, start= 1+chr_ends[chromo-1], stop=chr_ends[chromo])
    np.append(w_pi,temp)
  pi = st.describe(w_pi)
  save_moments_2_dict(pi,sumstats,name,sep,"Pi")
  if quiet is False: print("π: "+ str(total_pi) + " ; " + str(pi))
  # Watterson theta (from number of segregating sites)
  w_W_theta, _, _, _ = allel.windowed_watterson_theta(pos, ac, size=w_size, start=1, stop=chr_ends[0])
  for chromo in range(1,nchr):
    temp, _, _, _ = allel.windowed_watterson_theta(pos, ac, size=w_size, start=1+chr_ends[chromo-1], stop=chr_ends[chromo])
    np.append(w_W_theta,temp)
  W_theta = st.describe(w_W_theta)
  save_moments_2_dict(W_theta,sumstats,name,sep,"WT")
  if quiet is False: print("Watterson θ: "+ str(W_theta))
  # Tajima's D
  total_Taj_D = allel.tajima_d(ac, pos, start=1, stop=chr_ends[nchr-1])
  sumstats[name+sep+"TD"]=total_Taj_D
  w_Taj_D, _, _ = allel.windowed_tajima_d(pos, ac, size=w_size, start=1, stop=chr_ends[0])
  for chromo in range(1,nchr):
    temp, _, _ = allel.windowed_tajima_d(pos, ac, size=w_size, start=1+chr_ends[chromo-1], stop=chr_ends[chromo])
    np.append(w_Taj_D,temp)
  Taj_D = st.describe(w_Taj_D, nan_policy='omit')
  save_moments_2_dict(Taj_D,sumstats,name,sep,"TD")
  if quiet is False: print("Tajima's D: "+ str(total_Taj_D)+ " ; " + str(Taj_D))
  # distribution of sizes for naive runs of homozygosity (distance between heterozygous positions)
  roh_distribution = np.full(int(round(np.log10(w_size)))+1, 0)
  roh_distribution = roh_distribution + windowed_distribution_roh(ga, pos, w_size, start=1, stop=chr_ends[0])
  for chromo in range(1,nchr):
    roh_distribution + windowed_distribution_roh(ga, pos, w_size, start=1+chr_ends[chromo-1], stop=chr_ends[chromo])
  for i in range(0,len(roh_distribution)):
    sumstats[name+sep+"RoHD"+sep+str(i)]=roh_distribution[i]
  if quiet is False: print("Distribution of Runs of Homozygosity: "+ str(roh_distribution))
  return

@pytest.mark.parametrize("ga,pos,nchr,chr_ends,w_size,expected_S", testdata_S)
def test_single_sample_sumstats_S(ga,pos,nchr,chr_ends,w_size,expected_S):
  test_sumstats = {}
  single_sample_sumstats(ga, pos, nchr, chr_ends, w_size, test_sumstats)
  assert test_sumstats["_S"] == expected_S
  assert test_sumstats["_S"] >= 0

@pytest.mark.parametrize("ga,pos,nchr,chr_ends,w_size,expected_Pi", testdata_Pi)
def test_single_sample_sumstats_Pi(ga,pos,nchr,chr_ends,w_size,expected_Pi):
  test_sumstats = {}
  single_sample_sumstats(ga, pos, nchr, chr_ends, w_size, test_sumstats)
  assert test_sumstats["_Pi"] == pytest.approx(expected_Pi)
  assert test_sumstats["_Pi"] >= 0
  assert test_sumstats["_mPi"] == pytest.approx(expected_Pi)

@pytest.mark.parametrize("ga,pos,nchr,chr_ends,w_size,expected_WT", testdata_WT)
def test_single_sample_sumstats_WT(ga,pos,nchr,chr_ends,w_size,expected_WT):
  test_sumstats = {}
  single_sample_sumstats(ga, pos, nchr, chr_ends, w_size, test_sumstats)
  assert test_sumstats["_mWT"] >= 0
  assert test_sumstats["_mWT"] == pytest.approx(expected_WT)



def test_single_sample_sumstats():
  test_nchr = 1
  test_chr_end = [400]
  test_w_s = 50
  # 4 individuals (columns), 20 loci (rows)
  test_ga = allel.GenotypeArray([[[ 0, 1],[ 0, 1],[ 1, 1],[ 0, 0]], #  0
                                 [[ 0, 1],[ 0, 1],[ 0, 1],[ 1, 1]], #  1
                                 [[ 0, 1],[ 0, 1],[-1,-1],[ 1, 1]], #  2
                                 [[ 0, 1],[ 0, 0],[ 0, 1],[ 1, 1]], #  3
                                 [[ 0, 1],[ 0, 1],[ 0, 0],[ 1, 1]], #  4
                                 [[ 0, 1],[ 0, 1],[ 0, 0],[ 1, 1]], #  5
                                 [[ 1, 1],[ 0, 1],[ 0, 1],[ 1, 1]], #  6
                                 [[ 0, 1],[ 0, 1],[ 0, 0],[ 1, 1]], #  7
                                 [[ 1, 1],[ 0, 1],[ 0, 0],[ 0, 1]], #  8
                                 [[ 0, 1],[ 0, 1],[ 0, 0],[ 1, 1]], #  9
                                 [[ 1, 1],[ 0, 1],[ 0, 1],[ 1, 1]], # 10
                                 [[ 0, 1],[ 0, 1],[ 0, 0],[ 1, 1]], # 11
                                 [[ 0, 0],[ 0, 1],[ 0, 1],[ 0, 1]], # 12
                                 [[-1,-1],[-1,-1],[-1,-1],[ 1, 1]], # 13
                                 [[ 1, 1],[ 1, 1],[ 0, 0],[ 1, 1]], # 14
                                 [[ 1, 1],[ 1, 1],[ 0, 0],[ 1, 1]], # 15
                                 [[ 1, 1],[ 0, 1],[ 0, 0],[ 0, 1]], # 16
                                 [[ 1, 1],[ 1, 1],[ 1, 0],[ 1, 1]], # 17
                                 [[ 1, 1],[ 1, 1],[ 0, 0],[ 1, 1]], # 18
                                 [[ 0, 1],[ 0, 0],[ 0, 1],[-1,-1]]],# 19
                                 dtype='i1')
  #           0  1  2  3  4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19
  test_pos = (4,10,50,77,99,123,134,150,178,201,209,234,256,270,299,311,315,340,358,378)
  test_sumstats = {}
  test_name = "test"
  single_sample_sumstats(test_ga, test_pos, test_nchr, test_chr_end, test_w_s, test_sumstats, test_name)
  assert test_sumstats["test_maxTD"] == pytest.approx(1.71931236)
  #assert (test_sumstats["test_RoHD"] == [3,8,0]).all()
  assert test_sumstats["test_RoHD_0"] == 3
  assert test_sumstats["test_RoHD_1"] == 8
  assert test_sumstats["test_RoHD_2"] == 0
  test_sumstats = {}
  single_sample_sumstats(test_ga, test_pos, test_nchr, test_chr_end, 390, test_sumstats, test_name)
  #assert (test_sumstats["test_RoHD"] == [3,25,1,0]).all()
  assert test_sumstats["test_RoHD_0"] == 3
  assert test_sumstats["test_RoHD_1"] == 25
  assert test_sumstats["test_RoHD_2"] == 1
  assert test_sumstats["test_RoHD_3"] == 0



### end CALCULATE SUMMARY STATISTICS : SINGLE SAMPLE  ····························


### CALCULATE SUMMARY STATISTICS : ROH  ····························

def roh(gv,pos,missing_threshold=3):
  num_of_alt = gv.to_n_alt(fill=-1)
  missing_sites = np.where(num_of_alt==-1)[0]
  if missing_threshold > np.size(missing_sites) :
    het_sites = np.where(num_of_alt==1)[0]
    roh_limits = list(map(pos.__getitem__, list(het_sites)))
    roh = np.diff(roh_limits)
  else :
    roh = np.array([])
  return roh

def test_roh():
  # one individual, 11 loci
  test_ga = allel.GenotypeArray([[[0,1]],
                                 [[0,1]],
                                 [[0,0]],
                                 [[0,1]],
                                 [[0,1]],
                                 [[0,0]],
                                 [[0,1]],
                                 [[0,1]],
                                 [[1,1]],
                                 [[0,1]],
                                 [[0,1]]], dtype='i1')
  test_pos = (5,10,15,20,100,200,300,1300,3000,10000,20000)
  test_roh = roh(test_ga,test_pos)
  assert (test_roh == [5,10,80,200,1000,8700,10000]).all()
  test_ga = allel.GenotypeArray([[[ 0, 1]],
                                 [[ 0, 1]],
                                 [[ 0, 0]],
                                 [[-1,-1]],
                                 [[ 0, 1]],
                                 [[ 0, 0]],
                                 [[ 0, 1]],
                                 [[-1,-1]],
                                 [[-1,-1]],
                                 [[ 0, 1]],
                                 [[ 0, 1]]], dtype='i1')
  test_roh = roh(test_ga,test_pos)
  assert np.size(test_roh) == 0
  test_roh = roh(test_ga,test_pos,4)
  assert (test_roh == [5,90,200,9700,10000]).all()


def distribution_roh(ga, pos, w_start, w_stop, number_of_bins):
  # w_start and w_stop are the limits of the windows as positions in the genotype array, not positions in the genome
  if number_of_bins>1:
    roh_distribution = np.full(number_of_bins, 0)
  else:
    msg = "Negative value or zero. Number of bins has to be a positive integer"
    raise ValueError(msg)
  for ind in range(0,ga.n_samples):
    gv = ga[w_start:w_stop,ind]
    roh_lenghts = roh(gv, pos[w_start:w_stop])
    for roh_size in range(0,number_of_bins):
      roh_distribution[roh_size] += sum( roh_lenghts >= 10**(roh_size) ) - sum( roh_lenghts >= 10**(roh_size+1) )
  if (roh_distribution < 0).any():
    msg = "Negative value. Number of observations of RoH bin sizes has to be zero or higher"
    raise ValueError(msg)
  return roh_distribution

def test_distribution_roh():
  # 3 individuals, 11 loci
  test_ga = allel.GenotypeArray([[[ 0, 1],[ 0, 1],[ 0, 0]], #  0
                                 [[ 0, 1],[ 0, 1],[ 0, 1]], #  1
                                 [[ 0, 0],[-1,-1],[ 1, 1]], #  2
                                 [[ 0, 1],[ 1, 1],[ 0, 1]], #  3
                                 [[ 0, 1],[ 0, 1],[ 0, 1]], #  4
                                 [[ 0, 0],[ 0, 0],[-1,-1]], #  5
                                 [[ 0, 1],[ 0, 1],[ 0, 1]], #  6
                                 [[ 0, 1],[ 0, 1],[ 0, 1]], #  7
                                 [[ 1, 1],[ 0, 0],[ 0, 1]], #  8
                                 [[ 0, 1],[ 0, 1],[ 0, 1]], #  9
                                 [[ 0, 1],[ 0, 1],[ 0, 0]]],# 10
                                 dtype='i1')
  #            0   1   2   3    4    5    6     7     8      9     10
  test_pos = [ 5, 10, 15, 20, 100, 200, 300, 1300, 3000, 10000, 20000]
  test_start = 0
  test_end = 11
  number_of_bins = 5
  d_roh = distribution_roh(test_ga, test_pos, test_start, test_end, number_of_bins)
  assert (d_roh == [2,5,3,7,2]).all()


def windowed_distribution_roh(ga, pos, size, start, stop):
  # start and stop are the limts (in bp) of the genome whre the RoH are computed
  number_of_bins = int(round(np.log10(size)))+1
  roh_distribution = np.full(number_of_bins, 0)
  # verify that positions are monotonically increasing
  if np.all(pos[1:] > pos[:-1]):
    windows = allel.position_windows(pos, size=size, start=start, stop=stop, step=size)
    locs = allel.window_locations(pos, windows)
    for window_start, window_stop in locs :
      roh_distribution = roh_distribution + distribution_roh(ga, pos, window_start, window_stop, number_of_bins)
  else:
    msg = "Wrong order. Vector of positions has to be monotonically increasing"
    raise ValueError(msg)
  return roh_distribution

def test_windowed_distribution_roh():
  # 3 individuals, 11 loci
  test_ga = allel.GenotypeArray([[[ 0, 1],[ 0, 1],[ 0, 0]], #  0
                                 [[ 0, 1],[ 0, 1],[ 0, 1]], #  1
                                 [[ 0, 0],[-1,-1],[ 1, 1]], #  2
                                 [[ 0, 1],[ 1, 1],[ 0, 1]], #  3
                                 [[ 0, 1],[ 0, 1],[ 0, 1]], #  4
                                 [[ 0, 0],[ 0, 0],[-1,-1]], #  5
                                 [[ 0, 1],[ 0, 1],[ 0, 1]], #  6
                                 [[ 0, 1],[ 0, 1],[ 0, 1]], #  7
                                 [[ 1, 1],[ 0, 0],[ 0, 1]], #  8
                                 [[ 0, 1],[ 0, 1],[ 0, 1]], #  9
                                 [[ 0, 1],[ 0, 1],[ 0, 0]]],# 10
                                 dtype='i1')
  #            0   1   2   3    4    5    6     7     8      9     10
  test_pos = [ 5, 10, 15, 20, 100, 200, 300, 1300, 3000, 10000, 20000]
  test_size = 5000
  d_roh = windowed_distribution_roh(test_ga, test_pos, test_size, 1, 20000)
  assert (d_roh == [2,5,3,4,0]).all()

### CALCULATE SUMMARY STATISTICS : TWO SAMPLE  ····························
def two_samples_sumstats(ga,pair_of_groups,pos,nchr,chr_ends,w_size,sumstats,name="",sep="_"):
  a, b, c = allel.weir_cockerham_fst(g       = ga,
                                     subpops = pair_of_groups)
  fst = np.sum(a) / (np.sum(a) + np.sum(b) + np.sum(c))
  sumstats[name+sep+"Fst"]=fst
  fst_per_variant = (np.sum(a, axis=1) / (np.sum(a, axis=1) + np.sum(b, axis=1) + np.sum(c, axis=1)))
  moments_fst_per_variant = st.describe(fst_per_variant, nan_policy='omit')
  save_moments_2_dict(moments_fst_per_variant,sumstats,name,sep,"l_Fst")
  fst_per_window, _, _ = allel.windowed_weir_cockerham_fst(pos     = pos,
                                                           g       = ga,
                                                           subpops = pair_of_groups,
                                                           size    = w_size,
                                                           start   = 1,
                                                           stop    = chr_ends[0])
  for chromo in range(1,nchr):
    temp, _, _ = allel.windowed_weir_cockerham_fst(pos     = pos,
                                                   g       = ga,
                                                   subpops = pair_of_groups,
                                                   size    = w_size,
                                                   start   = 1+chr_ends[chromo-1],
                                                   stop    = chr_ends[chromo])
    np.append(fst_per_window,temp)
  moments_fst_per_window = st.describe(fst_per_window)
  save_moments_2_dict(moments_fst_per_window,sumstats,name,sep,"w_Fst")
  return

