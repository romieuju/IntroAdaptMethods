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

#Script used to calculate the summary statistics used by MaLAdapt (and other summary statistics). The functions for the various summary statistics can be found in sim_summary_stat_fin.py.  
#Summary statistics are calculated from genotypic and phenotypic matrices in scikit-allel package format, obtained via tree sequences. 
#Most of the functions used to calculate the summary statistics in this script are functions derived from those in the scikit allel packages and the MaLadapt method scripts (Zhang et al., 2023).  

#Import python package and script :
import os, allel, sys
import pandas as pd
import numpy as np
from introadapt import *
from sum_stat_function import *
from scipy.stats import kurtosis
from scipy.stats import skew
import time
import tskit

sum_stat_start = time.time()

# get options for project and simulation:
options = get_options(proj_options_file = sys.argv[1], sim_options_file = sys.argv[2])

#Save the results_dir path :
project_dir = options["results_dir"]+"/"+options["analysis"]+"/"+options["project"]
job_dir     = project_dir+"/"+options["analysis"]+"_"+options["project"]+"_sim_"+options["sim"]

#input and output files :
treeseq            = tskit.load(job_dir+"/Mut_TreeSeq_"+options["analysis"]+"_"+options["project"]+"_sim_"+options["sim"]+".trees")
outputname         = job_dir+"/Sum_Stat_Mut_TreeSeq_"+options["analysis"]+"_"+options["project"]+"_sim_"+options["sim"]+".csv"
csv_file_name      = job_dir+"/sum_stat_"+options["analysis"]+"_"+options["project"]+"_sim_"+options["sim"]+".csv"
duration_file_name = job_dir+"/duration_"+options["analysis"]+"_"+options["project"]+"_sim_"+options["sim"]+".csv"

#Save samples size for each population 
n_popA1 = options["n_samples"][0]
n_popA2 = options["n_samples"][1]
n_popB1 = options["n_samples"][2]
n_popB2 = options["n_samples"][3]
n_popO  = options["n_samples"][4]

#Save sampling time of the donor pop
sample_time_d = options["sample_time_d"]
#Save end time of forward (slim) generation
g_end_sim     = options["generations_forward"]
#Backward generation between end of sim and don pop sampling :
g_samp_don    = g_end_sim - sample_time_d

#idem for donor sister pop :
sample_time_do_s = options["sample_time_d_s"]
g_samp_don_s     = g_end_sim - sample_time_do_s

#Summary statistique parameters:
chr_lenght     = options['window_end']
size_window    = options["window_size"]
start_wind     = options["window_start"]
end_wind       = options["window_end"]
step_between_w = options["window_step"]

if end_wind == 0 :
    end_wind       = chr_lenght
if step_between_w == 0 :
    step_between_w = None 

zns_calculation = options["zns"]
donor_pop_stat  = options["donor_pop_stat"]
s_size          = options["sample_size"]
is_phased       = options['is_phased']
sum_stat_string = []



#Save variants positions in list 
positions = [variant.position for variant in treeseq.variants()]

# Export the genotype data to allel. Unfortunately there's a slight mismatch in the 
# terminology here where genotypes and haplotypes mean different things in the two
# libraries.
tree_genotype_matrix = treeseq.genotype_matrix()
#tree seq genotype matrix to scikit haplotype matrix :  
haplo_all            = allel.HaplotypeArray(tree_genotype_matrix)

if len(tree_genotype_matrix) == 0 :
    counts_dxy, windows_dxy = allel.windowed_count(positions, size=size_window, start=start_wind, stop=end_wind, step=step_between_w)
    nbr_stat                = ["S_rec","S_rec_sister","S_don","exp_het_rec","exp_het_rec_sister","exp_het_don","pi_rec", "pi_rec_sister", "pi_don", "h_theta_rec", "h_theta_rec_sister", "h_theta_don", "w_theta_rec", "w_theta_rec_sister", "w_theta_don", "dxy", "dxsy", "rnd", "rd", "w_c_fst_don_rec", "w_c_fst_don_recsis", "w_c_fst_recsis_rec", "patterson_f3", "patterson_f4", "green_d", "fd", "Q_10_100_q95","Q_10_100_q90","Q_10_100_max","Q_1_100_q95","Q_1_100_q90","Q_1_100_max", "U_10_0_100", "U_10_20_100", "U_10_50_100", "U_10_80_100", "U_1_0_100", "U_1_20_100", "U_1_50_100", "U_1_80_100","h1_rec", "h2_rec", "h12_rec", "h2_h1_rec","h1_rec_sister", "h2_rec_sister", "h12_rec_sister", "h2_h1_rec_sister","h1_don", "h2_don", "h12_don", "h2_h1_don", "zns_rec", "zns_rec_sister", "zns_do" ]
    nbr_stat                = len(nbr_stat)
    nbre_wind               = len(windows_dxy)
    sum_stat                = [[0.0]*nbre_wind]*nbr_stat
else :
    #Allel count of each site
    allel_count_all   = haplo_all.count_alleles()
    # keep only biallelic variant
    biall_filter      = allel_count_all.max_allele()<=1
    haplo_isbiall_all = haplo_all.subset(biall_filter)
    biall_positions   = np.array(positions)[biall_filter].tolist()
    ################################################################################
    #all population genotype :
    geno_isbiall_all               = haplo_isbiall_all.to_genotypes(ploidy=2)
    if is_phased == 'Yes' :
        geno_isbiall_all.is_phased = np.ones(geno_isbiall_all.shape[:-1])

    #Position list in integer :
    int_position = [int(x) for x in biall_positions]

    #Nodes for each population :
    nodes_popO  = [x.id for x in treeseq.nodes() if ((x.time == 1.0) and (x.population == 1))]
    nodes_popA1 = [x.id for x in treeseq.nodes() if ((x.time == g_samp_don_s) and (x.population == 5))]
    nodes_popA2 = [x.id for x in treeseq.nodes() if ((x.time == g_samp_don) and (x.population == 6))]
    nodes_popB1 = [x.id for x in treeseq.nodes() if ((x.time == 1.0) and (x.population == 7))]
    nodes_popB2 = [x.id for x in treeseq.nodes() if ((x.time == 1.0) and (x.population == 8))]

    #Scikit allel haplotype matrix for each pop :
    haplo_popO  = haplo_isbiall_all[:,nodes_popO]
    haplo_popA1 = haplo_isbiall_all[:,nodes_popA1]
    haplo_popA2 = haplo_isbiall_all[:,nodes_popA2]
    haplo_popB1 = haplo_isbiall_all[:,nodes_popB1]
    haplo_popB2 = haplo_isbiall_all[:,nodes_popB2]
    
    #Size range by pop : 
    range_size_popO  = list(set([int(i/2) for i in nodes_popO]))
    range_size_popA1 = list(set([int(i/2) for i in nodes_popA1]))
    range_size_popA2 = list(set([int(i/2) for i in nodes_popA2]))
    range_size_popB1 = list(set([int(i/2) for i in nodes_popB1]))
    range_size_popB2 = list(set([int(i/2) for i in nodes_popB2]))

    #Define donor pop and sister donor pop :
    if options["donor"] == '"A1"' :
        haplo_donor             = haplo_popA1
        haplo_donor_sister      = haplo_popA2
        range_size_donor        = range_size_popA1
        range_size_donor_sister = range_size_popA2
        n_don                   = n_popA1
        n_don_sis=n_popA2
    elif options["donor"] == '"A2"' :
        haplo_donor             = haplo_popA2
        haplo_donor_sister      = haplo_popA1
        range_size_donor        = range_size_popA2
        range_size_donor_sister = range_size_popA1
        n_don                   = n_popA2
        n_don_sis               = n_popA1
    elif options["donor"] == '"B1"' :
        haplo_donor             = haplo_popB1
        haplo_donor_sister      = haplo_popB2
        range_size_donor        = range_size_popB1
        range_size_donor_sister = range_size_popB2
        n_don                   = n_popB1
        n_don_sis               = n_popB2
    elif options["donor"] == '"B2"' :
        haplo_donor             = haplo_popB2
        haplo_donor_sister      = haplo_popB1
        range_size_donor        = range_size_popB2
        range_size_donor_sister = range_size_popB1
        n_don                   = n_popB2
        n_don_sis               = n_popB1

    #Define recipient pop and sister recipient pop :
    if options["recipient"] == '"A1"' :
        haplo_recipient             = haplo_popA1
        range_size_recipient        = range_size_popA1
        n_rec                       = n_popA1
        haplo_recipient_sister      = haplo_popA2
        range_size_recipient_sister = range_size_popA2
        n_rec_sis                   = n_popA2
    elif options["recipient"] == '"A2"' :
        haplo_recipient = haplo_popA2
        range_size_recipient        = range_size_popA2
        n_rec                       = n_popA2
        haplo_recipient_sister      = haplo_popA1
        range_size_recipient_sister = range_size_popA1
        n_rec_sis                   = n_popA1
    elif options["recipient"] == '"B1"' :
        haplo_recipient            = haplo_popB1
        range_size_recipient        = range_size_popB1
        n_rec                       = n_popB1
        haplo_recipient_sister      = haplo_popB2
        range_size_recipient_sister = range_size_popB2
        n_rec_sis                   = n_popB2
    elif options["recipient"] == '"B2"' :
        haplo_recipient             = haplo_popB2
        range_size_recipient        = range_size_popB2
        n_rec                       = n_popB2
        haplo_recipient_sister      = haplo_popB1
        range_size_recipient_sister = range_size_popB1
        n_rec_sis                   = n_popB1

    #Define outgroup population :
    haplo_outgroup      = haplo_popO
    range_size_outgroup = range_size_popO

    #Use donor sister population for summary stat calculation :
    if donor_pop_stat != options["donor"] :
        haplo_donor      = haplo_donor_sister
        range_size_donor = range_size_donor_sister


    #Genotype for each population :
    geno_recipient        = haplo_recipient.to_genotypes(ploidy=2)
    geno_recipient_sister = haplo_recipient_sister.to_genotypes(ploidy=2)
    geno_donor            = haplo_donor.to_genotypes(ploidy=2)
    geno_outgroup         = haplo_outgroup.to_genotypes(ploidy=2)
    if is_phased == 'Yes' :
        geno_recipient.is_phased        = np.ones(geno_recipient.shape[:-1])
        geno_recipient_sister.is_phased = np.ones(geno_recipient_sister.shape[:-1])
        geno_donor.is_phased            = np.ones(geno_donor.shape[:-1])
        geno_outgroup.is_phased         = np.ones(geno_outgroup.shape[:-1])


    #Allel count for each population :
    allel_count_donor            = geno_donor.count_alleles()
    allel_count_recipient        = haplo_recipient.count_alleles()
    allel_count_recipient_sister = geno_recipient_sister.count_alleles()
    allel_count_outgroup         = geno_outgroup.count_alleles()


    #Allel frequency :
    allel_freq_donor            = allel_count_donor.to_frequencies()
    allel_freq_recipient        = allel_count_recipient.to_frequencies()
    allel_freq_recipient_sister = allel_count_recipient_sister.to_frequencies()
    allel_freq_outgroup         = allel_count_outgroup.to_frequencies()

    ##############################################################
    ###################Summary statistics########################
    ##############################################################

    #list with all summary statistics :
    sum_stat = []

    ##Diversity :
    #Segregating site 
    S_don        = windowed_count_segregating(int_position, allel_count_donor, size=size_window, start=start_wind, stop=end_wind, step=step_between_w)
    S_rec        = windowed_count_segregating(int_position, allel_count_recipient, size=size_window, start=start_wind, stop=end_wind, step=step_between_w)
    S_rec_sister = windowed_count_segregating(int_position, allel_count_recipient_sister, size=size_window, start=start_wind, stop=end_wind, step=step_between_w)

    #Expected heterozygotity :
    exp_het_don, windows_exp_het_don, n_bases_exp_het_don, counts_expt_het_don                             = allel.windowed_diversity(int_position,allel_count_donor, size=size_window, start=start_wind, stop=end_wind, step=step_between_w)
    exp_het_don                                                                                            = exp_het_don*(2*n_don-1)/(2*n_don)
    exp_het_rec, windows_exp_het_rec, n_bases_exp_het_rec, counts_expt_het_rec                             = allel.windowed_diversity(int_position,allel_count_recipient, size=size_window, start=start_wind, stop=end_wind, step=step_between_w)
    exp_het_rec                                                                                            = exp_het_rec*(2*n_rec-1)/(2*n_rec)
    exp_het_rec_sister, windows_exp_het_rec_sister, n_bases_exp_het_rec_sister, counts_expt_het_rec_sister = allel.windowed_diversity(int_position,allel_count_recipient_sister, size=size_window, start=start_wind, stop=end_wind, step=step_between_w)
    exp_het_rec                                                                                            = exp_het_rec*(2*n_rec-1)/(2*n_rec)

    # Nucleotide diversity (Pi)
    #start time pi :
    t_pi_start = time.time()
    #
    pi_rec, windows_pi_rec, n_bases_pi_rec, counts_pi_rec                             = allel.windowed_diversity(int_position, allel_count_recipient, size=size_window, start=start_wind, stop=end_wind, step=step_between_w, is_accessible=None, fill=np.nan)
    pi_rec                                                                            = pi_rec*n_bases_pi_rec
    pi_don, windows_pi_don, n_bases_pi_don, counts_pi_don                             = allel.windowed_diversity(int_position, allel_count_donor, size=size_window, start=start_wind, stop=end_wind, step=step_between_w, is_accessible=None, fill=np.nan)
    pi_don                                                                            = pi_don*n_bases_pi_don
    pi_rec_sister, windows_pi_rec_sister, n_bases_pi_rec_sister, counts_pi_rec_sister = allel.windowed_diversity(int_position, allel_count_recipient_sister, size=size_window, start=start_wind, stop=end_wind, step=step_between_w, is_accessible=None, fill=np.nan)
    pi_rec_sister                                                                     = pi_rec_sister*n_bases_pi_rec_sister
    # end time pi :
    t_pi_end = time.time()
    #duration pi
    duration_pi = t_pi_end - t_pi_start
    print("(summary stat) theta pi duration : "+str(duration_pi)+" s ")

    #Fay and Wu theta (theta H) :
    t_H_start = time.time()
    h_theta_rec, h_theta_windows_rec, h_theta_counts_rec                      = windowed_thetah_mal(int_position, allel_count_recipient, size=size_window, start=start_wind, stop=end_wind, step=step_between_w, is_accessible=None, fill=np.nan)
    h_theta_rec_sister, h_theta_windows_rec_sister, h_theta_counts_rec_sister = windowed_thetah_mal(int_position, allel_count_recipient_sister, size=size_window, start=start_wind, stop=end_wind, step=step_between_w, is_accessible=None, fill=np.nan)
    h_theta_don, h_theta_windows_don, h_theta_counts_don                      = windowed_thetah_mal(int_position, allel_count_donor, size=size_window, start=start_wind, stop=end_wind, step=step_between_w, is_accessible=None, fill=np.nan)
    t_H_end = time.time()
    duration_H = t_H_end - t_H_start
    print("(summary stat) theta H duration : "+str(duration_H)+" s ")

    #Watterson theta (theta w) :
    t_W_start = time.time()
    w_theta_rec, windows_w_theta_rec, n_bases_w_theta_rec, counts_w_theta_rec                             = allel.windowed_watterson_theta(int_position, allel_count_recipient, size=size_window, start=start_wind, stop=end_wind, step=step_between_w)
    w_theta_rec                                                                                           = w_theta_rec*n_bases_w_theta_rec
    w_theta_rec_sister, windows_w_theta_rec_sister, n_bases_w_theta_rec_sister, counts_w_theta_rec_sister = allel.windowed_watterson_theta(int_position, allel_count_recipient_sister, size=size_window, start=start_wind, stop=end_wind, step=step_between_w)
    w_theta_rec_sister                                                                                    = w_theta_rec_sister*n_bases_w_theta_rec_sister
    w_theta_don, windows_w_theta_don, n_bases_w_theta_don, counts_w_theta_don                             = allel.windowed_watterson_theta(int_position, allel_count_donor, size=size_window, start=start_wind, stop=end_wind, step=step_between_w)
    w_theta_don                                                                                           = w_theta_don*n_bases_w_theta_don
    t_W_end = time.time()
    duration_W = t_W_end - t_W_start
    print("(summary stat) theta W duration : "+str(duration_W)+" s ")

    ##Population divergence :
    #Dxy, x = recipient , s = donor (divergence between two taxon)
    t_dxy_start = time.time()
    dxy, windows_dxy, n_bases_dxy, counts_dxy     = allel.windowed_divergence(int_position, allel_count_recipient, allel_count_donor, size=size_window, start=start_wind, stop=end_wind, step=step_between_w)
    dxy                                           = dxy*n_bases_dxy
    #Dxsy, x = recipient sister, y = donor : 
    dxsy, windows_dxsy, n_bases_dxsy, counts_dxsy = allel.windowed_divergence(int_position, allel_count_recipient_sister, allel_count_donor, size=size_window, start=start_wind, stop=end_wind, step=step_between_w)
    dxsy                                          = dxsy*n_bases_dxsy
    #Dxo, x = recipient sister, y = outgroup : 
    dxo, windows_dxo, n_bases_dxo, counts_dxo     = allel.windowed_divergence(int_position, allel_count_recipient, allel_count_outgroup, size=size_window, start=start_wind, stop=end_wind, step=step_between_w)
    dxo = dxo*n_bases_dxo
    #Dyo, x = donor, y = outgroup : 
    dyo, windows_dyo, n_bases_dyo, counts_dyo     = allel.windowed_divergence(int_position, allel_count_donor, allel_count_outgroup, size=size_window, start=start_wind, stop=end_wind, step=step_between_w)
    dyo                                           = dyo*n_bases_dyo
    dout                                          = (dxo + dyo)/2
    #RND
    rnd                                           = dxy/dout
    #Racimo Rd (naive way) 
    rd                                            = dxy/dxsy
    t_dxy_end    = time.time()
    duration_dxy = t_dxy_end - t_dxy_start
    print("(summary stat) divergence stat (dxy) duration : "+str(duration_dxy)+" s ")

    ##Admixture :
    #f3 Patterson (2012) (ac recipient, ac donor, ac recipient sister) (test population, first source population, second source population)
    t_f3_start   = time.time()
    patterson_f3 = windowed_patterson_f3(int_position, allel_count_recipient, allel_count_donor, allel_count_recipient_sister, size=size_window, start=start_wind, stop=end_wind, step=step_between_w, normed=True)
    t_f3_end     = time.time()
    duration_f3  = t_f3_end - t_f3_start
    print("(summary stat) f3 patterson duration : "+str(duration_f3)+" s ")

    #D-patterson (f-4) : 
    t_f4_start   = time.time()
    patterson_f4 = windowed_patterson_f4(int_position, allel_count_outgroup, allel_count_donor,allel_count_recipient_sister, allel_count_recipient, size=size_window, start=start_wind, stop=end_wind, step=step_between_w, normed=True)
    t_f4_end     = time.time()
    duration_f4  = t_f4_end - t_f4_start
    print("(summary stat) f4 patterson duration : "+str(duration_f4)+" s ")

    #D_Green or ABBA BABA test 
    G_D_start    = time.time()
    green_d      = windowed_Green_d(int_position, allel_count_recipient_sister, allel_count_recipient, allel_count_donor, size=size_window, start=start_wind, stop=end_wind, step=step_between_w)
    G_D_end      = time.time()
    duration_G_D = G_D_end - G_D_start
    print("(summary stat) Green D duration : "+str(duration_G_D)+" s ")

    #fD (Martin et al. 2015)
    fd_start     = time.time()
    fd           = windowed_f_d(int_position, allel_count_recipient_sister, allel_count_recipient, allel_count_donor, size=size_window, start=start_wind, stop=end_wind, step=step_between_w)
    fd_end       = time.time()
    duration_fd  = fd_end - fd_start
    print("(summary stat) fd duration : "+str(duration_fd)+" s ")


    #Q's Racimo :
    Q_start    = time.time()
    Q_10_100_q95, Q_10_100_q90, Q_10_100_max, Q_1_100_q95, Q_1_100_q90, Q_1_100_max = windowed_q_racimo(int_position, allel_count_donor, allel_count_recipient_sister, allel_count_recipient, size=size_window, start=start_wind, stop=end_wind, step=step_between_w)
    Q_end      = time.time()
    duration_Q = Q_end - Q_start
    print("(summary stat) Racimo's Q duration : "+str(duration_Q)+" s ")
    

    #U's Racimo :
    U_start    = time.time()
    U_10_0_100, U_10_20_100, U_10_50_100, U_10_80_100, U_1_0_100, U_1_20_100, U_1_50_100, U_1_80_100 = windowed_u_racimo(int_position, allel_count_donor, allel_count_recipient_sister, allel_count_recipient, size=size_window, start=start_wind, stop=end_wind, step=step_between_w)
    U_end      = time.time()
    duration_U =  U_end - U_start
    print("(summary stat) Racimo's U duration : "+str(duration_U)+" s ")

    #F-statistics
    #Weir and Cockerham's Fst :
    fst_start = time.time()
    w_c_fst_don_rec, window_w_c_fst_don_rec, counts_w_c_fst_don_rec          = allel.windowed_weir_cockerham_fst(int_position, geno_isbiall_all,[range_size_donor, range_size_recipient], size=size_window, start=start_wind, stop=end_wind, step=step_between_w)
    w_c_fst_recsis_rec, window_w_c_fst_recsis_rec, counts_w_c_fst_recsis_rec = allel.windowed_weir_cockerham_fst(int_position, geno_isbiall_all,[range_size_recipient_sister, range_size_recipient], size=size_window, start=start_wind, stop=end_wind, step=step_between_w)
    w_c_fst_don_recsis, window_w_c_fst_don_recsis, counts_w_c_fst_don_recsis = allel.windowed_weir_cockerham_fst(int_position, geno_isbiall_all,[range_size_recipient_sister, range_size_donor], size=size_window, start=start_wind, stop=end_wind, step=step_between_w)
    fst_end = time.time()
    duration_fst = fst_end - fst_start
    print("(summary stat) Fst duration : "+str(duration_fst)+" s ")

    #test Rd Jules/Raph version :
    rd_jules_start = time.time()
    rd_jules_window = rd_window_jules(int_position, allel_count_donor, allel_count_recipient, allel_count_recipient_sister, size=size_window, start=start_wind, stop=end_wind, step=None)
    rd_jules_end = time.time()
    duration_rd_jules = rd_jules_end - rd_jules_start
    print("(summary stat) Racimo's Rd (jules version) duration : "+str(duration_rd_jules)+" s ")

    #Save stat and stat names in list
    sum_stat        = [S_rec, S_rec_sister, S_don, exp_het_rec, exp_het_rec_sister, exp_het_don, pi_rec, pi_rec_sister, pi_don, h_theta_rec, h_theta_rec_sister, h_theta_don, w_theta_rec, w_theta_rec_sister, w_theta_don, dxy, dxsy, rnd, w_c_fst_don_rec, w_c_fst_don_recsis, w_c_fst_recsis_rec, patterson_f3, patterson_f4, green_d, fd, Q_10_100_q95,Q_10_100_q90,Q_10_100_max,Q_1_100_q95,Q_1_100_q90,Q_1_100_max, U_10_0_100, U_10_20_100, U_10_50_100, U_10_80_100, U_1_0_100, U_1_20_100, U_1_50_100, U_1_80_100, rd_jules_window]
    sum_stat_string = ["S_rec","S_rec_sister","S_don","exp_het_rec","exp_het_rec_sister","exp_het_don","pi_rec", "pi_rec_sister", "pi_don", "h_theta_rec", "h_theta_rec_sister", "h_theta_don", "w_theta_rec", "w_theta_rec_sister", "w_theta_don", "dxy", "dxsy", "rnd", "w_c_fst_don_rec", "w_c_fst_don_recsis", "w_c_fst_recsis_rec", "patterson_f3", "patterson_f4", "green_d", "fd", "Q_10_100_q95","Q_10_100_q90","Q_10_100_max","Q_1_100_q95","Q_1_100_q90","Q_1_100_max", "U_10_0_100", "U_10_20_100", "U_10_50_100", "U_10_80_100", "U_1_0_100", "U_1_20_100", "U_1_50_100", "U_1_80_100", "rd_jules"]
    
    duration_stat        = [duration_pi, duration_H, duration_W, duration_dxy, duration_f3, duration_f4, duration_G_D, duration_Q, duration_U, duration_fst, duration_rd_jules]
    duration_stat_string = ["duration_pi", "duration_H", "duration_W", "duration_dxy", "duration_f3", "duration_f4", "duration_G_D", "duration_Q", "duration_U", "duration_fst", "duration_rd_jules"]


    ##Summary stat calculated if data are phased (use the haplotype information) : 
    if is_phased == 'Yes' :
        #racimo rd (from MaLAdapt) :
        rd_racimo_start    = time.time()
        rd_window = windowed_rd_racimo(int_position, geno_recipient, geno_donor, geno_recipient_sister, size_window, start_wind, end_wind, w_step = step_between_w)
        rd_racimo_end      = time.time()
        duration_rd_racimo = rd_racimo_end - rd_racimo_start
        print("(summary stat) Racimo's Rd (Maladapt version) duration : "+str(duration_rd_racimo)+" s ")
        
        ##Haplotype homozygosity
        haplotyp_stat_start = time.time()
        h1_rec, h2_rec, h12_rec, h123_rec, h2_h1_rec                                    = windowed_garud_h(int_position, haplo_recipient, size=size_window, start=start_wind, stop=end_wind, step=step_between_w)
        h1_rec_sister, h2_rec_sister, h12_rec_sister, h123_rec_sister, h2_h1_rec_sister = windowed_garud_h(int_position, haplo_recipient_sister, size=size_window, start=start_wind, stop=end_wind, step=step_between_w)
        h1_don, h2_don, h12_don, h123_don, h2_h1_don                                    = windowed_garud_h(int_position, haplo_donor, size=size_window, start=start_wind, stop=end_wind, step=step_between_w)
        h1_out, h2_out, h12_out, h123_out, h2_h1_out                                    = windowed_garud_h(int_position, haplo_outgroup, size=size_window, start=start_wind, stop=end_wind, step=step_between_w)
        haplotyp_stat_end   = time.time()
        duration_hap_stat   = haplotyp_stat_end - haplotyp_stat_start
        print("(summary stat) Haplotype homozy stats duration : "+str(duration_hap_stat)+" s ")

        sum_stat        = sum_stat + [rd_window, h1_rec, h2_rec, h12_rec, h2_h1_rec, h1_rec_sister, h2_rec_sister, h12_rec_sister, h2_h1_rec_sister, h1_don, h2_don, h12_don, h2_h1_don]
        sum_stat_string = sum_stat_string + ["rd", "h1_rec", "h2_rec", "h12_rec", "h2_h1_rec","h1_rec_sister", "h2_rec_sister", "h12_rec_sister", "h2_h1_rec_sister","h1_don", "h2_don", "h12_don", "h2_h1_don"]

        duration_stat        = duration_stat + [duration_rd_racimo, duration_hap_stat]
        duration_stat_string = duration_stat_string + ["duration_rd_racimo", "duration_hap_stat"]

        if zns_calculation == 'Yes' :
            ##DL
            #ZnS's Kelly :
            zns_start = time.time()
            zns_window_rec        = windowed_zns(int_position, geno_recipient, size_window, start_wind, end_wind)
            zns_window_rec_sister = windowed_zns(int_position, geno_recipient_sister, size_window, start_wind, end_wind)
            zns_window_do         = windowed_zns(int_position, geno_donor, size_window, start_wind, end_wind)
            zns_end = time.time()
            duration_zns = zns_end - zns_start

            sum_stat        = sum_stat + [zns_window_rec, zns_window_rec_sister, zns_window_do]
            sum_stat_string = sum_stat_string + ["zns_rec", "zns_rec_sister", "zns_do" ]
            print("(summary stat) Zns duration : "+str(duration_zns)+" s ")

            duration_stat        = duration_stat + [duration_zns]
            duration_stat_string = duration_stat_string + ["duration_zns"]

duration_stat         = [options["sim"]] + duration_stat 
duration_stat_string  = ["sim"] + duration_stat_string
data_duration         = pd.DataFrame(data=duration_stat)
data_duration         = data_duration.transpose()
data_duration.columns = duration_stat_string

if os.path.exists(duration_file_name) :
    data_duration.to_csv(duration_file_name, mode = 'a', header = False, index = False)
else :
    data_duration.to_csv(duration_file_name, mode = 'a', header = True, index = False)

#create a data frame with summary statistics :
start                = np.transpose(windows_dxy)[0].tolist()
end                  = np.transpose(windows_dxy)[1].tolist()
nbr_var              = np.transpose(counts_dxy).tolist()
sum_stat             = sum_stat
sum_stat_dat         = pd.DataFrame(data=sum_stat)
sum_stat_dat         = sum_stat_dat.transpose()
sum_stat_dat.columns = sum_stat_string
sum_stat_dat.insert(0,'start', start)
sum_stat_dat.insert(1,'end', end)
sum_stat_dat.insert(2,'num_var', nbr_var)
sum_stat_dat         = sum_stat_dat.replace(np.inf, np.nan)
sum_stat_dat_final2  = sum_stat_dat.fillna(value="NaN")
sum_stat_dat_final2.to_csv(outputname, index=False)


#Create dataframe with the mean, variance, skevness, kurtosis, min, man, quant 5% and quand 95 % summarised for all windows of a simulation: 
sum_stat_sum_data = calcul_stat_ref_table(sum_stat_dat, options["sim"])
if os.path.exists(csv_file_name) :
    sum_stat_sum_data.to_csv(csv_file_name, mode = 'a', header = False, index = False)
else :
    sum_stat_sum_data.to_csv(csv_file_name, mode = 'a', header = True, index = False)


sum_stat_end      = time.time()
duration_sum_stat = sum_stat_end - sum_stat_start
print("(summary stat) summary stat duration : "+str(duration_sum_stat)+" s ")

