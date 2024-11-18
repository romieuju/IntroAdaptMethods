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

import itertools, allel
import numpy as np
import pandas as pd
from scipy.stats import kurtosis
from scipy.stats import skew

##Script with the functions used to calculate the various summary statistics for the pipeline. These functions have been taken and adapted from the python package scikit-allel and from the scripts present on the git of the MaLAdapt method (Zhang et al., 2023). 
#Scikit-allel : https://scikit-allel.readthedocs.io/en/stable/stats.html
#MaLAdapt git : https://github.com/xzhang-popgen/maladapt

##Parameters names descriptions :
#position_array, pos                                                                 : list or 1D array of variants positions
#geno_array_scikit,                                                                  : genotype array scikit-allel format
#ga_scikit_recipient,                                                                : recipient's genotype array (from scikit allel) 
#ga_scikit_donor,                                                                    : donor's genotype array (from scikit allel)
#ga_sckit_nonintro,                                                                  : non-introgressed's genotype array (from scikit allel)
#h,                                                                                  : haplotype array (from scikit allel)
#hap_donor,                                                                          : donor's haplotype array
#hap_recipient,                                                                      : recipient's haplotype array 
#hap_non_intro_recipient                                                             : non-introgressed haplotype array
#ac, acc, aca, acb, ac_poprecni, ac_poprec, ac_popd                                  : population's allel count (scikit allel); c,a,b : ac from pop ((a,b),c), poprecni, poprec, popd : ac from non-intro, recipient and donor pop. 
#allel_count_donor, allel_count_recipient, allel_count_recipient_sister : Idem donor : donor pop, recipient : recipient pop, recipient sister : recipient sister (non-into)
#w_size, size                                                                        : window size (int or float)
#w_start, start                                                                      : start position to define windows 
#w_stop, stop                                                                        : end position to define windows
#w_step, step                                                                       : for sliding windo, step size between each window



###Racimo's Rd, adapted from Racimo et al., 2017 and Zhang et al., 2023 scripts for MaLAdapt method (https://doi.org/10.1093/sysbio/syad033)
#Adapt to be calculated from scikit-allel matrix
#Input : donor, recipient and recipient sister (non-introgressed) haplotype array (need phased genotype matrix) 
def rd_racimo(hap_donor, hap_recipient, hap_non_intro_recipient):
    divratio = []
    for archi in range(0, hap_donor.shape[0]):
        divarchintro    = np.apply_along_axis(lambda x: sum(np.squeeze(np.asarray(hap_donor[archi,])) != x), 1, hap_recipient )
        divarchintro    = divarchintro.astype("float")
        divarchnonintro = np.apply_along_axis(lambda x: sum(np.squeeze(np.asarray(hap_donor[archi,])) != x), 1, hap_non_intro_recipient)
        divarchnonintro = 1 / divarchnonintro.astype("float")
        for comb in itertools.product(divarchintro,divarchnonintro):
            if comb[1] != 0:
                divratio.append( comb[0]*comb[1] )
    divratioavg = float(sum(divratio)) / float(len(divratio))
    return divratioavg

#Calcul Racimo's Rd on windows 
def windowed_rd_racimo(position_array, ga_scikit_recipient, ga_scikit_donor, ga_sckit_nonintro, w_size, w_start = None, w_stop = None , w_step = None) :
    count_rd, windo_rd = allel.windowed_count(position_array, start= w_start, stop = w_stop, size = w_size, step = w_step)
    #transpose genotype array
    t_ga_popr_scikit   = np.transpose(ga_scikit_recipient)
    t_ga_popd_scikit   = np.transpose(ga_scikit_donor)
    t_ga_popni_scikit  = np.transpose(ga_sckit_nonintro)
    #create variables for new shape
    shape_0_r  = t_ga_popr_scikit.shape[0]*t_ga_popr_scikit.shape[1]
    shape_1_r  = t_ga_popr_scikit.shape[2]
    shape_0_d  = t_ga_popd_scikit.shape[0]*t_ga_popd_scikit.shape[1]
    shape_1_d  = t_ga_popd_scikit.shape[2]
    shape_0_ni = t_ga_popni_scikit.shape[0]*t_ga_popni_scikit.shape[1]
    shape_1_ni = t_ga_popni_scikit.shape[2]
    #reshape 3D genotype array in 2D array : (haplotype array from maladapt)
    twod_ga_popr_scikit  = t_ga_popr_scikit.reshape(shape_0_r,shape_1_r)
    twod_ga_popd_scikit  = t_ga_popd_scikit.reshape(shape_0_d,shape_1_d)
    twod_ga_popni_scikit = t_ga_popni_scikit.reshape(shape_0_ni,shape_1_ni)
    #dtype = float64
    haplo_r    = twod_ga_popr_scikit.astype('float64')
    haplo_d    = twod_ga_popd_scikit.astype('float64')
    haplo_ni_r = twod_ga_popni_scikit.astype('float64')
    Rd_list    = []
    #for each window :
    for w in range(len(windo_rd)):
        #windows :
        #variants present in the window :
        thesepos = [pos for pos in position_array if (pos <= windo_rd[w][1] and pos >= windo_rd[w][0])]
        #variant index 
        if len(thesepos) > 0 : 
            start   = position_array.index(thesepos[0])
            end     = position_array.index(thesepos[len(thesepos)-1])+1
            #window's lenght 
            len_seg = windo_rd[w][1]-windo_rd[w][0]
            #windows haplotypes by pop 
            pop_r_haplo, pop_d_haplo, pop_ni_r_haplo = haplo_r[:,start:end], haplo_d[:,start:end], haplo_ni_r[:,start:end]
            #Racimo's rd calculation
            rd = rd_racimo(pop_d_haplo, pop_r_haplo, pop_ni_r_haplo)
            Rd_list.append(rd)
        else :
            rd = float("NaN")
            Rd_list.append(rd)
    return np.array(Rd_list)

##New way to calculate Racimo's rd from dxy (no need of phased data):
#Need variant position list (pos), donor allel count in scikit-allel format (allel_count_donor), recipient allel count (allel_count_recipient), recipient sister allel count (allel_count_recipient_sister), window size (size), start position to calculate rd (sart), end position to calculate rd (end), if sliding window step between window (step)
def rd_window_jules(pos, allel_count_donor, allel_count_recipient, allel_count_recipient_sister, size, start=None, stop=None, step=None) :
  allel_count_donor = allel_count_donor
  allel_count_recipient = allel_count_recipient
  allel_count_recipient_sister = allel_count_recipient_sister
  # Donor's haplotype number for each site 
  anD = np.sum(allel_count_donor, axis=1)
  # recipient's haplotype number for each site 
  anR = np.sum(allel_count_recipient, axis=1)
  # Possible number of haplotype pairs for each site
  n_pairs_dr = anD * anR
  # Number of pairwise comparisons where there is no difference for each site
  n_same_dr = np.sum(allel_count_donor * allel_count_recipient, axis=1)
  # Number of pairwise differences 
  n_diff_dr = n_pairs_dr - n_same_dr
  #pairwise difference normalized by the number of haplotype in donor pop
  n_diff_dr_sum, windows, counts = allel.windowed_statistic(pos, values=n_diff_dr, statistic=np.sum, size=size, start=start, stop=stop, step=step, fill=0)
  anR_sum, windows, counts = allel.windowed_statistic(pos, values=anR, statistic=np.sum, size=size, start=start, stop=stop, step=step, fill=0)
  new_dxy_bis = n_diff_dr_sum/anR_sum
  # donor's haplotype number for each site 
  anD = np.sum(allel_count_donor, axis=1)
  # recipient's haplotype number for each site 
  anNI = np.sum(allel_count_recipient_sister, axis=1)
  # Possible number of haplotype pairs for each site
  n_pairs_dni = anD * anNI
  # Number of pairwise comparisons where there is no difference for each site
  n_same_dni = np.sum(allel_count_donor * allel_count_recipient_sister, axis=1)
  # Number of pairwise differences 
  n_diff_dni = n_pairs_dni - n_same_dni
  #pairwise difference normalized by the number of haplotype in donor pop
  #new_dxys = sum(n_diff_dni)/sum(anNI)
  n_diff_dni_sum, windows, counts = allel.windowed_statistic(pos, values=n_diff_dni, statistic=np.sum, size=size, start=start, stop=stop, step=step, fill=0)
  anNI_sum, windows, counts = allel.windowed_statistic(pos, values=anNI, statistic=np.sum, size=size, start=start, stop=stop, step=step, fill=0)
  new_dxys_bis = n_diff_dni_sum/anNI_sum
  #Calcul divergence ratio :
  rd_jules=new_dxy_bis/new_dxys_bis
  return rd_jules


###Theta pi
#theta_pi, adapted from allel.mean_pairwise_difference window from scikit-allel:
def windowed_theta_pi(pos, ac, size=None, start=None, stop=None, step=None, windows=None, is_accessible=None, fill=np.nan):
    # assume number of chromosomes sampled is constant for all variants
    n   = ac.sum(axis=1).max()
    # calculate mean pairwise difference
    mpd = allel.mean_pairwise_difference(ac, fill=0)
    theta_pi_abs, windows, counts = allel.windowed_statistic(pos, values=mpd, statistic=np.sum, size=size, start=start, stop=stop, step=step, windows=windows, fill=np.nan)
    # theta per base
    theta_pi, n_bases = allel.per_base(theta_pi_abs, windows=windows, is_accessible=is_accessible, fill=fill)
    return theta_pi, windows, n_bases, counts

###Theta H (Fay and Wu, 2000)
## Adapted from MaLAdapt script
def theta_h_mal(ac) :
    #number of haplotypes
    n = ac.sum(axis=1).max()
    #for each variant number of derived allel
    if ac.shape[1] != 2 :
        thetaH = np.nan 
    else :
        Si = ac[:,1]
        #square number of derived allele*the number of time of this number of derived allel is find (x = number of seq)
        num = [x**2*sum(Si == x) for x in list(range(1,n))]
        #sum of numerator * combination(n,2)
        thetaH = sum(num)*2/(n*(n-1))
        return thetaH

#Windowed Theta H (caluclated from scikit-allel genotype matrix)
def windowed_thetah_mal(pos, ac, size=None, start=None, stop=None, step=None, windows=None, is_accessible=None, fill=np.nan) :
    theta_h_abs, windows, counts = allel.windowed_statistic(pos, values=ac, statistic=theta_h_mal, size=size, start=start, stop=stop, step=step, windows=windows, fill=np.nan)
    return theta_h_abs, windows, counts


##Patterson F3 (Patterson, 2012)
#Taken from scikit-allel (patterson_f3()) :
#acc : recipient pop, aca : donor pop 1, acb :  sister pop pop
def patterson_f3_jules(acc, aca, acb):
    # check inputs
    aca = allel.AlleleCountsArray(aca, copy=False)
    if aca.shape[1] == 1 :
        nbre_variant = aca.shape[0]
        aca = np.append(aca, np.array([[0]]*nbre_variant),axis=1)
        aca = allel.AlleleCountsArray(aca, copy=False)
    acb = allel.AlleleCountsArray(acb, copy=False)
    if acb.shape[1] == 1 :
        nbre_variant = acb.shape[0]
        acb = np.append(acb, np.array([[0]]*nbre_variant),axis=1)
        acb = allel.AlleleCountsArray(acb, copy=False)
    acc = allel.AlleleCountsArray(acc, copy=False)
    if acc.shape[1] == 1 :
        nbre_variant = acc.shape[0]
        acc = np.append(acc, np.array([[0]]*nbre_variant),axis=1)
        acc = allel.AlleleCountsArray(acc, copy=False)
    allel.util.check_dim0_aligned(aca, acb, acc)
    # compute allele number and heterozygosity in test population
    sc = acc.sum(axis=1)
    hc = allel.stats.admixture.h_hat(acc)
    # compute sample frequencies for the alternate allele
    a = aca.to_frequencies()[:, 1]
    b = acb.to_frequencies()[:, 1]
    c = acc.to_frequencies()[:, 1]
    # compute estimator
    T = ((c - a) * (c - b)) - (hc / sc)
    B = 2 * hc
    return T, B


#Windowed f3, taken from scikit-allel (moving_patterson_f3()):
def windowed_patterson_f3 (pos, acc, aca, acb, size, start=None, stop=None, step=None, normed=True) :
    T, B = patterson_f3_jules(acc, aca, acb)
    # calculate value of statistic within each block
    if normed:
        T_bsum, windows_t, counts_t = allel.windowed_statistic(pos, T, statistic=np.nansum, size=size,
                                  start=start, stop=stop, step=step)
        B_bsum, windows_b , counts_b = allel.windowed_statistic(pos, B, statistic=np.nansum, size=size,
                                  start=start, stop=stop, step=step)
        f3 = T_bsum / B_bsum
    else:
        f3, windows_f3, counts_f3 = allel.windowed_statistic(pos, T, statistic=np.nanmean, size=size,
                              start=start, stop=stop, step=step)
    return f3

##Patterson'D 
#Taken from scikit-allel (patterson_d())
def patterson_d_jules(aca, acb, acc, acd):
    # check inputs
    aca = allel.AlleleCountsArray(aca, copy=False)
    if aca.shape[1] == 1 :
        nbre_variant = aca.shape[0]
        aca = np.append(aca, np.array([[0]]*nbre_variant),axis=1)
        aca = allel.AlleleCountsArray(aca, copy=False)
    acb = allel.AlleleCountsArray(acb, copy=False)
    if acb.shape[1] == 1 :
        nbre_variant = acb.shape[0]
        acb = np.append(acb, np.array([[0]]*nbre_variant),axis=1)
        acb = allel.AlleleCountsArray(acb, copy=False)
    acc = allel.AlleleCountsArray(acc, copy=False)
    if acc.shape[1] == 1 :
        nbre_variant = acc.shape[0]
        acc = np.append(acc, np.array([[0]]*nbre_variant),axis=1)
        acc = allel.AlleleCountsArray(acc, copy=False)
    acd = allel.AlleleCountsArray(acd, copy=False)
    if acd.shape[1] == 1 :
        nbre_variant = acd.shape[0]
        acd = np.append(acd, np.array([[0]]*nbre_variant),axis=1)
        acd = allel.AlleleCountsArray(acd, copy=False)
    allel.util.check_dim0_aligned(aca, acb, acc, acd)
    # compute sample frequencies for the alternate allele
    a = aca.to_frequencies()[:, 1]
    b = acb.to_frequencies()[:, 1]
    c = acc.to_frequencies()[:, 1]
    d = acd.to_frequencies()[:, 1]
    # compute estimator
    num = (a - b) * (c - d)
    den = (a + b - (2 * a * b)) * (c + d - (2 * c * d))
    return num, den

#Windowed Patterson D, adapted from scikit-allel (average_patterson_d())
#(f-4)  (A;(B;(C,D))) : A = aca (out or anc), B = acb (don), c = acc (noni), d =acd (rec)
def windowed_patterson_f4 (pos, aca, acb, acc, acd, size, start=None, stop=None, step=None, normed=True) :
    # calculate per-variant values
    num, den = patterson_d_jules(aca, acb, acc, acd)
    # N.B., nans can occur if any of the populations have completely missing
    # genotype calls at a variant (i.e., allele number is zero). Here we
    # assume that is rare enough to be negligible.
    # compute the numerator and denominator within each window
    num_sum, windows_num, counts_num = allel.windowed_statistic(pos, num, statistic=np.nansum, size=size, start=start, stop=stop, step=step)
    den_sum, windows_den, counts_den = allel.windowed_statistic(pos, den, statistic=np.nansum, size=size, start=start, stop=stop, step=step)
    # calculate the statistic values in each block
    d = num_sum / den_sum
    return d

#D-statistics from Zhang et al.,2023 : MaLAdapt:
#Calcul D-Statistic (Green et al., 2010) :
#ac = allel count (ni = sister pop of intro pop, rec = intro pop, d = donor pop with (((pop_ni,pop_rec),pop_d,pop_o)))
def windowed_Green_d (pos, ac_poprecni, ac_poprec, ac_popd, size, start=None, stop=None, step=None):
    #Only biallelic variant control :
    ac_recni = allel.AlleleCountsArray(ac_poprecni, copy=False)
    if ac_recni.shape[1] == 1 :
        nbre_variant = ac_recni.shape[0]
        ac_recni = np.append(ac_recni, np.array([[0]]*nbre_variant),axis=1)
        ac_recni = allel.AlleleCountsArray(ac_recni, copy=False)
    ac_rec = allel.AlleleCountsArray(ac_poprec, copy=False)
    if ac_rec.shape[1] == 1 :
        nbre_variant = ac_rec.shape[0]
        ac_rec = np.append(ac_rec, np.array([[0]]*nbre_variant),axis=1)
        ac_rec = allel.AlleleCountsArray(ac_rec, copy=False)
    ac_don = allel.AlleleCountsArray(ac_popd, copy=False)
    if ac_don.shape[1] == 1 :
        nbre_variant = ac_don.shape[0]
        ac_don = np.append(ac_don, np.array([[0]]*nbre_variant),axis=1)
        ac_don = allel.AlleleCountsArray(ac_don, copy=False)
    allel.util.check_dim0_aligned(ac_recni, ac_rec, ac_don)
    #derived allel frequency : 
    freq_recni = ac_recni.to_frequencies()[:, 1]#p1
    freq_rec = ac_rec.to_frequencies()[:, 1]#p2
    freq_don = ac_don.to_frequencies()[:, 1]#p3
    #ABBA and BABA calcul
    abba_vec = (1.0 - freq_recni)*freq_rec*freq_don
    baba_vec = freq_recni*(1.0 - freq_rec)*freq_don
    abba_sum, windows_abba, counts_abba = allel.windowed_statistic(pos, abba_vec, statistic=np.nansum, size=size, start=start, stop=stop, step=step)
    baba_sum, windows_baba, counts_baba = allel.windowed_statistic(pos, baba_vec, statistic=np.nansum, size=size, start=start, stop=stop, step=step)
    #Green et al, 2010's D-statistic calcul : 
    D_stat = (abba_sum - baba_sum) / (abba_sum + baba_sum)
    return D_stat

# Racimo'Q, Adapted from Racimo et al., 2017 and Zhang et al., 2023 (MaLAdapt)
def windowed_q_racimo(pos, ac_popd, ac_poprecni, ac_poprec, size, start=None, stop=None, step=None):
    ac_don = allel.AlleleCountsArray(ac_popd, copy=False)
    if ac_don.shape[1] == 1 :
        nbre_variant = ac_don.shape[0]
        ac_don = np.append(ac_don, np.array([[0]]*nbre_variant),axis=1)
        ac_don = allel.AlleleCountsArray(ac_don, copy=False)
    ac_recni = allel.AlleleCountsArray(ac_poprecni, copy=False)
    if ac_recni.shape[1] == 1 :
        nbre_variant = ac_recni.shape[0]
        ac_recni = np.append(ac_recni, np.array([[0]]*nbre_variant),axis=1)
        ac_recni = allel.AlleleCountsArray(ac_recni, copy=False)
    ac_rec = allel.AlleleCountsArray(ac_poprec, copy=False)
    if ac_rec.shape[1] == 1 :
        nbre_variant = ac_rec.shape[0]
        ac_rec = np.append(ac_rec, np.array([[0]]*nbre_variant),axis=1)
        ac_rec = allel.AlleleCountsArray(ac_rec, copy=False)
    allel.util.check_dim0_aligned(ac_don, ac_recni, ac_rec)
    freq_don = ac_don.to_frequencies()[:, 1]
    freq_recni = ac_recni.to_frequencies()[:, 1]
    freq_rec = ac_rec.to_frequencies()[:, 1]
    donor_100 = (freq_don == 1)
    non_intro_1 = (freq_recni < 0.01) 
    non_intro_10 = (freq_recni < 0.1)
    donor_100_non_intro_1 = (donor_100 & non_intro_1)
    donor_100_non_intro_10 = (donor_100 & non_intro_10)
    freqs_donor_100_non_intro_1 = freq_rec[donor_100_non_intro_1]
    freqs_donor_100_non_intro_10 = freq_rec[donor_100_non_intro_10]
    posi_donor_100_non_intro_1=np.array(pos)[donor_100_non_intro_1]
    posi_donor_100_non_intro_10=np.array(pos)[donor_100_non_intro_10]
    if freqs_donor_100_non_intro_10.size > 0:
        Q_10_100_q95, windows_Q_10_100_q95, counts_Q_10_100_q95 = windowed_statistic_percentil(posi_donor_100_non_intro_10, freqs_donor_100_non_intro_10, 95, size=size, start=start, stop=stop, step=step)
        Q_10_100_q90, windows_Q_10_100_q90, counts_Q_10_100_q90 = windowed_statistic_percentil(posi_donor_100_non_intro_10, freqs_donor_100_non_intro_10, 90, size=size, start=start, stop=stop, step=step)
        Q_10_100_max, windows_Q_10_100_max, counts_Q_10_100_max = allel.windowed_statistic(posi_donor_100_non_intro_10, freqs_donor_100_non_intro_10,statistic=np.max, size=size, start=start, stop=stop, step=step)
    else:
        count_Q_10, windows_Q_10 = allel.windowed_count(pos, size = size, start = start, stop = stop , step = step)
        Q_10_100_q95 = np.empty(len(windows_Q_10))
        Q_10_100_q95.fill(float('nan'))
        Q_10_100_q90 = np.empty(len(windows_Q_10))
        Q_10_100_q90.fill(float('nan'))
        Q_10_100_max = np.empty(len(windows_Q_10))
        Q_10_100_max.fill(float('nan'))
    if freqs_donor_100_non_intro_1.size > 0:
        Q_1_100_q95, windows_Q_1_100_q95, counts_Q_1_100_q95 = windowed_statistic_percentil(posi_donor_100_non_intro_1, freqs_donor_100_non_intro_1, 95, size=size, start=start, stop=stop, step=step)
        Q_1_100_q90, windows_Q_1_100_q90, counts_Q_1_100_q90 = windowed_statistic_percentil(posi_donor_100_non_intro_1, freqs_donor_100_non_intro_1, 90, size=size, start=start, stop=stop, step=step)
        Q_1_100_max, windows_Q_1_100_max, counts_Q_1_100_max = allel.windowed_statistic(posi_donor_100_non_intro_1, freqs_donor_100_non_intro_1,statistic=np.max, size=size, start=start, stop=stop, step=step)
    else:
        count_Q_1, windows_Q_1 = allel.windowed_count(pos, size=size, start= start, stop = stop , step=step)
        Q_1_100_q95 = np.empty(len(windows_Q_1))
        Q_1_100_q95.fill(float('nan'))
        Q_1_100_q90 = np.empty(len(windows_Q_1))
        Q_1_100_q90.fill(float('nan'))
        Q_1_100_max = np.empty(len(windows_Q_1))
        Q_1_100_max.fill(float('nan'))
    return Q_10_100_q95, Q_10_100_q90, Q_10_100_max, Q_1_100_q95, Q_1_100_q90, Q_1_100_max

###Use percentil() with windows derived from scikit-allel (allel.windowed_statistic()) (used in windowed_q_racimo)
def windowed_statistic_percentil(pos, values, arg, size=None, start=None, stop=None, step=None, windows=None, fill=np.nan):
    # assume sorted positions
    if not isinstance(pos, allel.SortedIndex):
        pos = allel.SortedIndex(pos, copy=False)
    # check lengths are equal
    if isinstance(values, tuple):
        # assume multiple values arrays
        allel.util.check_equal_length(pos, *values)
    else:
        # assume a single values array
        allel.util.check_equal_length(pos, values)
    # setup windows
    if windows is None:
        windows = allel.position_windows(pos, size, start, stop, step)
    else:
        windows = asarray_ndim(windows, 2)
    # find window locations
    locs = allel.window_locations(pos, windows)
    # setup outputs
    out = []
    counts = []
    # iterate over windows
    for start_idx, stop_idx in locs:
        # calculate number of values in window
        n = stop_idx - start_idx
        if n == 0:
            # window is empty
            s = fill
        else:
            if isinstance(values, tuple):
                # assume multiple values arrays
                wv = [v[start_idx:stop_idx] for v in values]
                s = np.percentile(*wv, arg)
            else:
                # assume a single values array
                wv = values[start_idx:stop_idx]
                s = np.percentile(wv, arg)
        # store outputs
        out.append(s)
        counts.append(n)
    # convert to arrays for output
    return np.asarray(out), windows, np.asarray(counts)

##Racimo's U (Racimo et al., 2017), adapted from Racimo et al., 2017 and Zhang et al., 2023 (MaLAdapt)
#Windowed Racimo's U need variant's position, scikit allel count from donor (popd), non-introgressed (poprecni) and recipient (poprec), size is the side of each window, start and stop the start and the end of the interest chromosome portion and step is the step between windows for slidding window. 
def windowed_u_racimo(pos, ac_popd, ac_poprecni, ac_poprec, size, start=None, stop=None, step=None):
    #allel count for each pop
    ac_don = allel.AlleleCountsArray(ac_popd, copy=False)
    if ac_don.shape[1] == 1 :
        nbre_variant = ac_don.shape[0]
        ac_don = np.append(ac_don, np.array([[0]]*nbre_variant),axis=1)
        ac_don = allel.AlleleCountsArray(ac_don, copy=False)
    ac_recni = allel.AlleleCountsArray(ac_poprecni, copy=False)
    if ac_recni.shape[1] == 1 :
        nbre_variant = ac_recni.shape[0]
        ac_recni = np.append(ac_recni, np.array([[0]]*nbre_variant),axis=1)
        ac_recni = allel.AlleleCountsArray(ac_recni, copy=False)
    ac_rec = allel.AlleleCountsArray(ac_poprec, copy=False)
    if ac_rec.shape[1] == 1 :
        nbre_variant = ac_rec.shape[0]
        ac_rec = np.append(ac_rec, np.array([[0]]*nbre_variant),axis=1)
        ac_rec = allel.AlleleCountsArray(ac_rec, copy=False)
    allel.util.check_dim0_aligned(ac_don, ac_recni, ac_rec)
    #derived allel frequency for each pop
    freq_don = ac_don.to_frequencies()[:, 1]
    freq_recni = ac_recni.to_frequencies()[:, 1]
    freq_rec = ac_rec.to_frequencies()[:, 1]
    # freq in donor pop == 100%
    donor_100 = (freq_don == 1)
    # freq in sister recipient pop (non_intro) < 1%
    non_intro_1 = (freq_recni < 0.01) 
    # freq in sister recipient pop (non_intro) < 10%
    non_intro_10 = (freq_recni < 0.1)
    # site with freq of derived allele ==100% in donor pop and < 1% in sister pop (= TRUE)
    donor_100_non_intro_1 = (donor_100 & non_intro_1)
    # site with freq of derived allele ==100% in donor pop and < 10% in sister pop (= TRUE)
    donor_100_non_intro_10 = (donor_100 & non_intro_10)
    # derived allel freq in recipient pop where freq donor == 100% and sister pop < 1%
    freqs_donor_100_non_intro_1 = freq_rec[donor_100_non_intro_1]
    # derived allel freq in recipient pop where freq donor == 100% and sister pop < 1%
    freqs_donor_100_non_intro_10 = freq_rec[donor_100_non_intro_10]
    # site with allel freq in donor pop == 100%, sister pop < 10 % and rec > 0 (= TRUE)
    U_10_0_100 = ( donor_100_non_intro_10 & (freq_rec > 0) )
    U_10_20_100 = ( donor_100_non_intro_10 & (freq_rec > 0.2) )
    U_10_50_100 = ( donor_100_non_intro_10 & (freq_rec > 0.5) )
    U_10_80_100 = ( donor_100_non_intro_10 & (freq_rec > 0.8) )
    U_1_0_100 = ( donor_100_non_intro_1 & (freq_rec > 0) )
    U_1_20_100 = ( donor_100_non_intro_1 & (freq_rec > 0.2) )
    U_1_50_100 = ( donor_100_non_intro_1 & (freq_rec > 0.5) )
    U_1_80_100 = ( donor_100_non_intro_1 & (freq_rec > 0.8) )
    #sum for each windows : Racimo's U 
    U_10_0_100, windows_U_10_0_100, counts_U_10_0_100 = allel.windowed_statistic(pos, U_10_0_100, statistic=np.nansum, size=size, start=start, stop=stop, step=step)
    U_10_20_100, windows_U_10_20_100, counts_U_10_20_100 = allel.windowed_statistic(pos, U_10_20_100, statistic=np.nansum, size=size, start=start, stop=stop, step=step)
    U_10_50_100, windows_U_10_50_100, counts_U_10_50_100 = allel.windowed_statistic(pos, U_10_50_100, statistic=np.nansum, size=size, start=start, stop=stop, step=step)
    U_10_80_100, windows_U_10_80_100, counts_U_10_80_100 = allel.windowed_statistic(pos, U_10_80_100, statistic=np.nansum, size=size, start=start, stop=stop, step=step)
    U_1_0_100, windows_U_1_0_100, counts_U_1_0_100 = allel.windowed_statistic(pos, U_1_0_100, statistic=np.nansum, size=size, start=start, stop=stop, step=step)
    U_1_20_100, windows_U_1_20_100, counts_U_1_20_100 = allel.windowed_statistic(pos, U_1_20_100, statistic=np.nansum, size=size, start=start, stop=stop, step=step)
    U_1_50_100, windows_U_1_50_100, counts_U_1_50_100 = allel.windowed_statistic(pos, U_1_50_100, statistic=np.nansum, size=size, start=start, stop=stop, step=step)
    U_1_80_100, windows_U_1_80_100, counts_U_1_20_100 = allel.windowed_statistic(pos, U_1_80_100, statistic=np.nansum, size=size, start=start, stop=stop, step=step)
    return U_10_0_100, U_10_20_100, U_10_50_100, U_10_80_100, U_1_0_100, U_1_20_100, U_1_50_100, U_1_80_100

##Windowed fD (Martin et al. 2015) 
#Adapted from Zhang et al., 2023 (MaLAdapt): 
def windowed_f_d (pos, ac_poprecni, ac_poprec, ac_popd, size, start=None, stop=None, step=None):
    #Only biallelic variant control :
    ac_recni = allel.AlleleCountsArray(ac_poprecni, copy=False)
    if ac_recni.shape[1] == 1 :
        nbre_variant = ac_recni.shape[0]
        ac_recni = np.append(ac_recni, np.array([[0]]*nbre_variant),axis=1)
        ac_recni = allel.AlleleCountsArray(ac_recni, copy=False)
    ac_rec = allel.AlleleCountsArray(ac_poprec, copy=False)
    if ac_rec.shape[1] == 1 :
        nbre_variant = ac_rec.shape[0]
        ac_rec = np.append(ac_rec, np.array([[0]]*nbre_variant),axis=1)
        ac_rec = allel.AlleleCountsArray(ac_rec, copy=False)
    ac_don = allel.AlleleCountsArray(ac_popd, copy=False)
    if ac_don.shape[1] == 1 :
        nbre_variant = ac_don.shape[0]
        ac_don = np.append(ac_don, np.array([[0]]*nbre_variant),axis=1)
        ac_don = allel.AlleleCountsArray(ac_don, copy=False)
    allel.util.check_dim0_aligned(ac_don, ac_recni, ac_rec)
    #derived allel frequency : 
    freq_recni = ac_recni.to_frequencies()[:, 1]#p1
    freq_rec = ac_rec.to_frequencies()[:, 1]#p2
    freq_don = ac_don.to_frequencies()[:, 1]#p3
    #ABBA and BABA calcul
    abba_vec = (1.0 - freq_recni)*freq_rec*freq_don
    baba_vec = freq_recni*(1.0 - freq_rec)*freq_don
    # if p2 > p3 : 
    check_fd1 = (freq_rec > freq_don)
    #ABBA = (1-p1)*p2*p2
    abba_fd1 = (1.0 - freq_recni)*freq_rec*freq_rec
    #BABA = p1*(1-p2)*p2
    baba_fd1 = freq_recni*(1.0 - freq_rec)*freq_rec
    #if p3 > p2 :
    check_fd2 = (freq_rec < freq_don)
    #ABBA = (1-p1)*p3*p3
    abba_fd2 = (1.0 - freq_recni)*freq_don*freq_don
    #BABA = p1*(1-p3)*p3
    baba_fd2 = freq_recni*(1.0 - freq_don)*freq_don
    abba_fd = check_fd1 * abba_fd1 + check_fd2 * abba_fd2
    baba_fd = check_fd1 * baba_fd1 + check_fd2 * baba_fd2
    abba_sum, windows_abba, counts_abba = allel.windowed_statistic(pos, abba_vec, statistic=np.nansum, size=size, start=start, stop=stop, step=step)
    baba_sum, windows_baba, counts_baba = allel.windowed_statistic(pos, baba_vec, statistic=np.nansum, size=size, start=start, stop=stop, step=step)
    abba_fd_sum, windows_abba_fd, counts_abba_fd = allel.windowed_statistic(pos, abba_fd, statistic=np.nansum, size=size, start=start, stop=stop, step=step)
    baba_fd_sum, windows_baba_fd, counts_baba_fd = allel.windowed_statistic(pos, baba_fd, statistic=np.nansum, size=size, start=start, stop=stop, step=step)
    f_D = (abba_sum - baba_sum) / (abba_fd_sum - baba_fd_sum)
    return f_D

####HAPLOTYPE HOMOZYGOSITY###############
##Haplotype homozygosity (Garud et al. (2015) , H1, H12, H123, H2_H1)
# Taken from scikit-allel (garud_h()) 
def modify_garud_h(h):
    # check inputs
    h = allel.HaplotypeArray(h, copy=False)
    # compute haplotype frequencies
    f = h.distinct_frequencies()
    # compute H1
    h1 = np.sum(f**2)
    # compute H12
    h12 = np.sum(f[:2])**2 + np.sum(f[2:]**2)
    # compute H123
    h123 = np.sum(f[:3])**2 + np.sum(f[3:]**2)
    # compute H2/H1
    h2 = h1 - f[0]**2
    h2_h1 = h2 / h1
    return h1, h2, h12, h123, h2_h1

# Adapted from scikit-allel (moving_garud_h()):
def windowed_garud_h(pos, h, size, start=None, stop=None, step=None):
    gh, windows_gh, counts_gh = allel.windowed_statistic(pos, values=h, statistic=modify_garud_h, size=size, start=start, stop=stop, step=step)
    if len(gh.shape) == 1 :
        newgh  = np.where([type(x)==float for x in gh],"00000",gh)
        newgh  = np.array([list(x) for x in newgh])
        h1     = np.array([float(x) for x in newgh[:,0]])
        h2     = np.array([float(x) for x in newgh[:,1]])
        h12    = np.array([float(x) for x in newgh[:,2]])
        h123   = np.array([float(x) for x in newgh[:,3]])
        h2_h1  = np.array([float(x) for x in newgh[:,4]])
    else :
        h1     = gh[:, 0]
        h2     = gh[:, 1]
        h12    = gh[:, 2]
        h123   = gh[:, 3]
        h2_h1  = gh[:, 4]
    return h1, h2, h12, h123, h2_h1



#Number of segregating sites by window :
def windowed_count_segregating(pos, ac, size, start=None, stop=None, step=None):
    is_segreg = ac.is_segregating()
    s_site, windows_s_site, counts_s_site = allel.windowed_statistic(pos, values=is_segreg, statistic=np.sum, size=size, start=start, stop=stop, step=step)
    return s_site

####LINKAGE DESEQUILIBRIUM SUMMARY STATISTICS###############
##r2 (Hill et Robertson, 1968) :
#Taken from Zhang et al., 2023 script:
def r2(pAB,paB,pAb,pab):
    pA  = (pAB+pAb)
    pB  = (pAB+paB)
    pa  = (paB+pab)
    pb  = (pAb+pab)
    D   = pAB-pA*pB
    Dab = pab-pa*pb
    if (pA*(1-pA)*pB*(1-pB) != 0 ):
        r2 = (D**2)/(pA*(1-pA)*pB*(1-pB))
    else:
        r2 = (Dab**2)/(pa*(1-pa)*pb*(1-pb))
    return r2

##ZnS (Kelly, 1997)
#Adapt ZnS from Zhang et al., 2023 fuction to be calculated on haplotype matrix from scikit-allel:
def Zns_good(geno_array_scikit):
    #genotype array
    genotype_array_pop = geno_array_scikit
    #allel count
    allel_count_pop = genotype_array_pop.count_alleles()
    #derived allel frequency
    freq_pop_anc = allel_count_pop.to_frequencies()[:, 0]
    freq_pop_deriv = allel_count_pop.to_frequencies()[:, 1]
    #transpose genotype array
    t_geno = np.transpose(np.array(genotype_array_pop))
    #select only segegating sites 
    S_geno_1 = np.asarray([g[(freq_pop_deriv>0)&(freq_pop_deriv<1)] for g in t_geno[0]])
    S_geno_2 = np.asarray([g[(freq_pop_deriv>0)&(freq_pop_deriv<1)] for g in t_geno[1]])
    #Concatenate new gen array with only segregating sites 
    S_geno = np.array([S_geno_1, S_geno_2])
    #number of genotype
    S_gen = S_geno.shape[1]
    #number of variant
    S_var = S_geno.shape[2]
    #if segregating sites exist : 
    if S_var>1 :
        ls = range(1,S_var+1)
        # create list of all site's pairs possible (position)
        pairs = list(itertools.combinations(ls,2))
        liste_liste = []
        #create genotype array 2D (racimo and maladapt like)
        for i in range(0, S_gen) : 
            liste = list(zip(S_geno[0][i],S_geno[1][i]))
            liste_string = [''.join(map(str, x)) for x in liste]
            liste_liste.append(liste_string)
        #create all derived and ancestral combinaison based on pairs 
        genos_combos = [[list(zip(g[x-1],g[y-1])) for g in liste_liste] for x,y in pairs]
        geno = [[''.join(i) for n in g for i in n] for g in genos_combos]
        #calcul the number of AB (11)combinaison 
        AB = [np.sum([y == '11' for y in x]) for x in geno]
        #calcul the number of aB (01)combinaison 
        aB = [np.sum([y == '01' for y in x]) for x in geno]
        #calcul the number of Ab (10)combinaison
        Ab = [np.sum([y == '10' for y in x]) for x in geno]
        #calcul the number of ab (00)combinaison
        ab = [np.sum([y == '00' for y in x]) for x in geno]
        ps = [[w/float((w+x+y+z)),x/float((w+x+y+z)),y/float((w+x+y+z)),z/float((w+x+y+z))] for w,x,y,z in list(zip(AB,aB,Ab,ab))]
        r2s = [r2(w,x,y,z) for w,x,y,z in ps]
        #calcul zns
        Zns = np.sum(r2s)*2./(S_var*(S_var-1))
        return Zns


##Windowed ZnS:
def windowed_zns(position_array, geno_array_scikit, w_size, w_start = None , w_stop = None , w_step = None) :
    geno_array = geno_array_scikit
    count_zns, windo_zns = allel.windowed_count(position_array, size = w_size, start= w_start, stop = w_stop , step = w_step)
    zns_list = []
    #for each window :
    for w in range(len(windo_zns)):
        #windows :
        #variants present in the window :
        thesepos = [pos for pos in position_array if (pos <= windo_zns[w][1] and pos >= windo_zns[w][0])]
        #variant index 
        if len(thesepos) > 0 :
            start = position_array.index(thesepos[0])
            end = position_array.index(thesepos[len(thesepos)-1])+1
            #window's lenght 
            len_seg = windo_zns[w][1]-windo_zns[w][0]
            #windows haplotypes by pop 
            pop_geno = geno_array[start:end]
            list_geno_freq = [len(x) for x in pop_geno.count_alleles().to_frequencies()]
            #and np.all(np.array(list_geno_freq) == 1)
            if len(list_geno_freq) == 1 or np.all(np.array(list_geno_freq) == 1) : 
                zns_final = float("NaN")
                zns_list.append(zns_final)
            else :
                zns_final = Zns_good(pop_geno)
                zns_list.append(zns_final)
        else :
            zns_final = float("NaN")
            zns_list.append(zns_final)
    return np.array(zns_list)

#Moment metrics and other summary statistics from window summary statistics dataframe  :
def calcul_stat_ref_table (df_summary_stat, sim_index) :
    data_without_interv = df_summary_stat.drop(['start', 'end'], axis=1)
    sum_stat_names      = np.array(data_without_interv.columns)
    stat_names          = np.array(["mean","variance","skewness","kurtosis","min","max","quant_5","quant_95"])
    all_col_names       = np.array([''.join([x +"_"+ y])for x in stat_names for y in sum_stat_names])
    mean_array          = np.array(np.mean(data_without_interv, axis=0))
    var_array           = np.array(np.var(data_without_interv, axis=0))
    skewness_array      = np.array(skew(data_without_interv, axis=0, bias=True,  nan_policy='omit'))
    kurtosis_array      = np.array(kurtosis(data_without_interv, axis=0, bias=True,  nan_policy='omit'))
    min_array           = np.array(np.amin(data_without_interv, axis=0))
    max_array           = np.array(np.amax(data_without_interv, axis=0))
    quant_5_array       = np.array(data_without_interv.quantile(0.05))
    quant_95_array      = np.array(data_without_interv.quantile(0.95))
    all_stat_array      = np.hstack((mean_array, var_array, skewness_array, kurtosis_array, min_array, max_array, quant_5_array, quant_95_array))
    ref_table           = pd.DataFrame(data=all_stat_array)
    ref_table           = ref_table.transpose()
    ref_table.columns   = all_col_names
    ref_table.insert(0, "num_sim", sim_index)
    ref_table           = ref_table.fillna(value="NaN")
    ref_table           = ref_table.replace("NaN", 0.0)
    return ref_table
