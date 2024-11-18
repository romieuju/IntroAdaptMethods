#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#    IntroAdapt: joint inference of introgression and selection
#    Copyright (C) 2024 Ghislain Camarata and Jules Romieu, ISEM/CBGP, CNRS/INRAE. 
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

import allel
import numpy as np

#creates two tables that will be turned into txt files used as input for VolcanoFinder (SFS and AFF)

#Version with SFS for each chromosome :
def volcano_tables(recipient_indices, outgroup_indices, var_pos_tot, chr_length, haplo_all) : #add haplo_chro and position 
    volc_samp             = outgroup_indices+recipient_indices
    haplo_volc_allsites   = haplo_all[:,volc_samp]#keeps only the outgroup and recipient populations haplotypes, which are the only ones used by this method
    haplo_popO_allsites   = haplo_all[:,outgroup_indices]#keeps only the outgroup
    allel_count_out_full  = haplo_popO_allsites.count_alleles()#counts in the outgroup population (one row per SNP)
    volc_pos_out          = allel_count_out_full.is_non_variant()#check that the outgroup has an ancestral allele
    allel_count_volc_full = haplo_volc_allsites.count_alleles()#counts in the outgroup and recipient populations (one row per SNP)
    volc_pos_bia          = allel_count_volc_full.is_biallelic()#check that there are 2 alleles in this sample
    volc_pos_keep         = np.logical_and(volc_pos_out,volc_pos_bia)#combines the two filters
    haplo_rec_volc        = haplo_all.subset(sel0=volc_pos_keep,sel1=recipient_indices)#keeps only the SNP that are relevant for this method for the recipient population
    col_positions_volc    = np.array(var_pos_tot)[np.array(volc_pos_keep)]#keeps the relevant positions
    allel_count_volc      = haplo_rec_volc.count_alleles()#counts once more, but anly for the informative loci in the recipient population
    col_x_volc            = np.sum(allel_count_volc[:,1:],axis=1)#counts the number of variants for each SNP
    col_n_volc            = np.sum(allel_count_volc,axis=1)#counts the number of samples for each SNP
    col_folded_volc       = np.array([0]*len(col_positions_volc))
    AFF_columns           = np.column_stack((col_positions_volc,col_x_volc,col_n_volc,col_folded_volc))
    SFS_rec               = allel.sfs(col_x_volc)#creates the SFS for the recipient population
    SFS_rec               = np.delete(SFS_rec,0)#removes the monomorphic ancestral sites (neither used nor calculated)
    col_SFS_volc          = SFS_rec/chr_length#turns our SFS into proportions of callable sites
    con_nrow_volc         = np.array(range(len(recipient_indices)))+1
    SFS_columns           = np.column_stack((con_nrow_volc,col_SFS_volc))
    return AFF_columns, SFS_columns

#Version with SFS from all the genome :
def volcano_tables_all(recipient_indices, outgroup_indices, var_pos_tot, haplo_chr, ge_length, haplo_all) : #add haplo_chro and position 
    volc_samp             = outgroup_indices+recipient_indices
    haplo_volc_chrosites  = haplo_chr[:,volc_samp]#keeps only the outgroup and recipient populations haplotypes, which are the only ones used by this method
    haplo_popO_chrosites  = haplo_chr[:,outgroup_indices]#keeps only the outgroup
    allel_count_out_full  = haplo_popO_chrosites.count_alleles()#counts in the outgroup population (one row per SNP)
    volc_pos_out          = allel_count_out_full.is_non_variant()#check that the outgroup has an ancestral allele
    allel_count_volc_full = haplo_volc_chrosites.count_alleles()#counts in the outgroup and recipient populations (one row per SNP)
    volc_pos_bia          = allel_count_volc_full.is_biallelic()#check that there are 2 alleles in this sample
    volc_pos_keep         = np.logical_and(volc_pos_out,volc_pos_bia)#combines the two filters
    haplo_rec_volc        = haplo_chr.subset(sel0=volc_pos_keep,sel1=recipient_indices)#keeps only the SNP that are relevant for this method for the recipient population
    col_positions_volc    = np.array(var_pos_tot)[np.array(volc_pos_keep)]#keeps the relevant positions
    allel_count_volc      = haplo_rec_volc.count_alleles()#counts once more, but anly for the informative loci in the recipient population
    col_x_volc            = np.sum(allel_count_volc[:,1:],axis=1)#counts the number of variants for each SNP
    col_n_volc            = np.sum(allel_count_volc,axis=1)#counts the number of samples for each SNP
    col_folded_volc       = np.array([0]*len(col_positions_volc))
    AFF_columns           = np.column_stack((col_positions_volc,col_x_volc,col_n_volc,col_folded_volc))
    #SFS all genome :
    haplo_volc_allsites       = haplo_all[:,volc_samp]#keeps only the outgroup and recipient populations haplotypes, which are the only ones used by this method
    haplo_popO_allsites       = haplo_all[:,outgroup_indices]#keeps only the outgroup
    allel_count_out_full_all  = haplo_popO_allsites.count_alleles()#counts in the outgroup population (one row per SNP)
    volc_pos_out_all          = allel_count_out_full_all.is_non_variant()#check that the outgroup has an ancestral allele
    allel_count_volc_full_all = haplo_volc_allsites.count_alleles()#counts in the outgroup and recipient populations (one row per SNP)
    volc_pos_bia_all          = allel_count_volc_full_all.is_biallelic()#check that there are 2 alleles in this sample
    volc_pos_keep_all         = np.logical_and(volc_pos_out_all,volc_pos_bia_all)#combines the two filters
    haplo_rec_volc_all        = haplo_all.subset(sel0=volc_pos_keep_all,sel1=recipient_indices)#keeps only the SNP that are relevant for this method for the recipient population
    allel_count_volc_all      = haplo_rec_volc_all.count_alleles()#counts once more, but anly for the informative loci in the recipient population
    col_x_volc_all            = np.sum(allel_count_volc_all[:,1:],axis=1)#counts the number of variants for each SNP
    SFS_rec_all               = allel.sfs(col_x_volc_all)#creates the SFS for the recipient population
    SFS_rec_all               = np.delete(SFS_rec_all,0)#removes the monomorphic ancestral sites (neither used nor calculated)
    col_SFS_volc_all          = SFS_rec_all/ge_length#turns our SFS into proportions of callable sites
    con_nrow_volc_all         = np.array(range(len(recipient_indices)))+1
    SFS_columns_all           = np.column_stack((con_nrow_volc_all,col_SFS_volc_all))
    return AFF_columns, SFS_columns_all