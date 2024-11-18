#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#    IntroAdapt: joint inference of introgression and selection
#    Copyright (C) 2024  Ghislain Camarata and Jules Romieu, ISEM/CBGP, CNRS/INRAE. 
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

#Script to obtain appropriate files for existing AI inference methods (VolcanoFinder and genomatnn) from tree seq after addition of neutral mutations with msprime (allel package)

#Python package used :
import pyslim, os, allel, sys, csv, tskit
import numpy as np
import pandas as pd
from introadapt import *
from sum_stat_function import *
from volcanofinder_input import *


#Obtain information from ini files (timeadapt.py script)
options = get_options(proj_options_file = sys.argv[1], sim_options_file = sys.argv[2])

#Save results folder path in a variable (results_dir/analyse/project/project_sim)
project_dir = options["results_dir"]+"/"+options["analysis"]+"/"+options["project"]
job_dir     = project_dir+"/"+options["analysis"]+"_"+options["project"]+"_sim_"+options["sim"]


#load parameters from ini files 
#input and output files :
treeseq          = tskit.load(job_dir+"/Mut_TreeSeq_"+options["analysis"]+"_"+options["project"]+"_sim_"+options["sim"]+".trees")#treeseq = tskit.load("Mut_TreeSeq_1.trees")
#check which methods are suppoded to be used
volcano_switch   = options["volcano_switch"]
genomatnn_switch = options["genomatnn_switch"]

#Save samples size for each population 
n_popA1 = options["n_samples"][0]
n_popA2 = options["n_samples"][1]
n_popB1 = options["n_samples"][2]
n_popB2 = options["n_samples"][3]
n_popO  = options["n_samples"][4]

#genome length 
ge_l = options['L']
#Save chromosome data 
nchr        = options['nchr']
chr_ends    = options['chr_ends']
chr_starts  = [0]+chr_ends[:-1]
chr_lengths = np.array(chr_ends)-np.array(chr_starts)

#Save samples size
s_size = options["sample_size"]

#Save end time of forward (slim) generation
g_end_sim        = options["generations_forward"]
#Save sampling time of the donor pop
sample_time_d    = options["sample_time_d"]
#Backward generation between end of sim and don pop sampling :
g_samp_don       = g_end_sim - sample_time_d
#idem for donor sister pop :
sample_time_do_s = options["sample_time_d_s"]
g_samp_don_s     = g_end_sim - sample_time_do_s

#Variable used if user want add error and coverage informations in his samples sequences :
#If the genotype are phased or not :
is_phased = options['is_phased']

#which sample should represent the donor
donor_pop_stat = options["donor_pop_stat"]

################################################################################

#Save variants positions in list 
positions =[variant.position for variant in treeseq.variants()]

# Export the genotype data to allel. Unfortunately there's a slight mismatch in the 
# terminology here where genotypes and haplotypes mean different things in the two
# libraries.
tree_genotype_matrix = treeseq.genotype_matrix()#This is a haplotype matrix (each row is a SNP and each column is an individual)
#tree seq genotype matrix to scikit haplotype matrix :  
haplo_all = allel.HaplotypeArray(tree_genotype_matrix)#we just changed the format

if len(tree_genotype_matrix) == 0 :
    raise ValueError('Genotype matrix is empty !')

#save the id of the nodes of each population in a specific variable
nodes_outgroup = [x.id for x in treeseq.nodes() if ((x.time == 1.0) and (x.population == 1))]
nodes_popA1    = [x.id for x in treeseq.nodes() if ((x.time == g_samp_don_s) and (x.population == 5))]
nodes_popA2    = [x.id for x in treeseq.nodes() if ((x.time == g_samp_don) and (x.population == 6))]
nodes_popB1    = [x.id for x in treeseq.nodes() if ((x.time == 1.0) and (x.population == 7))]
nodes_popB2    = [x.id for x in treeseq.nodes() if ((x.time == 1.0) and (x.population == 8))]

#define donor pop and sister donor pop :
if options["donor"] == '"A1"' :
    nodes_donor = nodes_popA1
    nodes_don_sis = nodes_popA2
    n_don=n_popA1
    n_don_sis=n_popA2
elif options["donor"] == '"A2"' :
    nodes_donor = nodes_popA2
    nodes_don_sis = nodes_popA1
    n_don=n_popA2
    n_don_sis=n_popA1
elif options["donor"] == '"B1"' :
    nodes_donor = nodes_popB1
    nodes_don_sis = nodes_popB2
    n_don=n_popB1
    n_don_sis=n_popB2
elif options["donor"] == '"B2"' :
    nodes_donor = nodes_popB2
    nodes_don_sis = nodes_popB1
    n_don=n_popB2
    n_don_sis=n_popB1

#define recipient pop and sister recipient pop :
if options["recipient"] == '"A1"' :
    nodes_recipient = nodes_popA1
    n_rec=n_popA1
    nodes_rec_sis = nodes_popA2
    n_rec_sis=n_popA2
elif options["recipient"] == '"A2"' :
    nodes_recipient = nodes_popA2
    n_rec=n_popA2
    nodes_rec_sis = nodes_popA1
    n_rec_sis=n_popA1
elif options["recipient"] == '"B1"' :
    nodes_recipient = nodes_popB1
    n_rec=n_popB1
    nodes_rec_sis = nodes_popB2
    n_rec_sis=n_popB2
elif options["recipient"] == '"B2"' :
    nodes_recipient = nodes_popB2
    n_rec=n_popB2
    nodes_rec_sis = nodes_popB1
    n_rec_sis=n_popB1

#to use donor sister population for summary stat calculation :
if donor_pop_stat != options["donor"] :
    nodes_donor = nodes_don_sis
    n_don        = n_don_sis

#deduces the individuals' id from the nodes' 
ind_outgroup  = list(set([int(i/2) for i in nodes_outgroup]))
ind_don_sis   = list(set([int(i/2) for i in nodes_don_sis]))
ind_donor    = list(set([int(i/2) for i in nodes_donor]))
ind_rec_sis   = list(set([int(i/2) for i in nodes_rec_sis]))
ind_recipient = list(set([int(i/2) for i in nodes_recipient]))

##############################################################################################################################################

if volcano_switch=="On":
    vol_dir=job_dir+"/VolcanoFinder/"
    print("Check_vol_1")
    if not os.path.exists(vol_dir):
        os.mkdir(vol_dir)
        print("Check_vol_2")
    for chromosome in range(1,nchr+1):
        print("Check_vol_3"+"_chr_"+str(chromosome))

        #saves the caracteristics of the specific chromosome we're interested in
        chr_l = chr_lengths[chromosome-1]
        chr_s = chr_starts[chromosome-1]
        chr_e = chr_ends[chromosome-1]

        #isolates the positions that are part of this chromosome
        right_chr = np.logical_and(np.array(positions)>=chr_s,np.array(positions)<chr_e)
        chr_pos   = np.array(positions)[np.array(right_chr)]
        haplo_chr = haplo_all.subset(sel0=right_chr)

        ##For use SFS of each chromosome
        AFF_columns, SFS_columns = volcano_tables(nodes_recipient, nodes_outgroup, chr_pos, chr_l, haplo_chr)#defined in "fonctions_ghislain" #add position, haplo_all, ge_l
        #For use SFS of all the genome
        #AFF_columns, SFS_columns = volcano_tables_all(nodes_recipient, nodes_outgroup, chr_pos, haplo_chr, ge_l, haplo_all)

        #writes VolcanoFinder's input files
        with open(vol_dir+"AFF_volcano_"+options["sim"]+"_chr_"+str(chromosome)+".txt",'w') as f:#options = {"sim": "6"}
            np.savetxt(f, AFF_columns,delimiter='\t',header="position\tx\tn\tfolded",comments="",fmt='%i')
        with open(vol_dir+"SFS_volcano_"+options["sim"]+"_chr_"+str(chromosome)+".txt",'w') as f:
            np.savetxt(f, SFS_columns,delimiter='\t',fmt='%i %e')

if genomatnn_switch=="On":
    print("Check_gnn_condition")

    #saves gnn results path in a variable and creates the directory if it doesn't already exist
    gnn_dir = job_dir+"/genomatnn/"
    if not os.path.exists(gnn_dir):
        os.mkdir(gnn_dir)
        print("Check_gnn_dir")

    #saves a list of each relevant population's individals id in a specific file, so that genomatnn can actually read the vcf 
    vcf_id_donor    = ["tsk_"+str(i) for i in ind_donor]
    vcf_id_rec_sis   = ["tsk_"+str(i) for i in ind_rec_sis]
    vcf_id_recipient = ["tsk_"+str(i) for i in ind_recipient]
    with open(gnn_dir+"donor_"+options["sim"]+".indlist",'w') as f:#options = {"sim": "1"}
        f.write('\n'.join(str(i) for i in vcf_id_donor))
    with open(gnn_dir+"sister_"+options["sim"]+".indlist",'w') as f:#options = {"sim": "1"}
        f.write('\n'.join(str(i) for i in vcf_id_rec_sis))
    with open(gnn_dir+"recipient_"+options["sim"]+".indlist",'w') as f:#options = {"sim": "1"}
        f.write('\n'.join(str(i) for i in vcf_id_recipient))
    print("Check_gnn_indlists")

    #writes a first version of the vcf
    with open(gnn_dir+"don-rec-sis_"+options["sim"]+".vcf","w") as vcff :
        pyslim.convert_alleles(pyslim.generate_nucleotides(treeseq)).write_vcf(output=vcff,contig_id="UNKNOWN")#contig_id gives the chromosome's name, generate nucleotide and convert alleles make the vcf bcftools-compatible

    #breaks down the vcf into tiny little portions to add chromosome information
    data_header     = pd.read_csv(gnn_dir+"don-rec-sis_"+options["sim"]+".vcf", sep="\t",header=None,nrows=3)
    contig_template = pd.read_csv(gnn_dir+"don-rec-sis_"+options["sim"]+".vcf", sep="\t",header=None,skiprows=3,nrows=1).iloc[0,0]

    data_contig     = pd.DataFrame()#given as a replacement for the contig info in the vcf header
    chr_vcf_col     = []#given as a replacement for the first column of the vcf
    chr_enumeration = ""#given as an argument in the config file
    for chr_id in range(0,nchr) :
        chr_vcf_col += [str(chr_id+1)]*sum(np.logical_and(chr_starts[chr_id]<=np.array(positions),np.array(positions)<chr_ends[chr_id]))#sum([chr_starts[chr_id]<=snp<chr_ends[chr_id] for snp in positions])#chr_vcf_vol=[[chr_id+1]*sum(np.logical_and(chr_starts[chr_id]<=np.array(positions),np.array(positions)<chr_ends[chr_id])) for chr_id in range(0,nchr)]#sum([chr_starts[chr_id]<=snp<chr_ends[chr_id] for snp in positions])
        
        new_chr     = contig_template.replace("UNKNOWN",str(chr_id+1))
        new_chr     = new_chr.replace(str(chr_ends[-1]+1),str(chr_ends[chr_id]+1))
        new_chr     = pd.DataFrame([new_chr])
        data_contig = pd.concat([data_contig,new_chr])
        
        chr_enumeration += ","+str(chr_id+1)

    if chr_ends[-1] == positions[-1] :
        chr_vcf_col += str(nchr)

    data_format = pd.read_csv(gnn_dir+"don-rec-sis_"+options["sim"]+".vcf", sep="\t",header=None,skiprows=4,nrows=1)
    data_names  = pd.read_csv(gnn_dir+"don-rec-sis_"+options["sim"]+".vcf", sep="\t",header=None,skiprows=5,nrows=1)
    data_vcf    = pd.read_csv(gnn_dir+"don-rec-sis_"+options["sim"]+".vcf", sep="\t",header=None,skiprows=6)

    #replaces inapropriate column (resulting from tskit shortcomings)
    data_vcf[0] = chr_vcf_col

    #assembles the different component of the vcf and replaces the file
    vcf_full = pd.concat([data_header,data_contig,data_format,data_names,data_vcf])
    vcf_full.to_csv(gnn_dir+"don-rec-sis_"+options["sim"]+".vcf",header=False,index=None,sep="\t",mode='w',quoting=csv.QUOTE_NONE)

    #are the haplotypes phased ?
    phasing = "true" if is_phased == "Yes" else "false"

    #parse the configuration file using awk
    awk_command="gawk -f scripts/parse_gnn.awk -v path={} -v nosimu={} -v nochr={} -v phasing={} config_file/gnn_config.toml > {}"
    awk_command=awk_command.format(gnn_dir,options['sim'],chr_enumeration[1:],phasing,gnn_dir+'config_'+options['sim']+'.toml')
    os.system(awk_command)#awk_command='gawk -f scripts/parse_gnn.awk -v path=/home/camarata/Bureau/results/brouillon/awk_test/awk_test_1/genomatnn/ -v nosimu=1 -v nochr=1 tests/gnn_config.toml > /home/camarata/Bureau/results/brouillon/awk_test/awk_test_1/genomatnn/config_1.toml'

    print("Check_gnn_VCF")

