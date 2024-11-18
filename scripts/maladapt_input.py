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

import pandas as pd
from introadapt import *

#Script to rename pipeline column names in MaLAdapt column name, keep MaLAdapt columns for dataframe with and without Kelly's ZnS :

options     = get_project_options(proj_options_file = sys.argv[1])

project_dir = options["results_dir"]+"/"+options["analysis"]+"/"+options["project"]

ref_stat    = pd.read_csv(project_dir+"/reftable_wind_"+options["analysis"]+"_"+options["project"]+".csv", sep=" ")
if options["zns"]=="Yes" :
	column_convertor = {"IA":"classifier",	"num_sim":"sim_id",	"start":"start",	"end":"end",	"mean_p1":"mean_p1",	"fixation_g_r_back":"fixation_g_r_back",
	"green_d":"Dstat",	"fd":"fD",	"exp_het_rec":"Het",	"rd":"divratioavg",	"Q_10_100_q95":"Q_10_100_q95",	"Q_10_100_q90":"Q_10_100_q90",	
	"Q_10_100_max":"Q_10_100_max",	"Q_1_100_q95":"Q_1_100_q95",	"Q_1_100_q90":"Q_1_100_q90",	"Q_1_100_max":"Q_1_100_max",	"U_10_0_100":"U_10_0_100",
	"U_10_20_100":"U_10_20_100",	"U_10_50_100":"U_10_50_100",	"U_10_80_100":"U_10_80_100",	"U_1_0_100":"U_1_0_100",	"U_1_20_100":"U_1_20_100",	
	"U_1_50_100":"U_1_50_100",	"U_1_80_100":"U_1_80_100",	"S_don":"p1_S",	"S_rec_sister":"p2_S",	"S_rec":"p3_S",	"w_theta_don":"thetaW1",	
	"w_theta_rec_sister":"thetaW2",	"w_theta_rec":"thetaW3",	"pi_don":"thetapi1",	"pi_rec_sister":"thetapi2",	"pi_rec":"thetapi3",	"h_theta_don":"thetaH1",	
	"h_theta_rec_sister":"thetaH2",	"h_theta_rec":"thetaH3",	"h1_don":"p1_H1",	"h2_don":"p1_H2",	"h12_don":"p1_H12",	"h2_h1_don":"p1_H2H1",	
	"h1_rec_sister":"p2_H1",	"h2_rec_sister":"p2_H2",	"h12_rec_sister":"p2_H12",	"h2_h1_rec_sister":"p2_H2H1",	"h1_rec":"p3_H1",	
	"h2_rec":"p3_H2",	"h12_rec":"p3_H12",	"h2_h1_rec":"p3_H2H1",	"zns_do":"Zns1",	"zns_rec_sister":"Zns2",	"zns_rec":"Zns3",	"r":"r"}
else :
	column_convertor = {"IA":"classifier",	"num_sim":"sim_id",	"start":"start",	"end":"end",	"mean_p1":"mean_p1",	"fixation_g_r_back":"fixation_g_r_back",
	"green_d":"Dstat",	"fd":"fD",	"exp_het_rec":"Het",	"rd":"divratioavg",	"Q_10_100_q95":"Q_10_100_q95",	"Q_10_100_q90":"Q_10_100_q90",	
	"Q_10_100_max":"Q_10_100_max",	"Q_1_100_q95":"Q_1_100_q95",	"Q_1_100_q90":"Q_1_100_q90",	"Q_1_100_max":"Q_1_100_max",	"U_10_0_100":"U_10_0_100",
	"U_10_20_100":"U_10_20_100",	"U_10_50_100":"U_10_50_100",	"U_10_80_100":"U_10_80_100",	"U_1_0_100":"U_1_0_100",	"U_1_20_100":"U_1_20_100",	
	"U_1_50_100":"U_1_50_100",	"U_1_80_100":"U_1_80_100",	"S_don":"p1_S",	"S_rec_sister":"p2_S",	"S_rec":"p3_S",	"w_theta_don":"thetaW1",	
	"w_theta_rec_sister":"thetaW2",	"w_theta_rec":"thetaW3",	"pi_don":"thetapi1",	"pi_rec_sister":"thetapi2",	"pi_rec":"thetapi3",	"h_theta_don":"thetaH1",	
	"h_theta_rec_sister":"thetaH2",	"h_theta_rec":"thetaH3",	"h1_don":"p1_H1",	"h2_don":"p1_H2",	"h12_don":"p1_H12",	"h2_h1_don":"p1_H2H1",	
	"h1_rec_sister":"p2_H1",	"h2_rec_sister":"p2_H2",	"h12_rec_sister":"p2_H12",	"h2_h1_rec_sister":"p2_H2H1",	"h1_rec":"p3_H1",	
	"h2_rec":"p3_H2",	"h12_rec":"p3_H12",	"h2_h1_rec":"p3_H2H1",	"r":"r"}
ref_stat = ref_stat.rename(columns=column_convertor)#renames columns to match maladapt
ref_stat = ref_stat[list(column_convertor.values())]#keeps only the column named in the dictionnary, rearranges their order (summary statistics order is based on maladapt reftable, extra latent variables and parameters could be kept safely)
ref_stat.to_csv(project_dir+"/maladapt_tab_"+options["project"]+".csv", mode='w', header=True, index=False)
print("Check_mal_done")