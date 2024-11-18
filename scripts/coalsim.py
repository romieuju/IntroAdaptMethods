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

#Backward simulation with msprime : This script simule genealogy of all population without neutral mutation (add in mutsim.py) and create a tree sequence file which will be read by forward simulator (SLiM)

#Python packages used in this script
import sys
import msprime
import pyslim
import introadapt
import time 


def principal():
	msprime_start = time.time()

	# get options for project and simulation:
	options = introadapt.get_options(proj_options_file = sys.argv[1], sim_options_file = sys.argv[2])

	
	
	# print program name
	introadapt.print_info(sys.argv[0],options["verbose"],sim=options["sim"])
	
	# Save the simulation parameters : type of coalescent split (split or mass_migration) and if some split between population are will be made in SliM (0 or 1)
	coalescence_split = str(options["coalescence_split"])
	slim_split        = options["slim_split"]
	

	#Save the results_dir path :
	project_dir = options["results_dir"]+"/"+options["analysis"]+"/"+options["project"]
	job_dir = project_dir+"/"+options["analysis"]+"_"+options["project"]+"_sim_"+options["sim"]

	#If the user want do all population divergence in msprime (slim_plit = 0)
	if (slim_split == 0) :
		
		# save recombination map:
		rate_map = msprime.RateMap(position = options["msprime_r_map"]["positions"],
				     rate = options["msprime_r_map"]["rates"]) 
		
		# get demography
		demography = msprime.Demography()
		#Initialisy each population with their population size (N) :
		#Outgroup + other pop ancestral pop size
		demography.add_population(name = "popAncO", initial_size = options["popAncO"]) #p0
		#Outgroup pop size
		demography.add_population(name = "popO", initial_size = options["popO"]) #p1
		#donor/recipient ancestral pop size
		demography.add_population(name = "popAncAB", initial_size = options["popAncAB"]) #p2
		#donor/donor sister ancestral pop size
		demography.add_population(name = "popAncA", initial_size = options["popAncA"]) #p3
		#Recipient/Recipient sister ancestral pop size
		demography.add_population(name = "popAncB", initial_size = options["popAncB"]) #p4
		#Sister donor population size
		demography.add_population(name = "popA1", initial_size = options["popA1"]) #p5
		#donor population size
		demography.add_population(name = "popA2", initial_size = options["popA2"]) #p6
		#Recipient sister pop size
		demography.add_population(name = "popB1", initial_size = options["popB1"]) #p7
		#Recipient pop size
		demography.add_population(name = "popB2", initial_size = options["popB2"]) #p8
		

		# Divergence will be made by mass_migration option in Msprime
		if (coalescence_split == "mass_migration") :
			demography.add_mass_migration(time=options["generation_split_B"],source = "popB2", dest = "popAncB", proportion = 1)
			demography.add_mass_migration(time=options["generation_split_B"],source = "popB1", dest = "popAncB", proportion = 1)
			demography.add_mass_migration(time=options["generation_split_A"],source = "popA2", dest = "popAncA", proportion = 1)
			demography.add_mass_migration(time=options["generation_split_A"],source = "popA1", dest = "popAncA", proportion = 1)
			demography.add_mass_migration(time=options["generation_split_AB"],source = "popAncA", dest = "popAncAB", proportion = 1)
			demography.add_mass_migration(time=options["generation_split_AB"],source = "popAncB", dest = "popAncAB", proportion = 1)
			demography.add_mass_migration(time=options["generation_split_OA"],source = "popAncAB", dest = "popAncO", proportion = 1)
			demography.add_mass_migration(time=options["generation_split_OA"],source = "popO", dest = "popAncO", proportion = 1)
		# Divergence will be made by split option in Msprime
		elif (coalescence_split == "split") : 
			demography.add_population_split(time=options["generation_split_B"], derived=["popB1","popB2"], ancestral="popAncB")
			demography.add_population_split(time=options["generation_split_A"], derived=["popA1","popA2"], ancestral="popAncA")
			demography.add_population_split(time=options["generation_split_AB"], derived=["popAncA","popAncB"], ancestral="popAncAB")
			demography.add_population_split(time=options["generation_split_OA"], derived=["popO","popAncAB"], ancestral="popAncO")
			
			
		# simulate with msprime
		# https://tskit.dev/msprime/docs/latest/intro.html
		msp_ts = msprime.sim_ancestry(samples            = {"popO" : options["popO"],"popA1" : options["popA1"],"popB1" : options["popB1"],"popA2" : options["popA2"],"popB2" : options["popB2"]},
				        demography         = demography,
				        model              = "dtwf", # because sample size = pop size
				        recombination_rate = rate_map,
				        random_seed        = options["seed_coal"])
		
		# make tree-sequence a SLiM-tree-sequence
		slim_ts = pyslim.annotate(msp_ts, model_type = "WF", tick = 1)
		
		# save tree-sequence
		slim_ts.dump(job_dir+"/coalsim_"+options["analysis"]+"_"+options["project"]+"_sim_"+options["sim"]+".trees")
		
	#If the user wants the backward simulations to go as far as the divergence between the ancestral populations of the donor and recipient populations. 
	elif (slim_split == 1 ) :
		
		# get recombination map:
		rate_map = msprime.RateMap(position = options["msprime_r_map"]["positions"],
				     rate = options["msprime_r_map"]["rates"]) 
		
		# get demography
		demography = msprime.Demography()
		demography.add_population(name = "popAncO", initial_size = options["popAncO"]) #p0
		demography.add_population(name = "popO", initial_size = options["popO"]) #p1
		demography.add_population(name = "popAncAB", initial_size = options["popAncAB"]) #p2
		demography.add_population(name = "popAncA", initial_size = options["popAncA"]) #p3
		demography.add_population(name = "popAncB", initial_size = options["popAncB"]) #p4
		
		if (coalescence_split == "mass_migration") :
			demography.add_mass_migration(time=options["generation_split_AB"],source = "popAncA", dest = "popAncAB", proportion = 1)
			demography.add_mass_migration(time=options["generation_split_AB"],source = "popAncB", dest = "popAncAB", proportion = 1)
			demography.add_mass_migration(time=options["generation_split_OA"],source = "popAncAB", dest = "popAncO", proportion = 1)
			demography.add_mass_migration(time=options["generation_split_OA"],source = "popO", dest = "popAncO", proportion = 1)
		elif (coalescence_split == "split") : 
			demography.add_population_split(time=options["generation_split_AB"], derived=["popAncA","popAncB"], ancestral="popAncAB")
			demography.add_population_split(time=options["generation_split_OA"], derived=["popO","popAncAB"], ancestral="popAncO")
		
		
		# simulate with msprime
		# https://tskit.dev/msprime/docs/latest/intro.html
		msp_ts = msprime.sim_ancestry(samples            = {"popO" : options["popO"], "popAncA" : options["popAncA"], "popAncB" : options["popAncB"]},
				        demography         = demography,
				        model              = "dtwf",
				        recombination_rate = rate_map,
				        random_seed        = options["seed_coal"])
		
		# make tree-sequence a SLiM-tree-sequence
		slim_ts = pyslim.annotate(msp_ts, model_type = "WF", tick = 1)
		# slim_ts = pyslim.annotate(msp_ts, model_type = "WF", slim_generation = 1)
		# save tree
		slim_ts.dump(job_dir+"/coalsim_"+options["analysis"]+"_"+options["project"]+"_sim_"+options["sim"]+".trees")

		msprime_end = time.time()
		duration_msprime = msprime_end - msprime_start
		print("(msprime) backward simulation duration : "+str(duration_msprime)+" s ")

############################################################################################################
############################################################################################################
if __name__ == "__main__":
	principal()
