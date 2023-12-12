#!/usr/bin/env python2.7
# -*- coding: utf-8 -*
#
# Copyright (C) 2019 Vaiomer
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

__author__ = 'Sebastian Van Blerk - Vaiomer'
__copyright__ = 'Copyright (C) 2019 Vaiomer'
__license__ = 'GNU General Public License'
__version__ = '0.1'
__email__ = 'bioinfo@vaiomer.com'

import os
import sys
import pandas as pd
import skbio.diversity as sd
from skbio import TreeNode as st
from skbio.stats.distance import permanova, DistanceMatrix

import utilities as ut
import numpy as np
import GUniFrac2 as GUniFrac

##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################
def get_distance_matrix(biom, input_tree):

	#Create output directories
	if not os.path.exists("betaDiversityMatrices"):
		os.mkdir("betaDiversityMatrices")
	if not os.path.exists("betaDiversityPermanova"):
		os.mkdir("betaDiversityPermanova")
	
	#Get sample groupings
	groupingdf = ut.get_groups(biom, "Group")
	groupingdf = pd.DataFrame.from_dict(groupingdf, orient="index", columns=["Group"])
	
	#Get count ratio
	otutab = ut.get_otu_table(biom)
	ratio_otutab = otutab.apply(lambda x: (x*10000000)/sum(x), 0) #Multiplied by a big number to turn floats into int to bypass floating point errors in this version of scikit-bio. Doesn't change the resulting distances.
	ratio_otutab_array = ratio_otutab.values.astype(np.int64).transpose() #Transform to int array for skbio
	
	#Calculate Jaccard distance
	bool_otutab_array = np.array(ratio_otutab_array, dtype = bool) #Change abundance matrix to boolean.
	dist_data = sd.beta_diversity("jaccard", bool_otutab_array, ids = otutab.columns.values)
	dist_datadf = dist_data.to_data_frame()	
	dist_datadf.to_csv(path_or_buf = "betaDiversityMatrices/betaDiversity_jaccard.tsv", sep = "\t")
	compute_permanova(dist_data, groupingdf, permanova_output = "betaDiversityPermanova/jaccard_permanova.tsv")
	
	#Calculate Bray-Curtis distance
	dist_data = sd.beta_diversity("braycurtis", ratio_otutab_array, ids = otutab.columns.values)
	dist_datadf = dist_data.to_data_frame()	
	dist_datadf.to_csv(path_or_buf = "betaDiversityMatrices/betaDiversity_bray.tsv", sep = "\t")
	compute_permanova(dist_data, groupingdf, permanova_output = "betaDiversityPermanova/bray_permanova.tsv")
	
	#Calculate GUniFrac distance with alpha=0
	df_ratio_otutab_array = pd.DataFrame(ratio_otutab_array.transpose(),  index = otutab.index, columns = otutab.columns)
	dist_data = GUniFrac.compute_GUniFrac(df_ratio_otutab_array, input_tree, alpha = 0)
	dist_datadf = pd.DataFrame(data = dist_data, index = otutab.columns.values, columns = otutab.columns.values)
	dist_datadf.to_csv(path_or_buf = "betaDiversityMatrices/betaDiversity_gunifrac_alpha0.tsv", sep = "\t")
	dist_matrix = DistanceMatrix.read("betaDiversityMatrices/betaDiversity_gunifrac_alpha0.tsv")#NOTE : GUniFrac doesn't return a distance matrix object as needed by the permanova function. We need to reload the saved file as a distance matrix.
	compute_permanova(dist_matrix, groupingdf, permanova_output = "betaDiversityPermanova/gunifrac_alpha0_permanova.tsv")

def compute_permanova(dist_data, groupingdf, permanova_output):
	#Compute permanova and write to file
	permanova_result = permanova(distance_matrix = dist_data, grouping = groupingdf, column = "Group", permutations = 2000)
	permanova_result.to_csv(permanova_output, sep="\t", header=True)
	
##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################

if __name__ == "__main__":
	
	# Manage parameters
	input_biom = "PC0121.biom"
	comp_name = "Group"
	input_tree = "PC0121.nwk"
		
	biom = ut.load_and_filter_biom(input_biom, comp_name)
	dist_data = get_distance_matrix(biom, input_tree)
