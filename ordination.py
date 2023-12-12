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
import re
import matplotlib
matplotlib.use('svg')
from skbio.stats.ordination import pcoa
import utilities as ut
import plotnine as pn


##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################
def get_ordination(biom, input_distance_matrix):
	"""
	@summary: Removes filtered samples from distance matrix and performs PCoA ordination on it.
	@param biom: [Biom] Biom object
	return: [OrdinationResults] scipy OrdinationResults type object containing data from the ordination such as eigenvalues and eigenvectors.
	"""
	#Load and modify distance matrix to eliminate filtered samples and change it to condensed form.
	otutab = ut.get_otu_table(biom)
	dist_data = pd.read_csv(input_distance_matrix, sep = '\t', header = 0, index_col = 0)
	dist_data.index = dist_data.index.astype(str)
	
	#Perform ordination
	ord_data = pcoa(dist_data)
	
	return ord_data
	
def plot_ordination (biom, comp_name, ord_data, distance_measure):
	"""
	@summary: Plots ordination and write it to output image file.
	@param biom: [Biom] Biom object
	@param ord_data: [OrdinationResults] OrdinationResults object.
	"""
	#Get sample metadata.
	metadata = ut.get_metadata_sample(biom)
	metadata = pd.DataFrame(metadata)
	metadata.set_index("Sample_name", drop = False, inplace = True)
	del metadata.index.name
	
	#Add ordination values for axes 1 and 2 to the metadata dataframe.
	metadata["Axis.1"] = ord_data.samples["PC1"].values 
	metadata["Axis.2"] = ord_data.samples["PC2"].values
	#Plot values.
	group_order = ut.get_group_order(biom, comp_name)
	metadata[comp_name] = pd.Categorical(metadata[comp_name],group_order)
	ord_map = pn.aes(x = "Axis.1", y = "Axis.2", **{"color" : "Group"})
	p = pn.ggplot(metadata, ord_map) + pn.geom_point()
	
	#Get relative eigenvalues for axes names.
	function = lambda x: x/sum(x)
	eig = pd.DataFrame(ord_data.eigvals)
	rel_eig = eig.apply(function)
	ax1_string = "Axis.1 ["+ str(round(rel_eig.iloc[0][0]*100, 1))+"%]"
	ax2_string = "Axis.2 ["+ str(round(rel_eig.iloc[1][0]*100, 1))+"%]"
	p = p + pn.xlab(ax1_string) + pn.ylab(ax2_string) + pn.theme_bw() + pn.scale_color_manual(values = {"Colostrum":"#000000", "Mature":"#dc3912"})
	
	#Set title
	title = distance_measure+" + MDS"
	#Modify title if input_distances is gunifrac.
	if "alpha" in distance_measure:
		alpha = re.sub("^.+alpha([\.0-9]*[0-9]).+$","\\1", distance_measure) #Get the value after alpha
		
	p = p+pn.ggtitle(title)
	
	#Add ellipses
	p = p+ pn.stat_ellipse(level = 0.95)
		
	#Save to file.	
	if not os.path.exists("betaDiversityOrdination"):
		os.mkdir("betaDiversityOrdination")
	p.save(os.path.join("betaDiversityOrdination", "betaDiversity_ordination_"+ distance_measure +".svg"), width = 6, height = 4)

##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################

if __name__ == "__main__":
	#Load Biom
	input_biom = "PC0121.biom"
	comp_name = "Group"
	biom = ut.load_and_filter_biom(input_biom, comp_name)
	
	#calculate ordination for each distance matrix
	distance_matrix_directory = "betaDiversityMatrices"
	for matrix_file in os.listdir(distance_matrix_directory):
		if matrix_file.endswith(".tsv"):
			distance_measure = re.sub("betaDiversity_(.+)\.tsv","\\1",matrix_file)
			ord_data = get_ordination(biom, os.path.join(distance_matrix_directory, matrix_file)) #ord_data is a skbio.stats.ordination.OrdinationResults type object (see doc) 
			plot_ordination(biom, comp_name, ord_data, distance_measure)


