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
#import operator
import pandas as pd
import re
from frogsBiom import BiomIO, Biom

##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################

def load_and_filter_biom(biomfile, comp_name, selected_groups = None, discarded_groups = None, min_spl_seq = None, max_spl_seq = None):
	"""
	@summary: Load biom object from a .biom file and filter it as specified in arguments.
	@param biomfile: Path to the .biom file to load and filter.
	@return: [biom] The biom object.
	"""
	##Load file
	biom = BiomIO.from_json(biomfile)
	
	#Verify loaded biom contains comp_name metadata
	if Biom.has_metadata(biom, comp_name, "sample", True) == False :
		sys.exit("The BIOM file does not contain metadata for specified comp_name argument (default value : 'Group')")
	
	##Filter by groups
	#Some scripts don't contain each of these filter arguments. The try/except is done to ignore errors when that is the case.
	try:
		biom = filter_samples_by_meta(biom, comp_name, selected_groups, discarded_groups)
	except AttributeError: 
		pass	
	#Filter by count	
	try:
		biom = filter_samples_by_count(biom, min_spl_seq, max_spl_seq)
	except AttributeError: 
		pass
	return biom
	
def filter_samples_by_meta(biom, comp_name, selected_groups = None, discarded_groups = None):
	"""
	@summary: Keep only selected_groups or remove discarded_groups in Biomfile.
	@param biom: biom object to filter.
	@param selected_groups: List of selected sample groups to keep.
	@param discarded_groups: List of sample groups to discard.
	@return: [biom] The filtered  biom object.
	"""
	#Obtain sample names
	names = []
	discard = []
	for sname in Biom.get_samples_names(biom): #get_samples_names yields a generator. Multiple iterations are needed to get all names.
		names.append(sname)	
	
	#Remove samples in discarded groups		
	if discarded_groups != None :	
		for dsname in names :
			metadat = Biom.get_sample_metadata(biom, dsname)
			for dgroup in discarded_groups :
				if dgroup ==  metadat[comp_name]: #Check if sample is part of a discarded group.
					discard.append(dsname)
		Biom.remove_samples(biom, discard)
	
	#Remove samples not in selected groups
	elif selected_groups != None :
		for ssname in names :
			remove = True
			metadat = Biom.get_sample_metadata(biom, ssname)
			for sgroup in selected_groups :
				if sgroup == metadat[comp_name] : #If samples group is a selected group, it will not be removed.
					remove = False
			if remove == True :
				discard.append(ssname)
				
		Biom.remove_samples(biom, discard)
	return biom
	
def filter_samples_by_count(biom, min_count, max_count):
	"""
	@summary: Remove samples with less counts than min_count and more counts than max_count.
	@param biom: biom object to filter.
	@param min_count: Minimum accepted number of counts.
	@param max_count: Maximum accepted number of counts.
	@return: [biom] The filtered  biom object.
	"""
	discard = []
	for sname in Biom.get_samples_names(biom):
		if (min_count != None and Biom.get_sample_count( biom, sname) <= min_count) or (max_count != None and Biom.get_sample_count( biom, sname) > max_count) :
			discard.append(sname)
	Biom.remove_samples(biom, discard)
	return biom

def get_otu_table(biom) :
	"""
	@summary: Extract OTU table from biom object and return it as a pandas dataframe.
	@param biom: biom object containing the OTU table.
	@return: [df_otu_table] The OTU table.
	"""
	otu_table_generator = Biom.to_count_table(biom)
	otu_table = []
	for otu in otu_table_generator :
		otu_table.append(otu)
	df_otu_table = pd.DataFrame(otu_table, columns = otu_table[0]) #Preciser que première ligne = noms des colonnes
	df_otu_table.drop(0, axis = 0, inplace = True)  #Supprimer première ligne. inplace = True permet de faire le changement directement dans df.
	df_otu_table = df_otu_table.set_index('#Observation') #Nom des lignes correspond à la colonne #Observation.	
	return df_otu_table

def get_all_observations_taxonomy(biom):
	names = []
	obstax = []
	for obs in Biom.get_observations_names( biom ):
		names.append(obs)
	taxkey = dict.keys(Biom.get_observation_metadata(biom, names[0]))[0] #Gives key for taxonomy in metadata, needed for get_observation_taxonomy.
	for obs in names :
		obstax.append(Biom.get_observation_taxonomy( biom, obs, taxkey )) #List containing taxonomy of each OTU
		#["K","P","C","O","F","G","sp"]
	return names,obstax
	
def merge_by_rank(biom, rank, tax_ranks, checkprevious = True):
	"""
	@summary: Merge all OTU with the same taxonomic rank.
	@param biom: biom object with the OTU metadata on taxonomy.
	@param rank: Taxonomy rank.
	@param checkprevious: True -> Merge only if OTU also share previous taxonomic rank. 
	@return: [df_otu_table] The OTU table with all merges done.
	"""
	names = []
	obstax = []
	df_otu_table = get_otu_table(biom)
	rankidx = tax_ranks.index(rank)
	for obs in Biom.get_observations_names( biom ):
		names.append(obs) #List of all OTU names
	taxkey = "taxonomy" #Name for taxonomy in metadata, needed for get_observation_taxonomy.
	for obs in names :
		obstax.append(Biom.get_observation_taxonomy( biom, obs, taxkey )) #List containing taxonomy of each OTU
		#["K","P","C","O","F","G","sp"]
	if checkprevious == True:
		i = 0
		while i < len(obstax):
			j = i + 1
			while j < len(obstax):

				if obstax[i][rankidx] == obstax[j][rankidx] and obstax[i][rankidx-1] == obstax[j][rankidx-1] and obstax[i][rankidx]!= "Multi-affiliation" and obstax[i][rankidx]!= "Unknown" : #Check if any OTU j has the same [rank] as OTU i.

					df_otu_table = merge_otu(df_otu_table, names[i], names[j]) #If it does, merge them
					del obstax[j] #Then delete j OTU from taxonomy and name lists.
					del names[j]
				
				elif obstax[i][rankidx]== "Multi-affiliation" and obstax[i][rankidx] == obstax[j][rankidx]: #if OTU i and j are multiaffiliations

					for k in range(rankidx, -1, -1):	#Only merge OTU if previous non multiaffiliated rank is same for i and j (from species to kingdom)		
						if obstax[i][k] != "Multi-affiliation" and obstax[i][k] == obstax[j][k]:
							df_otu_table = merge_otu(df_otu_table, names[i], names[j])
							del obstax[j] 
							del names[j]
							break
						elif obstax[i][k] != obstax[j][k]:
							j = j+1
							break
				elif (obstax[i][rankidx]== "Unknown" or "unknown" in obstax[i][rankidx]) and obstax[i][rankidx] == obstax[j][rankidx]:
					for k in range(rankidx, -1, -1):
						if (obstax[i][k] != "Unknown" and "unknown" not in obstax[i][k]) and obstax[i][k] == obstax[j][k]:
							df_otu_table = merge_otu(df_otu_table, names[i], names[j])
							del obstax[j] 
							del names[j]
							break
						elif obstax[i][k] != obstax[j][k]:
							j = j+1
							break
				else :

					j = j+1
			df_otu_table.rename(index = {names[i]:obstax[i][rankidx]}, inplace=True) #Rename rows as the observed taxonomic rank. Not an obligation, done for clarity during testing.
			i += 1	
	else :
		i = 0
		while i < len(obstax):
			j = i + 1
			while j < len(obstax):
				if obstax[i][rankidx] == obstax[j][rankidx]:
					df_otu_table = merge_otu(df_otu_table, names[i], names[j]) #If it does, merge them
					del obstax[j] #Then delete j OTU from taxonomy and name lists.
					del names[j]

				else :
					j = j+1
			i += 1
						
	return df_otu_table

def merge_otu(df_otu_table, otuname1, otuname2):
	"""
	@summary: Merge the two selected OTU and sum their counts for each sample. Note : the first otu is kept, the second one is removed.
	@param df_otu_table: OTU table where the merge is done.
	@param otuname1: Name of the first OTU to merge.
	@param otuname2: Name of the second OTU to merge.
	@return: [df_otu_table] The OTU table with the specified merge done.
	"""
	df_otu_table.loc[otuname1] += df_otu_table.loc[otuname2] #add otu2 to otu 1
	df_otu_table.drop([otuname2], inplace = True) #drop otu2
	return df_otu_table

def color_set(numberofcolors):
	"""
	@summary: Returns a list containing as many different colors as specified
	@param numberofcolors: Number of colors to return.
	@return: [colorlist]: List containing the required number of different colors.
	
	"""
	colorlist = []
	allcolors = ["#3366CC", "#DC3912", "#FF9900", "#109618", "#990099", "#0099C6", "#DD4477", "#66AA00", "#B82E2E", "#316395", "#FFBB78", "#FF9896", "#C5B0D5", "#98DF8A", "#DBDB8D", "#7ac5cd", "#ee6a50", "#6495ed", "#d2691e", "#ee82ee", "#7fff00", "#ee3b3b", "#9932cc", "#9bcd9b", "#00bfff", "#ff6347", "#00ff00", "#8b8b00", "#b9d3ee", "#dda0dd", "#7fffd4", "#ffdab9", "#ff3030", "#b0e0e6", "#a52a2a", "#bdb76b", "#40e0d0", "#ffc0cb", "#cd5c5c", "#CCFF00FF", "#00FF66FF", "#0066FFFF", "#FF00E6FF", "#B3FFFFFF", "#E6FFFFFF", "#FFD9FFFF", "#FF99FFFF", "#bf3eff", "#00ced1", "#cdcd00", "#c0ff3e", "#ab82ff", "#ffa07a", "#c71585", "#ff34b3", "#00ffff", "#8b7355", "#00cd00", "#ffe4e1", "#ffff00", "#ff0000", "#000080", "#006400", "#ffa500", "#0000ff", "#ff1493", "#8b008b", "#ffd700", "#556b2f", "#ee7942", "#d8bfd8", "#ff3e96", "#ffb90f"]
	i = 0
	while i < numberofcolors: #Add a color from allcolors to the list for each group
		try: colorlist.append(allcolors[i]) #In the case where there are more groups than different colors : loop back to the beginning.
		except IndexError : 
			numberofcolors = numberofcolors-i
			i = 0		
			colorlist.append(allcolors[i])
		i += 1
	return colorlist

def get_group_colors(biom, groups_colors, comp_name):
	"""
	@summary: Returns a dictionnary associating each group in the biom object with a color for plotting.
	@param biom: biom object containing sample metadata on group.
	@return: [color_dict] Dictionnary containing each group's name (keys) and an associated color (values)
	
	"""
	color_dict = {}
	groupdict = get_groups(biom, comp_name)
	grouplist = []
	for key in groupdict.keys():
			if groupdict[key] not in grouplist:
				grouplist.append(groupdict[key])
	if groups_colors == None:		
		used_colors = color_set(len(grouplist))
		for i in range (0,len(grouplist)):
			color_dict[grouplist[i]] = used_colors[i]
	else:
		colordf = pd.read_csv(groups_colors, sep = "\t", names = ["Groups","Colors"], index_col = "Groups")
		color_dict = colordf.to_dict()['Colors'] #colordf dictionnary initially contains one key "Colors" and a value which is the dictionnary containing each group as a key and colors as a value.
		#As we only want that second dictionnary, we select the value of the key 'Colors'.
		missing_groups = False
		missing_groups_list = []
		for group in grouplist:
			if group not in color_dict.keys():
				missing_groups = True
				missing_groups_list.append(group)
		if missing_groups is True:		
			raise Exception('There are groups missing in file '+ groups_colors +' : '+ str(missing_groups_list))
	return color_dict

def get_groups(biom, comp_name):
	"""
	@summary: Returns dictionnary containing each sample(keys) and their groups (values)
	@param biom: biom object containing sample metadata on group.
	@return: [samplegroup] Dictionnary containing samples and their groups.
	
	"""
	samplegroup = {}
	for sname in Biom.get_samples_names(biom):
		metadat = Biom.get_sample_metadata(biom, sname)
		samplegroup[sname] = metadat[comp_name] #Assign its group to each sample.
	return samplegroup	

def get_metadata_sample(biom):
	names = []
	for sname in Biom.get_samples_names(biom): #get_samples_names yields a generator. Multiple iterations are needed to get all names.
		names.append(sname)	
	metadata = []
	for dsname in names :
		sampledata = Biom.get_sample_metadata(biom, dsname)
		sampledata["Sample_name"] = dsname
		metadata.append(Biom.get_sample_metadata(biom, dsname))
		
	return metadata

def filter_distance_matrix(otutab, dist_data, output_path = None):
	for sample in dist_data.columns.values:
		if sample not in otutab.columns.values:
			dist_data.drop(columns = sample, inplace = True)
			dist_data.drop(index = sample, inplace = True)
	if output_path is not None:
		dist_data.to_csv(path_or_buf = output_path, sep = '\t')
	return dist_data
	
def get_group_order(biom, comp_name):
	grouplist = list(set(get_groups(biom, comp_name).values())) #List of all groups in biom after filtering.
	if biom.comment is not None and comp_name in biom.comment:
		comp_name_regex = re.escape(comp_name)
		group_order_str=re.sub(".*__"+comp_name_regex+"_group_order=([^\s]*)__.*", "\\1",biom.comment) #Get a character string containing each group in order, separated by a ";"
		group_order_list = group_order_str.split(";") #Split into a list
		grouplist = list(set(get_groups(biom, comp_name).values())) #List of all groups in biom after filtering.
		if not set(grouplist).issubset(group_order_list):
			raise Exception("ERROR: There are groups missing in the biom comment field. Please add missing groups :"+str(set(grouplist).difference(set(group_order_list))))
		removed_groups = []
		for group in group_order_list:
			if group not in grouplist:
				removed_groups.append(group)
		for group in removed_groups :
			group_order_list.remove(group)
	else :
		group_order_list = sorted(grouplist)
	return group_order_list
