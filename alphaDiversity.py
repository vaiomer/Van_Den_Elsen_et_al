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
import matplotlib
import pandas as pd
import skbio.diversity as sd
import scipy.stats.mstats as sm
import scipy.stats as ss
import re
import utilities as ut
from skbio import TreeNode as st
from frogsBiom import BiomIO, Biom



def estimate_richness(otu_table, rank, group_order):
	div = []
	array = otu_table.values.astype(int) #Turn content of otutable into  nparray to int type. Needed with skbio.diversity.
	array = array.transpose() #Array needs to be transposed.
	
	simp = sd.alpha_diversity("simpson", array, ids = otu_table.columns.values)
	div.append(sd.alpha_diversity("simpson", array, ids = otu_table.columns.values).round(3))
	
	div.append(simp.apply(calculate_invsimpson).round(3))
		
	divdf = pd.concat(div, axis = 1) #concat series into a single dataframe.
	divdf.columns = ["Simpson", "InvSimpson"] #Rename columns.
	
	groups = ut.get_groups(biom, comp_name) #Get dictionnary containing each sample and its group.
	divdf = pd.concat([divdf,pd.Series(groups, name = "Groups")], axis = 1) #Add group column to the diversity_measures dataframe
	divdf["Groups"] = pd.Categorical(divdf["Groups"],group_order)
	divdf = divdf.sort_values("Groups")
	
	write_richness_to_file(divdf, rank)
	
	return divdf


def calculate_invsimpson(simpsonvalue):
	invsimpsonvalue = 1/(1-simpsonvalue)
	return invsimpsonvalue
	

def write_richness_to_file (measures_dataframe, rank):
	filename = re.sub("(\.[^\.]+)$", "_"+rank+"\\1", output_result)
	measures_dataframe.to_csv(filename, sep = "\t")


def kruskal_test(biom, diversity_measures,rank, kruskalwriter):
	diversity_measures = diversity_measures.set_index('Groups')
	kw_results = {}
	wilcox_pairwise_results = {}
	for indice in diversity_measures.columns:
		liste = [] #Will contain an array with alpha diversity for index "indice" of each group.
		smallgroup = []
		for group in set(diversity_measures.index): #Select each group for an indice
			groupdf = diversity_measures[diversity_measures.index == group]
			liste.append(groupdf[indice].values) #add each group dataframe to a list.			
			kw_results[indice] = sm.kruskal(*liste) #Perform the kruskal comparing each group (for one alpha index) and add it to the kwresults dict
			
		for i in range(0, len(set(diversity_measures.index))): #These loops get every pair of groups, do wilcoxon ranksums test on them and write them in the wilcox_pairwise_results dict. The key is ("AlphaIndice", "Group1", "Group2")
			group_i = list(set(diversity_measures.index))[i] #Name of group i
			if len(diversity_measures.loc[[group_i]]) < 5 and group_i not in  smallgroup: #Check if a group is a smallgroups (<5 samples) and add it to the smallgroup list if needed.
				smallgroup.append(group_i) #Smallgroups will later get a message indicating the kruskal and wilcoxon test aren't precise
			if i < len(set(diversity_measures.index)):
				for j in range(i+1, len(set(diversity_measures.index))): 
					group_j = list(set(diversity_measures.index))[j]
					if (group_i, group_j) not in wilcox_pairwise_results.keys():
						wilcox_pairwise_results[(group_i, group_j)] = {}					
					pairlist = [] #This list contains dataframes of pairs of groups to run wilcox test on them.
					pairlist.append(diversity_measures.loc[[group_i],[indice]])
					pairlist.append(diversity_measures.loc[[group_j],[indice]])
					wilcox_pairwise_results[(group_i, group_j)][indice] = ss.ranksums(*pairlist)

	write_kruskal_to_file(rank, smallgroup, kruskalwriter, kwresults = kw_results)	#Write the results of the kw test to output file.
	#Ecrire resultats posthoc si pvalue < 0.05
	excel_row = 12
	for kwkey in kw_results.keys():
		if kw_results[kwkey][1] < 0.05 : #If a signicative difference is found with the kw test for at least one alpha index, write the pairwise results to the output file.
			write_kruskal_to_file(rank, smallgroup, kruskalwriter, pairwiseresults = wilcox_pairwise_results, excel_row = excel_row)
			break
	
	

def write_kruskal_to_file(rank, smallgrouplist, kruskalwriter, kwresults = None, pairwiseresults = None, excel_row = None):
	indices = ["Simpson", "InvSimpson"]
	if rank != "OTU" and "PD" in indices: #PD is only calculated at OTU rank for now.
		indices.remove("PD")
	if kwresults != None :	
		kruskaldf = pd.DataFrame(kwresults, index = ["Kruskal-Wallis Statistic", "P-Value"], columns = indices).transpose()
		significativitylist = []
		for ind in kruskaldf.index:
			if kruskaldf.at[ind,'P-Value'] <= 0.001 :
				significativitylist.append('***')
			elif kruskaldf.at[ind,'P-Value'] <= 0.01 :
				significativitylist.append('**')
			elif kruskaldf.at[ind,'P-Value'] <= 0.05 :
				significativitylist.append('*')
			else :
				significativitylist.append(' ')
		kruskaldf['Significativity'] = significativitylist
		kruskaldf.to_excel(kruskalwriter, sheet_name= rank, startrow = 3)
		usedsheet = kruskalwriter.sheets[rank]
		### set column sizes
		usedsheet.column_dimensions['A'].width = 25
		usedsheet.column_dimensions['B'].width = 25
		usedsheet.column_dimensions['C'].width = 25
		usedsheet.column_dimensions['D'].width = 25
		usedsheet.cell(row = 1, column = 1).value = 'Significativity codes :'
		usedsheet.cell(row = 1, column = 2).value = '≤0.001 : ***'
		usedsheet.cell(row = 1, column = 3).value = '≤0.01 : **'
		usedsheet.cell(row = 1, column = 4).value = '≤0.05 : *'
		if smallgrouplist != []:
			usedsheet.cell(row = 1, column = 5).value = 'Warning (!) : 	Kruskal Wallis and Wilcoxon ranksums tests are inaccurate when used on groups with less than 5 samples. The symbol (!) will be used to show if a test is innacurate.'
			usedsheet.cell(row = 2, column = 5).value = 'The following groups contain less than 5 samples : '+ ', '.join(smallgrouplist)+'. '
			usedsheet.cell(row = 3, column = 2).value = '(!)'
			
		usedsheet.cell(row = 3, column = 1).value = 'All Comparisons'
		
	if pairwiseresults != None :				
		for pair in pairwiseresults.keys():
			group1 = pair[0]
			group2 = pair[1]
			wilcoxondf = pd.DataFrame(pairwiseresults[pair], index = ["Wilcoxon Statistic", "P-Value"], columns = indices).transpose()
			# ~ wilcoxondf = wilcoxondf.round(4)
			significativitylist = []
			for ind in wilcoxondf.index :
				if wilcoxondf.at[ind,'P-Value'] <= 0.001 :
					significativitylist.append('***')
				elif wilcoxondf.at[ind,'P-Value'] <= 0.01 :
					significativitylist.append('**')
				elif wilcoxondf.at[ind,'P-Value'] <= 0.05 :
					significativitylist.append('*')
				else :
					significativitylist.append(' ')
			wilcoxondf['Significativity'] = significativitylist
			wilcoxondf.to_excel(kruskalwriter, sheet_name = rank, startrow = excel_row)
			usedsheet = kruskalwriter.sheets[rank]
			usedsheet.cell(row = excel_row, column = 1).value = group1+" VS "+group2
			if group1 in smallgrouplist or group2 in smallgrouplist :
				usedsheet.cell(row = excel_row, column = 2).value = "(!)"
			excel_row += 9


##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################

if __name__ == "__main__":
	
	# Inputs
	input_biom = "PC0121.biom"
	comp_name = "Group"
	selected_ranks = "Genus"
	output_result = "alpha_diversity.tsv"
	output_kruskal = "alpha_diversity_kruskal.xlsx"
	taxrank = "Genus"
	
	# Load biom_file
	biom = ut.load_and_filter_biom(input_biom, comp_name)
	
	# Alpha diversity
	group_order = ut.get_group_order(biom, comp_name)
	rank_otu_table = ut.merge_by_rank(biom, taxrank, ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"])
	diversity_measures = estimate_richness(rank_otu_table, taxrank, group_order)
	
	# Kruskal-Wallis test
	kruskalwriter = pd.ExcelWriter(output_kruskal)
	kruskal_test(biom, diversity_measures, taxrank, kruskalwriter)
	kruskalwriter.save()
