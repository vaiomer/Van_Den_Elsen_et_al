#!/usr/bin/env python2.7
# -*- coding: utf-8 -*

# Original script by Koichi Higashi
# Found at : https://github.com/khigashi1987/GUniFrac/blob/master/GUniFrac.py
# This version of the script has been modified to improve execution time.

from optparse import OptionParser
import numpy as np
import pandas as pd
from ete2 import Tree
import itertools
from scipy.spatial.distance import squareform

def importData(filename):
	table = pd.read_table(filename,index_col=0)
	relative_abundance = preprocessing.relative_abundance(table.values, axis=0)
	new_table = pd.DataFrame(relative_abundance, index=table.index, columns=table.columns)
	return new_table

def compute_GUniFrac(abundance,treefile, alpha=0.5):
	n_samples = len(abundance.columns)
	n_distance = n_samples * (n_samples - 1) / 2
	d_array = np.zeros((n_distance))
	t = Tree(treefile,format=1)
	if set(t.get_leaf_names()) < set(abundance.index):
		print 'Error: OTU table contains unknown OTUs. All of OTU names in OTU table should be contained in tree file.'
		quit()
	#List that will contain all rows of the dataframe. Each row will contain every values calculated at a node.
	alldata = []
	
	i = 0
	for node in t.traverse():
		if node.is_root():
			continue
		else:
			nodevalues = {}			
			for leaf in node.get_leaf_names():
				if leaf in abundance.index:
					if "nodename" not in nodevalues.keys() and "nodedist" not in nodevalues.keys():
						nodevalues["nodename"] = "node"+str(i)
						nodevalues["nodedist"] = node.dist
					for sample in abundance.columns:
						if sample not in nodevalues.keys():
							nodevalues[sample] = 0
							nodevalues[sample] += abundance.loc[leaf,sample]
						else:
							nodevalues[sample] += abundance.loc[leaf,sample]
			i+=1
		if nodevalues != {}:
			alldata.append(nodevalues)	
	alldata = pd.DataFrame(alldata)
	alldata.set_index("nodename", True, inplace = True)  
	for i,(sample1, sample2) in enumerate(itertools.combinations(abundance.columns, 2)):
		numer = 0
		denom = 0
		for node in alldata.index.values:
			p_a = alldata.loc[node,sample1]
			p_b = alldata.loc[node,sample2]
			if p_a == 0 and p_b == 0:
				continue
			denom += alldata.loc[node,"nodedist"]*(p_a+p_b)**alpha
			numer += alldata.loc[node,"nodedist"]*(p_a+p_b)**alpha*abs(p_a-p_b)/(p_a+p_b)
		d_array[i] = numer/denom
	return squareform(d_array)

if __name__ == '__main__':
	usage = 'usage: python %prog [options]'
	parser = OptionParser(usage)
	parser.add_option( "-f", "--file", action="store", dest="data_file", help="TAB-separated OTU table. rows are OTUs, columns are samples.")
	parser.add_option( "-t", "--tree", action="store", dest="tree_file", help="Rooted phylogenetic tree. NEWICK format.")
	parser.add_option( "-a", "--alpha", action="store", type="float", dest="alpha", default=0.5, help="alpha parameter of generalized UniFrac metric.")
	options, args = parser.parse_args()
	if options.data_file == None or options.tree_file == None:
		print "ERROR: requires options"
		parser.print_help()
		quit()
	datafile = options.data_file
	treefile = options.tree_file
	alpha = options.alpha
	
	abundance = importData(datafile)
	distance_matrix = compute_GUniFrac(abundance, treefile, alpha=alpha)
	distance_table = pd.DataFrame(distance_matrix, index=abundance.columns, columns=abundance.columns)
	distance_table.to_csv('./result.dist',sep='\t')
