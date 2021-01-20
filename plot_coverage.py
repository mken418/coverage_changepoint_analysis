#!/usr/local/apps/bioapps/python/Python-3.5.3/bin/python3


import argparse
import os
import re
import random
from pprint import pprint
from collections import OrderedDict

import matplotlib.pyplot as plt



"""
Plot the perbase coverage values to better visualize the variant
Designed for a maximum of 25 individuals and will cut off at 25 individuals. 
Break up the input files into chunks of 25 you have more individuals per variant to plot. 
Min and max coverage is on each y axis
"""

parser = argparse.ArgumentParser()
#parser.add_argument('--input_sig_loci', required=True, help='file with parsed genes of most significant variants')
parser.add_argument('--variant_locus', required=True, help='locus in the form of example: chr1:10-100')
parser.add_argument('--coverage_output', required=True, help='dir with per base coverage output from sentieon')
parser.add_argument('--matlab_output', required=True, help='file for the matlab signal detection output')
parser.add_argument('--outfile', required=True, help='dir for output plots to go')
args = parser.parse_args()


def plot_coverage(var_locus, coverage_output, matlab_output, outfile):
	"""
	plots coverage from sentieon covg metrics output, given a line split into list from the input_sig_loci file
	line looks like: id      chr     start   stop    probands        parents biobank gene_info       p_values        fdr_adjusted
	"""
	match = re.match(r'(.+):(.+)-(.+)', var_locus)
	chrom, start, stop = match.group(1), match.group(2), match.group(3)
	var_coords = [int(start), int(stop)]
	var_size = abs(var_coords[1] - var_coords[0])


	#get coverage range from covg_metrics file (just check the first file)	
	covg = []
	for sample in os.listdir(coverage_output):
		fh = open(coverage_output+'/'+sample, 'r')
		for line in fh:
			if line.startswith('Locus'):
				continue

			line = line.strip()
			line = line.split('\t')

			match = re.match(r'chr.+:(\d+)', line[0])
			covg.append(match.group(1))
		fh.close()
		break

	covg_range = [min(covg), max(covg)]
	x_axis_label="{}: {}-{}".format(chrom, covg_range[0], covg_range[1])
	title = "{}: {}-{} {}bp".format(chrom, covg_range[0], covg_range[1], var_size)
		

	all_sample_covg = {}
	y_labels = []

	for sample in os.listdir(coverage_output): 
		match = re.match(r'(.+)\.covg_metrics_output', sample)
		person = match.group(1)


		y_labels.append(person)
		all_sample_covg[person] = OrderedDict()

		fh = open(coverage_output+'/'+sample, 'r')
		for line in fh:
			if line.startswith("Locus"):
				continue

			line = line.strip()
			line = line.split('\t')

			match = re.match(r'chr.+:(.+)', line[0])
			sent_locus = match.group(1)
			depth = line[1]
			all_sample_covg[person][int(sent_locus)] = int(depth)

		fh.close()

	y_labels.sort()

	#make dict to knoe which individual matlab callde the variant in
	has_dict = {}
	fh = open(matlab_output, 'r')
	for line in fh:
		if line.startswith('Sample'):
			continue

		line = line.strip()
		line = line.split('\t')	
		match = re.match(r'(.+)\.covg_metrics_output', line[0])
		ind = match.group(1)

		if line[1].startswith('chr'):
			has_dict[ind] = 'y'
		else:
			 has_dict[ind] = 'n'
	fh.close()
	

	x = []
	y = []

	for person in y_labels:
		x.append(list(all_sample_covg[person].keys()))
		y.append(list(all_sample_covg[person].values()))

	
	##hardcoded for now
        ##might need to make more flexible later
	size = 25
	pad = 55


	#fig, axs = plt.subplots(size, sharex=True, sharey=True, gridspec_kw={'hspace': 0}, figsize=(8,10))
	fig, axs = plt.subplots(size, sharex=True, gridspec_kw={'hspace': 0.3}, figsize=(12,20))
	for i in range(size):
		axs[i].plot(x[i], y[i], linewidth=1)

	i=0
	for ax in axs:
		ax.yaxis.set_ticks_position('right')
		ax.set_yticks([min(y[i]), max(y[i])])
		#print(min(y[i]), max(y[i]))
		y_lim = ax.get_ylim()
		#print(y_lim)
		if has_dict[y_labels[i]] == 'y':
			ax.set_ylabel(y_labels[i], labelpad=pad+5, rotation='horizontal', va='center', fontweight='bold', fontsize=12)
		else:
			ax.set_ylabel(y_labels[i], labelpad=pad, rotation='horizontal', va='center', fontsize=12)

		ax.set_xticks(var_coords)
		ax.set_xticklabels([str(x) for x in var_coords], fontsize=15)
		i += 1

	for ax in axs:
		ax.label_outer()

	plt.xticks(rotation=45, va='top', ha='right')

	axs[0].set_title(title, fontsize=20)
	axs[0].title.set_position([.5, 1.05])
	plt.xlabel(x_axis_label, fontsize=15)
	#plt.tight_layout()

	plt.savefig('{}.png'.format(outfile), format='png', dpi=500)
	plt.close('all')

#############################################


plot_coverage(args.variant_locus, args.coverage_output, args.matlab_output, args.outfile)
