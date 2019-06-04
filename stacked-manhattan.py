"""
Given list of pylmm results files, creates Manhattan plots for each, stacking
  	them on top of each other.
"""


import argparse, glob
import numpy as np
from scipy.stats import pearsonr, spearmanr
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure


def parseargs():    # handle user arguments
	parser = argparse.ArgumentParser(description="Make stacked Manhattan" +
				" plots for list of pylmm output files.")
	parser.add_argument('--pylmm', required = True, nargs = '+',
		help = 'List of pylmm output file(s), space-separated if multiple.')
	parser.add_argument('--dir', action = 'store_true',
		help = 'Use if first item in --pylmm list is a directory.')
	parser.add_argument('--all_strains',
		default = '/u/home/n/nlapier2/mousedata/all_strains.tped',
		help = 'Location of all_strains.tped. Default is known location.')
	parser.add_argument('--outname', default = 'manhattans.png',
		help = 'Where to save Manhattan plots figure.')
	parser.add_argument('--plink12', action = 'store_true',
		help = 'Use if this tped has been put into plink12 format.')
	args = parser.parse_args()
	return args


def parse_allstrains(all_strains_loc, is_plink12):
	jax2loc, snpsPerChrom, mm10order = {}, {}, {}
	with(open(all_strains_loc, 'r')) as infile:
		if not is_plink12:
			infile.readline()  # skip header
		for line in infile:
			if not is_plink12:
				splits = line[:200].split('\t')
			else:
				splits = line[:200].split(' ')
			jax, chromosome, mm10 = splits[1], splits[0], int(splits[3])
			if chromosome == 'X':
				chromosome = 20
			elif chromosome == 'Y':
				chromosome = 21
			else:
				chromosome = int(chromosome)
			jax2loc[jax] = [chromosome, mm10]
			if chromosome not in snpsPerChrom:
				snpsPerChrom[chromosome] = [mm10]
			else:
				snpsPerChrom[chromosome].append(mm10)

	for chromosome in snpsPerChrom:
		mm10order[chromosome] = {}
		sorted_mm10s = sorted(snpsPerChrom[chromosome])
		for i in range(len(sorted_mm10s)):
			mm10order[chromosome][sorted_mm10s[i]] = i
		snpsPerChrom[chromosome] = len(snpsPerChrom[chromosome])
	return jax2loc, snpsPerChrom, mm10order


def parse_pylmm(phenofile):
	"""
	Read SNP pvalues for a given trait from pylmm results.
	Argument: phenofile is the phenotype file to read
	Returns: dict mapping JAX SNP names to their pvalue for the trait
	"""
	pylmm_pvals = {}
	with(open(phenofile, 'r')) as infile:
		infile.readline()  # skip header
		for line in infile:
			splits = line.split()
			if splits[-1] == 'nan':
				continue
			jax, pval = splits[0], float(splits[-1])
			pylmm_pvals[jax] = pval
	return pylmm_pvals


def pack_pvals_and_locs(jax2loc, snpsPerChrom, mm10order, study_pvals):
	X_by_chrom = [[] for i in range(len(snpsPerChrom))]
	y_all_studies = [[] for i in range(len(snpsPerChrom))]
	# map chromosome to location (indexed mm10 scaled to 1) and list of study_pvals
	total_snps = float(sum([v for k,v in snpsPerChrom.items()]))
	chrom_proportions = {k:v/total_snps for k,v in snpsPerChrom.items()}
	chrom_start_pos = []
	for chromosome in range(len(chrom_proportions)):
		if chromosome == 0:
			chrom_start_pos.append(0.0)#chrom_proportions[1])
		else:
			chrom_start_pos.append(chrom_start_pos[-1] + chrom_proportions[chromosome+1])

	for jax in jax2loc:
		jax_pvals = []
		for study in study_pvals:
			if jax in study:
				if study[jax] == 0.0:
					study[jax] = 0.00000001
				jax_pvals.append(study[jax])
			else:
				jax_pvals.append(1.0)

		chrom, mm10 = jax2loc[jax]
		jax_x_val = chrom_start_pos[chrom-1] + chrom_proportions[chrom] * (
			float(mm10order[chrom][mm10])/float(snpsPerChrom[chrom]))
		X_by_chrom[chrom-1].append(jax_x_val)
		y_all_studies[chrom-1].append(jax_pvals)

	for chrom in range(len(y_all_studies)):
		y_all_studies[chrom] = -np.log10(np.array(y_all_studies[chrom]).T)
	return np.array(X_by_chrom), y_all_studies, chrom_start_pos


def make_manhattan_plots(studies, saveloc, X_by_chrom, y_all_studies, chrom_start_pos):
	num_plots = len(studies)
	colors = ['black', 'grey']
	x_labels = [str(i+1) for i in range(len(X_by_chrom))]

	figure(num=None, figsize=(8, len(studies)*4), dpi=80, facecolor='w', edgecolor='k')
	plt.title('Manhattan plots for pylmm results')
	plt.subplots_adjust(hspace=0.2)
	for i in range(len(studies)):
		ax = plt.subplot(len(studies), 1, i+1)
		for chrom in range(len(X_by_chrom)):
			this_color = colors[(chrom+1) % len(colors)]
			plt.scatter(X_by_chrom[chrom], y_all_studies[chrom][i], color = this_color, s=2.5)
			sig = -np.log10(4.1 * 10**(-6))
			plt.plot([0,1], [sig,sig], color='red', linewidth=0.05)  # significance line
			plt.ylabel('-log10(pvalues)')
			ax.set_title(studies[i])
			ax.set_xticks(chrom_start_pos)
			ax.set_xticklabels(x_labels)
			ax.set_xlim([0.0, 1.0])
			#ax.set_ylim([0, 7.5])
	plt.xlabel('Chromosome')
	plt.savefig(saveloc)


def main():
	args = parseargs()  # handle user arguments
	if args.dir:
		if not args.pylmm[0].endswith('/'):
			args.pylmm[0] += '/'
		dir_studies = sorted(glob.glob(args.pylmm[0]+'*'))
		args.pylmm = dir_studies + args.pylmm[1:]

	jax2loc, snpsPerChrom, mm10order = parse_allstrains(args.all_strains, args.plink12)
	study_pvals = []
	for study in args.pylmm:
		study_pvals.append(parse_pylmm(study))  # get pylmm SNPs to pvalues
	X_by_chrom, y_all_studies, chrom_start_pos = pack_pvals_and_locs(
		jax2loc, snpsPerChrom, mm10order, study_pvals)
	make_manhattan_plots(args.pylmm, args.outname, X_by_chrom, y_all_studies, chrom_start_pos)
	#compare_pvals(args, sql_pvals, pylmm_pvals)  # comparison metrics & plots


if __name__ == '__main__':
	main()
#
