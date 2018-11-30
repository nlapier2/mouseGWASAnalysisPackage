"""
Compares p-values between SQL QTL file and pylmm results. Shows corerlations
  	and significant hits, and makes some plots.
"""


import argparse
import numpy as np
from scipy.stats import pearsonr, spearmanr
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


def parseargs():    # handle user arguments
	parser = argparse.ArgumentParser(description="Compare pvals" +
				" from pylmm with those from sql database.")
	parser.add_argument('--pylmm', required = True, help = 'pylmm pheno file.')
	parser.add_argument('--qtls', required = True,
		help = 'Clinical QTLs from SQL DB.')
	parser.add_argument('--trait_name', required = True,
		help = 'Name of trait being studied, needed for SQL QTL file.')
	parser.add_argument('--fastlmm', action='store_true',
		help = 'Use this flag if --qtls is raw fastlmm output.')
	parser.add_argument('--out_basename', default='',
		help = 'Plot output base names.')
	args = parser.parse_args()
	return args


def get_jax2rsid(args):
	"""
	Read a pre-specfied "rsid2jax.txt" file that maps SNP names from rsID to
	  	JAX (Jackson labs). Needed for comparing the same SNPs between
		different files.
	Argument: args are the user arguments parsed with argparse.
	Returns: dict mapping rsID SNP names to JAX (Jackson Labs) SNP names
	"""
	rsid2jax = {}
	with(open('/u/home/n/nlapier2/mousedata/rsid2jax.txt', 'r')) as infile:
		for line in infile:
			splits = line.strip().split()
			rsid2jax[splits[0]] = splits[1]  # splits[0]=rsID, splits[1]=JAX
	return rsid2jax


def parse_sql(args, rsid2jax):
	"""
	Read SNP pvalues for a given trait from clinical QTL file from SQL database.
	Arguments:
	-- args are the user arguments parsed with argparse.
	-- rsid2jax is the rsID to JAX SNP name mapping from the get_rsid2jax method
	Returns: dict mapping JAX SNP names to their pvalue for the trait
	"""
	sql_pvals = {}
	with(open(args.qtls, 'r')) as infile:
		infile.readline()  # skip header
		for line in infile:
			splits = line.split('\t')
			if splits[0] != args.trait_name:  # QTL for different trait; ignore
				continue
			rsid, pval = splits[2], float(splits[-3])
			#if rsid not in rsid2jax:
			#	continue
			jax = rsid2jax[rsid]  # translate the rsID to JAX SNP name
			sql_pvals[jax] = pval
	return sql_pvals


def parse_pylmm(phenofile, fastlmm=False):
	"""
	Read SNP pvalues for a given trait from pylmm results.
	Arguments:
	-- phenofile is the phenotype file to read
	-- fastlmm flag parses fastlmm format file, which is similar format
	Returns: dict mapping JAX SNP names to their pvalue for the trait
	"""
	pylmm_pvals = {}
	with(open(phenofile, 'r')) as infile:
		infile.readline()  # skip header
		for line in infile:
			splits = line.split()
			if not fastlmm:
				if splits[-1] == 'nan':
					continue
				jax, pval = splits[0], float(splits[-1])
			else:
				jax, pval = splits[2], float(splits[6])
			pylmm_pvals[jax] = pval
	return pylmm_pvals


def plot_pvals(args, sql_pvals, pylmm_pvals):
	"""
	Make plots comparing the SQL (Fast-LMM) results with the pylmm results.
	Arguments:
	-- args are the user specified arguments handled by argparse
	-- sql_pvals are the pvalues for each SNP from the SQL database
	-- pylmm_pvals are the pvalues for each SNP from the pylmm results
	"""
	# This plot sorts pvalues for each method and plots their negative logs.
	# It is useful for checking for general inflation of pvalues for one method
	#  	versus the other.
	ord_sql_pvals = -np.log(sorted(sql_pvals))
	ord_pylmm_pvals = -np.log(sorted(pylmm_pvals))
	plt.plot(ord_sql_pvals, color='blue', label='fastlmm')
	plt.plot(ord_pylmm_pvals, color='red', label='pylmm')
	plt.legend()
	plt.title('Ordered p-values for each result')
	plt.xlabel('Most to least significant')
	plt.ylabel('-log10(pval)')
	plt.savefig(args.out_basename + 'both-methods-ordered-pvals.png')
	plt.close()

	# This plot zooms in on the top 1000 most significant SNP pvalues for SQL
	#  	and also plots the pvalues for those same SNPs for pylmm. This is useful
	#  	for getting a close look at the correlation on the most signifcant SNPs.
	sql_pvals, pylmm_pvals = np.log(sql_pvals), np.log(pylmm_pvals)
	pvals_w_order = [[i, sql_pvals[i]] for i in range(len(sql_pvals))]
	pvals_w_order.sort(key=lambda x: x[1])
	pvals_w_order = pvals_w_order[:1000]
	indices = [i[0] for i in pvals_w_order]
	pylmm_pvals_order = np.array(pylmm_pvals)[indices]
	plt.plot(-np.array((sorted(sql_pvals)[:1000])), 'bo', label='fastlmm')
	plt.plot(-pylmm_pvals_order, 'ro', label='pylmm')
	plt.legend()
	plt.title('Comparison of results on SQL most signifcant SNPs')
	plt.xlabel('Most to least significant')
	plt.ylabel('-log10(pval)')
	plt.savefig(args.out_basename + 'sql-ordered-snp-pvals.png')
	plt.close()

	# Same as above, but now using the top 1000 SNPs of pylmm instead
	pvals_w_order = [[i, pylmm_pvals[i]] for i in range(len(pylmm_pvals))]
	pvals_w_order.sort(key=lambda x: x[1])
	pvals_w_order = pvals_w_order[:1000]
	indices = [i[0] for i in pvals_w_order]
	sql_pvals_order = np.array(sql_pvals)[indices]
	plt.plot(-np.array((sorted(pylmm_pvals)[:1000])), 'ro', label='pylmm')
	plt.plot(-sql_pvals_order, 'bo', label='fastlmm')
	plt.legend()
	plt.title('Comparison of results on pylmm most signifcant SNPs')
	plt.xlabel('Most to least significant')
	plt.ylabel('-log10(pval)')
	plt.savefig(args.out_basename + 'pylmm-ordered-snp-pvals.png')


def compare_pvals(args, sql_results, pylmm_results):
	"""
	Compare pvalue results between the SQL database (Fast-LMM) and pylmm.
	Prints correlation & signifcance metrics and makes some plots.
	Arguments:
	-- args are the user specified arguments handled by argparse
	-- sql_pvals are the pvalues for each SNP from the SQL database
	-- pylmm_pvals are the pvalues for each SNP from the pylmm results
	"""
	# see how many SNPs are signifcant for each method, both with a basic 0.05
	#  	pvalue and a naively bonferroni corrected pvalue of about 2*(10**-7)
	bonf_cutoff = 0.05 / len(sql_results)
	sql_pvals, pylmm_pvals = [], []  # put SNPs in same order for each method
	raw_sql_sig, raw_pylmm_sig, raw_both_sig = 0, 0, 0
	bonf_sql_sig, bonf_pylmm_sig, bonf_both_sig = 0, 0, 0
	for jax in sql_results:
		if jax not in pylmm_results:  # only look at SNPs that both methods have
			continue
		sqlpval, pylpval = sql_results[jax], pylmm_results[jax]
		if sqlpval < 0.05:
			raw_sql_sig += 1
			if sqlpval < bonf_cutoff:
				bonf_sql_sig += 1
		if pylpval < 0.05:
			raw_pylmm_sig += 1
			if pylpval < bonf_cutoff:
				bonf_pylmm_sig += 1
		if sqlpval < 0.05 and pylpval < 0.05:
			raw_both_sig += 1
			if sqlpval < bonf_cutoff and pylpval < bonf_cutoff:
				bonf_both_sig += 1
		sql_pvals.append(sqlpval)
		pylmm_pvals.append(pylpval)

	plot_pvals(args, sql_pvals, pylmm_pvals)  # make some plots

	# Compute mean squared difference and Pearson/Spearman correlations between
	#  	pvalues for the different methods.
	sql_pvals, pylmm_pvals = np.array(sql_pvals), np.array(pylmm_pvals)
	mse = np.mean([(sql_pvals[i] - pylmm_pvals[i])**2
						for i in range(len(sql_pvals))])
	pearson = pearsonr(sql_pvals, pylmm_pvals)
	spearman = spearmanr(sql_pvals, pylmm_pvals)

	# Nice printing for user.
	print('Mean Squared Error: ' + str(mse))
	print('Pearson (correlation, p-value): ' + str(pearson))
	print('Spearman rank (correlation, p-value): ' + str(spearman))
	print('Total SNPs in sqlQTL: ' + str(len(sql_pvals)))
	print('Significant at pval<0.05 [sqlQTL, pylmm, both]: ' +
	 	str([raw_sql_sig, raw_pylmm_sig, raw_both_sig]))
	print('Bonferroni corrected p-value: ' + str(bonf_cutoff))
	print('Bonferroni significant [sqlQTL, pylmm, both]: ' +
		str([bonf_sql_sig, bonf_pylmm_sig, bonf_both_sig]))


def main():
	args = parseargs()  # handle user arguments
	if not args.fastlmm:
		rsid2jax = get_jax2rsid(args)  # map SNP names from rsID to JAX
		sql_pvals = parse_sql(args, rsid2jax)  # get SQL SNPs to pvalues
	else:
		sql_pvals = parse_pylmm(args.qtls, fastlmm=True)  #if raw fastlmm format
	pylmm_pvals = parse_pylmm(args.pylmm)  # get pylmm SNPs to pvalues
	compare_pvals(args, sql_pvals, pylmm_pvals)  # comparison metrics & plots


if __name__ == '__main__':
	main()
#
