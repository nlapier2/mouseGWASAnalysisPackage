import argparse
import numpy as np
from scipy.stats import pearsonr, spearmanr
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


def parseargs():    # handle user arguments
	parser = argparse.ArgumentParser(description="Compare pvals," +
				" from pylmm with those from sql database.")
	parser.add_argument('rsid2jax', help = 'rsID to jax SNP name file.')
	parser.add_argument('pylmm', help = 'pylmm pheno file.')
	parser.add_argument('qtls', help = 'Clinical QTLs from SQL DB.')
	args = parser.parse_args()
	return args


def get_jax2rsid(args):
	rsid2jax = {}
	with(open(args.rsid2jax, 'r')) as infile:
		for line in infile:
			splits = line.strip().split()
			rsid2jax[splits[0]] = splits[1]
	return rsid2jax


def parse_sql(args, rsid2jax):
	sql_phenos = {}
	with(open(args.qtls, 'r')) as infile:
		infile.readline() ; infile.readline()  # skip headers
		for line in infile:
			splits = line.split()
			rsid, pval = splits[2], float(splits[-3])
			jax = rsid2jax[rsid]
			sql_phenos[jax] = pval
	return sql_phenos


def parse_pylmm(args):
	pylmm_phenos = {}
	with(open(args.pylmm, 'r')) as infile:
		infile.readline()  # skip headers
		for line in infile:
			splits = line.split()
			jax, pval = splits[0], float(splits[-1])
			pylmm_phenos[jax] = pval
	return pylmm_phenos


def plot_pvals(sql_pvals, pylmm_pvals):
	ord_sql_pvals = -np.log(sorted(sql_pvals))
	ord_pylmm_pvals = -np.log(sorted(pylmm_pvals))
	plt.plot(ord_sql_pvals, color='blue', label='sqlQTL')
	plt.plot(ord_pylmm_pvals, color='red', label='pylmm')
	plt.legend()
	plt.title('Ordered p-values for each result')
	plt.xlabel('Most to least significant')
	plt.ylabel('-log10(pval)')
	plt.savefig('ordered-pvals.png')
	plt.close()

	sql_pvals, pylmm_pvals = np.log(sql_pvals), np.log(pylmm_pvals)
	pvals_w_order = [[i, sql_pvals[i]] for i in range(len(sql_pvals))]
	pvals_w_order.sort(key=lambda x: x[1])
	pvals_w_order = pvals_w_order[:1000]
	indices = [i[0] for i in pvals_w_order]
	pylmm_pvals_order = np.array(pylmm_pvals)[indices]

	#plt.plot(-np.log(sql_pvals[:30]), 'bo', label='sqlQTL')
	#plt.plot(-np.log(pylmm_pvals[:30]), 'ro', label='pylmm')
	plt.plot(-np.array((sorted(sql_pvals)[:1000])), 'bo', label='sqlQTL')
	plt.plot(-pylmm_pvals_order, 'ro', label='pylmm')
	plt.legend()
	plt.title('Comparison of results on most signifcant SNPs')
	plt.xlabel('Most to least significant')
	plt.ylabel('-log10(pval)')
	plt.savefig('ordered-snp-pvals.png')


def compare_phenos(sql_phenos, pylmm_phenos):
	# get pvals in sql (pylmm has extra) and ensure same order
	# also check significance with raw (<0.05) and bonferroni corrected
	bonf_cutoff = 0.05 / len(sql_phenos)
	sql_pvals, pylmm_pvals = [], []
	raw_sql_sig, raw_pylmm_sig, raw_both_sig = 0, 0, 0
	bonf_sql_sig, bonf_pylmm_sig, bonf_both_sig = 0, 0, 0
	for jax in sql_phenos:
		sqlpval, pylpval = sql_phenos[jax], pylmm_phenos[jax]
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

	# now compare
	plot_pvals(sql_pvals, pylmm_pvals)
	sql_pvals, pylmm_pvals = np.array(sql_pvals), np.array(pylmm_pvals)
	mse = np.mean([(sql_pvals[i] - pylmm_pvals[i])**2
						for i in range(len(sql_pvals))])
	pearson = pearsonr(sql_pvals, pylmm_pvals)
	spearman = spearmanr(sql_pvals, pylmm_pvals)
	print('Mean Squared Error: ' + str(mse))
	print('Pearson (correlation, p-value): ' + str(pearson))
	print('Spearman rank (correlation, p-value): ' + str(spearman))
	print('Total SNPs in sqlQTL: ' + str(len(sql_phenos)))
	print('Significant at pval<0.05 [sqlQTL, pylmm, both]: ' +
	 	str([raw_sql_sig, raw_pylmm_sig, raw_both_sig]))
	print('Bonferroni corrected p-value: ' + str(bonf_cutoff))
	print('Bonferroni significant [sqlQTL, pylmm, both]: ' +
		str([bonf_sql_sig, bonf_pylmm_sig, bonf_both_sig]))


def main():
	args = parseargs()
	rsid2jax = get_jax2rsid(args)
	sql_phenos = parse_sql(args, rsid2jax)
	pylmm_phenos = parse_pylmm(args)
	compare_phenos(sql_phenos, pylmm_phenos)


if __name__ == '__main__':
	main()
#

