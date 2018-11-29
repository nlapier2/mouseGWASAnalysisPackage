"""
This script preprocesses a given clinical traits file, and should always be run
  	upon downloading such a file from the SQL server. This file must be
	tab-delimited (see the wiki for more information).
The main idea is that phenotype and strain names will be standardized between
	all studies. This is done by reading in a manually curated mapping file
	(which we call pheno_map here) that maps names seen in the clinical trait
	files to the names we want to have.
Some other minor tweaks are made. The output is a new clinical_traits file with
  	all these tweaks and substitutions.
This script is not considered "analysis" and is thus not included in the
  	submit-pylmm-analysis.sh or simple-analysis.py wrapper scripts.
"""


import argparse, glob, subprocess, sys


def parseargs():    # handle user arguments
	parser = argparse.ArgumentParser(description="Prepare SQL clincal trait" +
				" and QTL files for further processing and pylmm analysis.")
	parser.add_argument('--clinical', required = True,
		help = 'Clinitical trait tsv file.')
	parser.add_argument('--qtls', default = 'NONE',
		help = 'Clinitical QTLs tsv file.')  # DEPRECATED
	parser.add_argument('--outname', default='AUTO', help = 'Output file name.')
	parser.add_argument('--pheno_map',
		default = '/u/home/n/nlapier2/mousedata/pheno_map.txt',
		help = 'Specify phenotype map file to read. Default to known location.')
	args = parser.parse_args()
	return args


def read_pheno_map(mapfile):
	"""
	Read in a file mapping original phenotype names to new names.
	This allows us to keep phenotype names consistent across traits.
	Argument: the args.pheno_map file specified by the user (or default "AUTO")
	Returns: dict mapping old phenotype names to new phenotype names
	"""
	pheno_map = {}
	# Parse mapfile contents into a dict. The file contains
	#  	tab separated original/new pairs, so splits[0] is the original and
	#  	splits[1] is the new name. Dict thus maps original name to new name.
	with(open(mapfile, 'r')) as f:
		for line in f:
			if line.startswith('#'):
				continue
			splits = line.strip().split('\t')
			pheno_map[splits[0]] = splits[1]
	return pheno_map


def check_trait_QTL_arg_order(args):
	"""
	A simple method that checks whether the args.clinical and args.qtls user
	  	arguments are in the correct order.
	Argument: args are the user arguments parsed from argparse.
	Returns: the strings that will be the final args.clinical and args.qtls
	"""
	with(open(args.clinical, 'r')) as infile:
		# check if user accidentally switched QTL and traits arguments
		if 'pvalue' in infile.readline():  # pvalue field only in QTLs file
			return args.qtls, args.clinical  # switched order
		else:
			return args.clinical, args.qtls  # don't switch


def preprocess_traits(args, pheno_map):
	"""
	Preprocess the clinical traits file for easier downstream usage.
	Mainly, replace the pheno_map strings and 'NULL' with 'NA' (for R).
	Arguments:
	-- args: the user arguments parsed by parseargs
	-- pheno_map: the phenotype map read in by the read_pheno_map method.
	"""
	with(open(args.clinical, 'r')) as traitfile:
		header = traitfile.readline().strip().split('\t')
		# remove leading non-ascii characters, replace col names using pheno_map
		header[0] = ''.join([ch for ch in header[0] if ord(ch) < 128])
		header = [i if i not in pheno_map else pheno_map[i] for i in header]
		# Determine columns that define mouse number and strain. If default
		#  	phenotype map file is used, columns guaranteed to have these names.
		mousecol, straincol=header.index('mouse_number'), header.index('Strain')

		with(open(args.outname + '_clinical_traits.tsv', 'w')) as outfile:
			outfile.write('\t'.join(header) + '\n')
			for line in traitfile:
				splits = line.strip().split('\t')
				if len(splits) < 2:  # trailing line at end of file
					break
				#splits[straincol] = splits[straincol].replace('/', '.')
				# replace strain name using pheno_map, if necessary
				if splits[straincol] in pheno_map:
					splits[straincol] = pheno_map[splits[straincol]]
				for i in range(len(splits)):  # replace NULL with NA, for R
					if splits[i] == 'NULL':
						splits[i] = 'NA'
				outfile.write('\t'.join(splits) + '\n')


def preprocess_qtls(args):
	"""
	DEPRECATED. Do not use this method unless you know what you're doing.
	The goal of this method is to take the massive QTL files and output them in
	  	a more compact form, avoiding the repetition of trait and SNP info.
	Reformats QTLs file in a more compact transposed form with one
	  	line per trait. First two columns of each line are the trait name and
		category. Subsequent columns are SNPs, with each entry corresponding to
		the pvalue,snp_weight,snp_odds_ratio for that SNP for the line's trait.
	Also outputs a snp2info file that maps SNP rsID to chromosome number and
	  	other info (location, LD block, etc.).
	Argument: args are the user arguments parsed from argparse.
	"""
	snpcols, allsnps = {}, []  # defines column for each SNP
	trait2cat, trait_pvals = {}, {}  # maps traits to their category & pvalues
	with(open(args.qtls, 'r')) as qtlfile:
		with(open('preprocessed_snp2info.txt', 'w')) as snpfile:
			with(open('preprocessed_clinical_QTL.tsv', 'w')) as outfile:
				header = qtlfile.readline().split('\t')

				# SNP file contains the third through ninth columns in the
				#  	original QTL file, which contain SNP-related info that is
				#  	unrelated to specific traits and can thus be separated out.
				snpfile.write('\t'.join(header[2:9]) + '\n')  # snps file header
				for line in qtlfile:
					splits = line.strip().split('\t')
					if len(splits) < 2:
						break
					trait, category, snp = splits[0:3]
					if snp not in snpcols:
						snpcols[snp] = len(snpcols)  # assign column to this SNP
						allsnps.append(snp)  # maintains SNP order
						snpfile.write('\t'.join(splits[2:9]) + '\n')

					# Enter NA for all trait/SNP pairs missing p-values,
					#  	and ensure that all traits have all columns specified.
					if trait not in trait_pvals:
						trait_pvals[trait] = ['NA' for i in range(len(snpcols))]
						trait2cat[trait] = category
					elif len(trait_pvals[trait]) <= snpcols[snp]:
						trait_pvals[trait].extend(['NA' for i in
							range(snpcols[snp] + 1 - len(trait_pvals[trait]))])
					trait_pvals[trait][snpcols[snp]] = ','.join(splits[-3:])

				# write out the traits with their category and p-values
				outfile.write('trait_name\ttrait_category\t' +
					'\t'.join(allsnps) + '\n')
				for trait in trait_pvals:
					# Ensure that all traits have all SNP columns specified.
					if len(trait_pvals[trait]) < len(snpcols):
						trait_pvals[trait].extend(['NA' for i in
							range(len(snpcols) - len(trait_pvals[trait]))])
					outfile.write(trait + '\t' + trait2cat[trait] + '\t' +
						'\t'.join(trait_pvals[trait]) + '\n')


def main():
	args = parseargs()
	if args.outname == 'AUTO':
		args.outname = 'preprocessed_clinical_traits.tsv'
	# ensure args.clinical and args.qtls were specified in correct order
	args.clinical, args.qtls = check_trait_QTL_arg_order(args)

	# Reads in a phenotype mapping file, uses this to replace phenotype names
	#  	in the clincal traits file.
	pheno_map = read_pheno_map(args.pheno_map)
	preprocess_traits(args, pheno_map)

	# DEPRECATED: performs QTL file reformatting into sparser form.
	if args.qtls != 'NONE':
		print('Warning: --qtls is deprecated; ' +
			'only use if you know what you are doing.')
		preprocess_qtls(args)


if __name__ == '__main__':
	main()
#
