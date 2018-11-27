import argparse, glob, subprocess, sys


def parseargs():    # handle user arguments
	parser = argparse.ArgumentParser(description="Prepare SQL clincal trait" +
				" and QTL files for further processing and pylmm analysis.")
	parser.add_argument('traits', help = 'Clinitical trait tsv file.')
	parser.add_argument('--qtls', default = 'NONE',
		help = 'Clinitical QTLs tsv file.')
	parser.add_argument('--outname', default='AUTO', help = 'Output file name.')
	parser.add_argument('--pheno_map', default='NONE',
		help = 'Specify phenotype mapping file to read.')
	args = parser.parse_args()
	return args


def read_pheno_map(mapfile):
	pheno_map = {}
	if mapfile == 'NONE':
		return pheno_map
	with(open(mapfile, 'r')) as f:
		for line in f:
			if line.startswith('#'):
				continue
			splits = line.strip().split('\t')
			pheno_map[splits[0]] = splits[1]
	return pheno_map


def check_trait_QTL_arg_order(args):
	with(open(args.traits, 'r')) as infile:
		# check if user accidentally switched QTL and traits arguments
		if 'pvalue' in infile.readline():  # pvalue field only in QTLs file
			return args.qtls, args.traits
		else:
			return args.traits, args.qtls


def preprocess_traits(args, pheno_map):
	with(open(args.traits, 'r')) as traitfile:
		header = traitfile.readline().strip().split('\t')
		# remove leading non-ascii characters, replace col names using pheno_map
		header[0] = ''.join([ch for ch in header[0] if ord(ch) < 128])
		header = [i if i not in pheno_map else pheno_map[i] for i in header]
		mousecol, straincol=header.index('mouse_number'), header.index('Strain')
		with(open(args.outname + '_clinical_traits.tsv', 'w')) as outfile:
			outfile.write('\t'.join(header) + '\n')
			for line in traitfile:
				splits = line.strip().split('\t')
				if len(splits) < 2:
					break
				#splits[straincol] = splits[straincol].replace('/', '.')
				for i in range(len(splits)):
					if splits[i] == 'NULL':
						splits[i] = 'NA'
				outfile.write('\t'.join(splits) + '\n')


def preprocess_qtls(args, pheno_map):
	snpcols, allsnps = {}, []  # defines column for each SNP
	trait2cat, trait_pvals = {}, {}  # maps traits to their category & pvalues
	with(open(args.qtls, 'r')) as qtlfile:
		with(open('preprocessed_snp2info.txt', 'w')) as snpfile:
			with(open('preprocessed_clinical_QTL.tsv', 'w')) as outfile:
				header = qtlfile.readline().split('\t')
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

					if trait not in trait_pvals:
						trait_pvals[trait] = ['NA' for i in range(len(snpcols))]
						trait2cat[trait] = category
					elif len(trait_pvals[trait]) <= snpcols[snp]:
						trait_pvals[trait].extend(['NA' for i in
							range(snpcols[snp] + 1 - len(trait_pvals[trait]))])
					trait_pvals[trait][snpcols[snp]] = ','.join(splits[-3:])

				outfile.write('trait_name\ttrait_category\t' +
					'\t'.join(allsnps) + '\n')
				for trait in trait_pvals:
					if len(trait_pvals[trait]) < len(snpcols):
						trait_pvals[trait].extend(['NA' for i in
							range(len(snpcols) - len(trait_pvals[trait]))])
					outfile.write(trait + '\t' + trait2cat[trait] + '\t' +
						'\t'.join(trait_pvals[trait]) + '\n')


def main():
	args = parseargs()
	if args.outname == 'AUTO':
		args.outname = 'preprocessed_clinical_traits.tsv'
	args.traits, args.qtls = check_trait_QTL_arg_order(args)

	pheno_map = read_pheno_map(args.pheno_map)
	preprocess_traits(args, pheno_map)
	if args.qtls != 'NONE':
		preprocess_qtls(args, pheno_map)


if __name__ == '__main__':
	main()
#
