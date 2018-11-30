"""
This script makes a phenotype file in pylmm format, which pylmm then uses to
  	read in the phenotype values for each mouse for GWAS analysis. Essentially,
	we are simply extracting relevant information from the clinical traits file.
One of the plink files called a tfam is required input, as it specifies which
  	mice are present and genotyped for this study.
This script also normalizes phenotypes by subtracting their mean and dividing by
  	their standard deviation. This can be turned off, but generally shouldn't
	be. Eventually this script will allow the user to specify covariates to
	regress out, but we currently do not have that option implemented.
"""


import argparse, math, sys

def parseargs():    # handle user arguments
	parser = argparse.ArgumentParser(description="Given phenotype file," +
				" get into proper format for pylmm.")
	parser.add_argument('--clinical', required = True,
		help = 'Name of existing pheno file.')
	parser.add_argument('--target', required=True,
		help = 'EXACT name of target trait to study.')
	parser.add_argument('--tfam', required=True,
		help = 'Plink tfam file; exclude mice not in it.')
	parser.add_argument('--output', default='pheno_file_pylmm.txt',
		help = 'Name of output file.')
	parser.add_argument('--no_normalization', action='store_true',
		help = 'Do not normalize target phenotype.')
	#parser.add_argument('--regress', action='store_true',
	#	help = 'Use this to regress target phenotype on other phenotypes.')
	args = parser.parse_args()
	return args


def parse_tfam(args):
	"""
	Reads the mice for this study from a tfam file. Used to ensure that mice
	  	that were excluded from analysis are also excluded from the pheno file.
	Argument: args are the user arguments parsed with argparse.
	Returns: tfams dict that maps strain to mouse number
	"""
	tfams = {}
	if args.tfam != 'NONE':
		with(open(args.tfam, 'r')) as tfamfile:
			for line in tfamfile:
				splits = line.split(' ')  # splits[0, 1] == [strain, number]
				if splits[0] not in tfams:
					tfams[splits[0]] = [splits[1]]
				else:
					tfams[splits[0]].append(splits[1])
	return tfams


def parse_clinical_file(args, tfams):
	"""
	Given a target phenotype, read the values for that phenotype.
	Arguments:
	-- args are the user arguments parsed with argparse
	-- tfams are the strain to mouse number mappings from the tfam file
	Returns:
	-- mouse2pheno: maps mouse to the value of the target phenotype and the
	  	the non-target phenos (other), which can be corrected for
	-- target_pheno: just the raw target phenotype values
	-- other_phenos: just the raw non-target (other) phenotype values
	"""
	mouse2pheno, target_pheno, other_phenos = [], [], []
	with(open(args.clinical, 'r')) as infile:
		header = infile.readline().strip().split('\t')
		targetcol = header.index(args.target)  # find column of target phenotype
		for line in infile:
			splits = line.split('\t')
			if len(splits) < 2:
				break
			if splits[1] not in tfams or splits[0] not in tfams[splits[1]]:
				continue  # this mouse not present in the tfam file
			mouse_fid_iid = splits[1] + ' ' + splits[0]  # family & indiv. ID
			# grab target phenotype and other phenotyes in this line
			tgt = splits[targetcol]
			if tgt != 'NA':
				tgt = float(tgt)
			others = splits[2 : targetcol] + splits[targetcol + 1 : ]
			# pack info together with mouse and also store target/other phenos
			mouse2pheno.append([mouse_fid_iid, tgt, others])
			if tgt != 'NA':
				target_pheno.append(tgt)
			other_phenos.append(others)
	return mouse2pheno, target_pheno, other_phenos


# NOT YET IMPLEMENTED
#def regress_pheno(args, target_pheno, other_phenos):
#	return target_pheno


def normalize_and_write(args, mouse2pheno, target_pheno):
	"""
	Normalize target phenotype values by subtracting mean and dividing by
	  	standard deviation, then write reults in pylmm format to args.output
	Arguments:
	-- mouse2pheno: maps mouse to the value of the target phenotype and the
	  	the non-target phenos (other), which can be corrected for
	-- target_pheno: just the raw target phenotype values
	-- other_phenos: just the raw non-target (other) phenotype values
	"""
	if not args.no_normalization:
		# calculate mean and standard deviation of target phenotype values
		mean = sum(target_pheno) / float(len(target_pheno))
		sq_diffs = [(i - mean)**2 for i in target_pheno]
		sample_sd = math.sqrt(sum(sq_diffs) / float(len(target_pheno) - 1))
		# normalize target phenotypes, stored in mouse2pheno[i][1]
		for i in range(len(mouse2pheno)):
			if mouse2pheno[i][1] != 'NA':
				mouse2pheno[i][1] = ((mouse2pheno[i][1] - mean) / sample_sd)

	with(open(args.output, 'w')) as outfile:
		for m in mouse2pheno:
			# m[0] is mouse strain and number, m[1] is target phenotype value
			outfile.write(m[0] + ' ' + str(m[1]) + '\n')


def main():
	args = parseargs()  # handle user arguments
	tfams = parse_tfam(args)  # read mice in this study from tfam file
	# read the phenotype values from the clinical trait file
	mouse2pheno, target_pheno, other_phenos = parse_clinical_file(args, tfams)
	#if args.regress:
	#	target_pheno = regress_pheno(args, target_pheno, other_phenos)
	# normalize phenotype values and write to args.output
	normalize_and_write(args, mouse2pheno, target_pheno)


if __name__ == '__main__':
	main()
#
