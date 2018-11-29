"""
This script takes the mice from a specified clinical traits file and matches
  	them to their genotype information using a file called all_strains.tped.
This information is then output in plink format, and then plink is used to
  	recode alleles into the numbers 1 and 2, for GWAS analysis purposes.
This is the first step of analysis, and importantly involves the decision of
  	whether to include sex chromosomes or not.
"""


import argparse, glob, subprocess, sys


def parseargs():    # handle user arguments
	parser = argparse.ArgumentParser(description="Given phenotype file," +
				" get into proper format for pylmm.")
	parser.add_argument('clinical', help = 'Clinical trait tsv file.')
	parser.add_argument('--all_strains',
		default = '/u/home/n/nlapier2/mousedata/all_strains.tped',
		help= 'Location of all_strains.tped file. Default is a known location.')
	parser.add_argument('--no_sex_chromosomes', action='store_true',
		help = 'Ignore sex chromosomes.')
	parser.add_argument('--plink_basename', default='AUTO',
		help = 'Base name of plink output.')
	args = parser.parse_args()
	return args


def get_fid_iid(args):
	"""
	Extract fid (family ID == mouse Strain) & iid (individual ID = mouse_number)
	  	from clinical traits file.
	Argument: args are the user arguments parsed from argparse.
	Returns: fid_iid, the strain / mouse number pairs from clinical traits file.
	"""
	# grab all strains in all_strains.tped; we ignore strains in the clinical
	#  	traits file that aren't in all_strains.tped because we don't have
	#  	genotype info for them
	with(open(args.all_strains, 'r')) as genofile:
		strains = genofile.readline().strip().split('\t')[4:]

	fid_iid = []  # fid == mouse strain, iid == individual mouse number
	with(open(args.clinical, 'r')) as infile:
		header = infile.readline().strip().split('\t')
		# determine columns with fid (strain) and iid (mouse number) info
		straincol, numcol = header.index('Strain'), header.index('mouse_number')
		for line in infile:
			splits = line.strip().split('\t')
			fid, iid = splits[straincol], splits[numcol]
			if fid not in strains:  # see above note on all_strains.tped
				print('WARNING: strain '+ str(fid) +' not in all_strains.tped.')
				continue
			fid_iid.append([fid, iid])
	return fid_iid


def write_tfam(fid_iid):
	"""
	Given the mouse strains (fid) and individual numbers (iid), write out a
	  	plink-format tfam file. In this case we don't specify the other tfam
		fields, so this is mostly just writing out the mice in the study.
	For more information on tfam files, see plink documentation online.
	Argument: fid_iid specifies mouse strains and numbers in this study.
	"""
	other_cols = '\t0\t0\t0\t-9\n'  # father/mother/sex/pheno unspecified
	with(open('TEMP.tfam', 'w')) as outfile:
		for pair in fid_iid:  # pair is [fid, iid] of a mouse
			outfile.write(pair[0] + '\t' + pair[1] + other_cols)


def write_tped(args, fid_iid):
	"""
	Writes plink-formatted tped file, which maps each SNP to genotype info for
	  	the mice in this study at that SNP.
	For more information on tped files, see plink documentation online.
	Arguments:
	-- args are the user arguments parsed from argparse.
	-- fid_iid contains the mouse strains (fids) & numbers (iid) in this study.
	"""
	with(open(args.all_strains, 'r')) as genofile:
		with(open('TEMP.tped', 'w')) as outfile:
			# all_strains.tped has the genotype for each SNP for each strain.
			# Here we take the mice from the clincal traits file and write
			#  	out their genotype at the SNP by indexing into the proper
			#  	header column & grabbing that strain's genotype in each line.
			# Then we write the SNP info and the genotypes for the study's mice.
			header = genofile.readline().strip().split('\t')
			genocols = [header.index(pair[0]) for pair in fid_iid]
			for line in genofile:
				splits = line.strip().split('\t')
				if args.no_sex_chromosomes and splits[0] in ['X', 'Y']:
					continue
				genos = [splits[col] for col in genocols]
				outfile.write('\t'.join(splits[:4]) + '\t' +
					'\t'.join(genos) + '\n')


def get_genotypes(args):
	"""
	Grabs genotype info for mice in a study & puts it in plink format.
	Argument: args are the user arguments parsed from argparse.
	"""
	# extract family and mouse IDs from clinical trait file, write tfam & tped
	fid_iid = get_fid_iid(args)
	write_tfam(fid_iid)
	write_tped(args, fid_iid)


def main():
	args = parseargs()
	if args.plink_basename == 'AUTO':
		# default basename is 'plink12' in same directory as args.clinical
		clinical_rpt_dir = ''
		splits = args.clinical.split('/')
		if len(splits) > 1:  # clinical report not in current dir
			clinical_dir = '/'.join(splits[:-1]) + '/'
		args.plink_basename = clinical_dir + 'plink12'

	get_genotypes(args)  # get genotypes for mice, format in plink format

	# load the plink module
	# loader = '. /u/local/Modules/default/init/modules.sh && module load plink'
	# subprocess.Popen(['bash', '-c', loader]).wait()

	# now recode (reformat) the tfam & tped files using plink
	subprocess.Popen(['plink', '--mouse', '--recode', '12', '--transpose',
						'--tfile', 'TEMP', '--out', args.plink_basename]).wait()
	subprocess.Popen(['rm', 'TEMP.tfam', 'TEMP.tped']).wait()


if __name__ == '__main__':
	main()
#
