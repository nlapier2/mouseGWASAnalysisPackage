import argparse, glob, subprocess, sys


def parseargs():    # handle user arguments
	parser = argparse.ArgumentParser(description="Given phenotype file," +
				" get into proper format for pylmm.")
	parser.add_argument('clinical', help = 'Clinitical trait tsv file.')
	parser.add_argument('trait_name', help = 'EXACT name of trait to study.')
	parser.add_argument('--clean_intermediates', action='store_true',
		help = 'Remove intermediate files.')
	parser.add_argument('--get_genotypes', default='get_genotypes_from_tped.R',
		help = 'Path to get_genotypes_from_tped.R file. Default: current dir.')
	parser.add_argument('--intermediate_basename', default='AUTO',
		help = 'Base name of intermediate files.')
	parser.add_argument('--no_sex_chromosomes', action='store_true',
		help = 'Ignore sex chromosomes.')
	parser.add_argument('--plink_basename', default='AUTO',
		help = 'Base name of plink output.')
	parser.add_argument('--study', default='AUTO',
		help = 'Specify study type. Default: read from clinical_rpt name.')
	args = parser.parse_args()
	return args


# write file mapping family id (fid) to individual id (iid), used by
#  	get_genotypes_from_tped.R. Removes mice with NULL value for trait.
def make_fid_iid(args):
	with(open(args.clinical, 'r')) as infile:
	    with(open(args.intermediate_basename, 'w')) as outfile:
			# parse header line to find column of intended trait
			splits = infile.readline().strip().split('\t')
			if args.trait_name not in splits:
				print('Error: specified trait name not seen in clinical trait file header.')
				sys.exit()
			trait_col = splits.index(args.trait_name)

			# write fid to iid file
			for line in infile:
				splits = line.strip().split('\t')
				if len(splits) < 2:
					break
				if  splits[trait_col] == 'NULL':
					continue
				if args.study == 'iron':
					fid, iid = splits[3], splits[0]
				elif args.study == 'chow':
					fid, iid = splits[1], splits[0]
				else:
					fid, iid = splits[0], splits[2]
				if fid == 'NZB/BlNJ' and iid == '1474':
					continue  # this specific mouse in one study has inconsistent values
				outfile.write(fid + '\t' + iid + '\n')


# filter out sex chromosomes
def filter_sex_chromosomes(args, fname):
	with(open(fname, 'r')) as infile:
		with(open('TRAIT_SETUP_TEMP', 'w')) as outfile:
			for line in infile:
				if line.startswith('X') or line.startswith('Y'):
					continue
				outfile.write(line)
	subprocess.Popen(['rm', fname]).wait()
	subprocess.Popen(['mv', 'TRAIT_SETUP_TEMP', fname]).wait()


# run get_genotypes_from_tped.R and plink; if user wants, cleans intermediate files
def run_get_genotypes_and_plink(args):
	subprocess.Popen(['Rscript', args.get_genotypes, args.intermediate_basename]).wait()
	if args.no_sex_chromosomes:
		filter_sex_chromosomes(args, args.intermediate_basename+'.tped')
		#filter_sex_chromosomes(args, args.intermediate_basename+'.map')
	subprocess.Popen(['plink', '--mouse', '--recode', '12', '--transpose',
						'--tfile', args.intermediate_basename, '--out', args.plink_basename]).wait()
	if args.clean_intermediates:
		for fname in glob.glob(args.intermediate_basename+'*'):
			subprocess.Popen(['rm', fname]).wait()


def main():
	args = parseargs()
	clinical_rpt_dir = ''
	splits = args.clinical.split('/')
	if len(splits) > 1:  # clinical report not in current dir
		clinical_rpt_dir = '/'.join(splits[:-1]) + '/'
	if args.intermediate_basename == 'AUTO':
		# not set by user -> auto named and placed in same dir as clinical rpt
		args.intermediate_basename = clinical_rpt_dir + 'fid_iid_' + args.trait_name.replace(' ', '_')
	if args.plink_basename == 'AUTO':
		args.plink_basename = clinical_rpt_dir + 'plink12_' + args.trait_name.replace(' ', '_')
	if args.study == 'AUTO':
		if 'iron' in args.clinical:
			args.study = 'iron'
		elif 'chow' in args.clinical:
			args.study = 'chow'

	make_fid_iid(args)
	run_get_genotypes_and_plink(args)


if __name__ == '__main__':
	main()
#