import argparse, subprocess, sys


def parseargs():    # handle user arguments
	parser = argparse.ArgumentParser(description="Given phenotype file," +
				" get into proper format for pylmm.")
	parser.add_argument('clinical_rpt', help = 'Clinitical trait rpt file.')
	parser.add_argument('trait_name', help = 'EXACT name of trait to study.')
	parser.add_argument('--clean_intermediates', action='store_true',
		help = 'Remove intermediate files.')
	parser.add_argument('--intermediate_basename', default='AUTO',
		help = 'Base name of intermediate files.')
	parser.add_argument('--plink_basename', default='AUTO',
		help = 'Base name of plink output.')
	parser.add_argument('--study', default='AUTO',
		help = 'Specify study type. Default: read from clinical_rpt name.')
	parser.add_argument('--get_genotypes', default='get_genotypes_from_tped.R',
		help = 'Path to get_genotypes_from_tped.R file. Default: current dir.')
	args = parser.parse_args()
	return args


# write file mapping family id (fid) to individual id (iid), used by
#  	get_genotypes_from_tped.R. Removes mice with NULL value for trait.
def make_fid_iid(args):
	with(open(args.clinical_rpt, 'r')) as infile:
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


# run get_genotypes_from_tped.R and plink; if user wants, cleans intermediate files
def run_get_genotypes_and_plink(args):
	subprocess.Popen(['Rscript', args.get_genotypes, args.intermediate_basename]).wait()
	subprocess.Popen(['plink', '--mouse', '--recode', '12', '--transpose',
						'--tfile', args.intermediate_basename, '--out', args.plink_basename]).wait()
	if args.clean_intermediates:
		subprocess.Popen(['rm', args.intermediate_basename+'*']).wait()


def main():
	args = parseargs()
	clinical_rpt_dir = ''
	splits = args.clinical_rpt.split('/')
	if len(splits) > 1:  # clinical report not in current dir
		clinical_rpt_dir = '/'.join(splits[:-1]) + '/'
	if args.intermediate_basename == 'AUTO':
		# not set by user -> auto named and placed in same dir as clinical rpt
		args.intermediate_basename = clinical_rpt_dir + 'fid_iid_' + args.trait_name.replace(' ', '_')
	if args.plink_basename == 'AUTO':
		args.plink_basename = clinical_rpt_dir + 'plink12_' + args.trait_name.replace(' ', '_')
	if args.study == 'AUTO':
		if 'iron' in args.clinical_rpt:
			args.study == 'iron'
		elif 'chow' in args.clinical_rpt:
			args.study == 'chow'

	make_fid_iid(args)
	run_get_genotypes_and_plink(args)


if __name__ == '__main__':
	main()
#
