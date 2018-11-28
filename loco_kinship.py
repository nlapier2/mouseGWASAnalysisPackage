import argparse, os, subprocess, sys


def parseargs():    # handle user arguments
	parser = argparse.ArgumentParser(description="Given plink files," +
				" generate LOCO tpeds/maps and kinship matrices.")
	parser.add_argument('--plink', required = True,
		help = 'Base name of plink files.')
	parser.add_argument('--outdir', default='loco_pylmm/',
		help = 'Name of output directory for LOCO and kinship files.')
	parser.add_argument('--pylmm',
		default = '/u/home/n/nlapier2/mousedata/pylmm/scripts/pylmmKinship.py',
		help = 'Executable path to pylmmKinship.py. Default is known location.')
	args = parser.parse_args()
	return args


def determine_num_chroms(tpedname):
	"""
	Figure out number of chromosomes in tped file by taking the highest seen
	  	chromosome number. Assumes no missing numbers.
	Argument: tpedname is path to plink tped file.
	Returns: number of chromosomes
	"""
	num_chroms = 0
	with(open(tpedname, 'r')) as tped:
		for line in tped:
			chrom = int(line.split()[0])
			if chrom > num_chroms:
				num_chroms = chrom
	return num_chroms


def gen_loco_files(args, num_chroms):
	"""
	Generate leave-one-chromosome-out (LOCO) files for later pylmm analysis.
	  	Essentially, these files are copies of the plink files with one
		chromosome left out each time (thus there is one copy per chromosome).
	Also generate "one-chromosome" (one_chr) files for downstream GWAS analysis,
	  	which are copies of the plink files with only one chromosome included.
		Thus, for each chromosome, the plink files are split into LOCO files
		excluding that chromosome and one_chr files with only that chromosome.
	pylmmGWAS is then run on each one_chr file taking as input the corresponding
		LOCO kinship matrix -- this is the definition of LOCO analysis.
	Arguments:
	-- args are the user specified arguments handled by argparse
	-- num_chroms is the number of chromosomes in the plink files
	Returns loco_fnames, a list of the LOCO file names
	"""
	# LOCO files for LOCO kinship analysis, one_chr files for subsequent GWAS
	loco_fnames = [args.outdir+'loco_chr'+str(i+1) for i in range(num_chroms)]
	onechr_fnames = [args.outdir+'one_chr'+str(i+1) for i in range(num_chroms)]
	tped_locos = [open(loco_fnames[i] +'.tped', 'w') for i in range(num_chroms)]
	map_locos = [open(loco_fnames[i] + '.map', 'w') for i in range(num_chroms)]
	tped_onechrs=[open(onechr_fnames[i]+'.tped','w') for i in range(num_chroms)]
	map_onechrs= [open(onechr_fnames[i]+'.map', 'w') for i in range(num_chroms)]

	#write full tped file to tped LOCO files, leaving appropriate chromosome out
	with(open(args.plink + '.tped', 'r')) as tped:
		with(open(args.plink + '.map', 'r')) as map:
			for line in tped:
				mapline = map.readline()
				chrom = int(line.split()[0])
				for i in range(len(tped_locos)):
					if i+1 != chrom:
						tped_locos[i].write(line)
						map_locos[i].write(mapline)
					else:
						tped_onechrs[i].write(line)
						map_onechrs[i].write(mapline)
	for i in range(len(tped_locos)):
		tped_locos[i].close()
		map_locos[i].close()
		tped_onechrs[i].close()
		map_onechrs[i].close()
	for fname in loco_fnames:  # tfam file doesn't change by chromosome
		subprocess.Popen(['cp', args.plink + '.tfam', fname + '.tfam']).wait()
	for fname in onechr_fnames:
		subprocess.Popen(['cp', args.plink + '.tfam', fname + '.tfam']).wait()
	return loco_fnames


def run_pylmm_kinship(args, num_chroms, loco_fnames):
	"""
	Run pylmmKinship.py on the LOCO files to generate all LOCO kinship matrices.
	Arguments:
	-- args are the user specified arguments handled by argparse
	-- num_chroms is the number of chromosomes in the plink files
	-- loco_fnames is a list of the LOCO file names
	"""
	# run pylmmKinship for each LOCO file
	for i in range(num_chroms):
		print('Running pylmm on ' + loco_fnames[i])
		subprocess.Popen(['python', args.pylmm, '--tfile', loco_fnames[i],
		 	loco_fnames[i] + '.kin']).wait()
		# clean loco files; no longer needed
		subprocess.Popen(['rm', loco_fnames[i] + '*']).wait()


def main():
	args = parseargs()
	if args.plink.endswith('.'):
		args.plink = args.plink[:-1]
	# make directories and open file handlers
	if not(args.outdir.endswith('/')):
		args.outdir += '/'
	if not os.path.isdir(args.outdir):
		os.makedirs(args.outdir)

	# determine number of chromosomes in file,
	#  	generate leave-one-chromosome-out (LOCO) files, run pylmmKinship.py
	num_chroms = determine_num_chroms(args.plink + '.tped')
	loco_fnames = gen_loco_files(args, num_chroms)
	run_pylmm_kinship(args, num_chroms, loco_fnames)


if __name__ == '__main__':
	main()
#
