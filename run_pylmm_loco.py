import os, subprocess, sys

if len(sys.argv) != 4:
	print('Usage: python run_pylmm_loco.py loco_dir/ myPhenofile.txt path/to/pylmmGWAS.py')
	print('Runs pylmmGWAS on each chromosome, then aggregates.')
	sys.exit()
loco_dir, phenofile, pylmm_path = sys.argv[1:]
if not loco_dir.endswith('/'):
	loco_dir += '/'
outname = loco_dir + 'aggregated-gwas-results.txt'
open(outname, 'w').close()  # clear file and check if writeable

# Run pylmm for all chromosomes, using loco_chr for the kinship and the
#  	non-left-out chromosome as the set of SNPs to test
num_chroms = 21  # for mice
for i in range(num_chroms):
	kin_name = loco_dir + 'loco_chr' + str(i+1) + '.kin'
	if not os.path.exists(kin_name):
		break  # subset of chromosomes available
	one_chr_name = loco_dir + 'one_chr' + str(i+1)
	subprocess.Popen(['python', pylmm_path, '--tfile', one_chr_name,
						'--kfile', kin_name, '--phenofile', phenofile,
						loco_dir + 'gwas-out-chrom' + str(i+1) + '.txt']).wait()

# Aggregate the files
with(open(outname, 'a')) as outfile:
	for i in range(num_chroms):
		gwas_outfile = loco_dir + 'gwas-out-chrom' + str(i+1) + '.txt'
		if not os.path.exists(gwas_outfile):
			break  # only subset of chromosomes available
		if i == 0:  # keep header from first file
			subprocess.Popen(['cat', gwas_outfile], stdout=outfile).wait()
		else:
			subprocess.Popen(['tail', '-n', '+2', gwas_outfile], stdout=outfile).wait()
#

