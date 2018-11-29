"""
This script is made to be a very simple wrapper to run the GWAS analysis steps
  	in get_genotypes_and_plink.py, make_pheno_file.py, loco_kinship.py, and
	run_pylmm_loco.py.
This script may be useful for users who are uncomfortable with editing files
  	from the command line, which is necessary for running the
	submit-pylmm-analysis.sh script. However, users who are comfortable editing
	files should use the shell script for more flexibility.
"""


import argparse, subprocess


# Here we define some set blocks of the shell script to run that are not
#  	dependent on user input.

SET_BLOCK_1 = "#!/bin/sh\n"+
"#$ -S /bin/bash\n"+
"#$ -N pylmm_analysis_qsub_job\n"+
"#$ -cwd\n"+
"#$ -o stdout-pylmm-analysis.out\n"+
"#$ -l h_data=10G,h_rt=24:00:00\n"+
". /u/local/Modules/default/init/modules.sh\n"+
"module load R\n"+
"module load plink\n"

SET_BLOCK_2 = "all_strains='/u/home/n/nlapier2/mousedata/all_strains.tped'\n"+
"plink_basename=$output_directory+'plink12'\n"+
"pheno_file_name=$output_directory+'pheno_file_pylmm.txt'\n"+
"no_normalization=1\n"+
"loco_outdir=$output_directory+'loco/'\n"+
"pylmm_loc='/u/home/n/nlapier2/mousedata/pylmm/'\n"+
"gwas_outfile=$output_directory+'pylmm_gwas_results.txt'\n"+
"if [ $include_sex_chromosomes = 1 ]\n"+
"then\n"+
"\tsex_chrom_opt='--no_sex_chromosomes'\n"+
"else\n"+
"\tsex_chrom_opt=''\n"
"fi\n"+
"if [ $no_normalization = 1 ]\n"+
"then\n"+
"\tnormalize_opt='--no_normalization'\n"+
"else\n"+
"\tnormalize_opt=''\n"+
"fi\n"+
"tfam=$plink_basename+'.tfam'\n"+
"pylmmKinship=$pylmm_loc+'pylmmKinship.py'\n"+
"pylmmGWAS=$pylmm_loc+'pylmmGWAS.py'\n"+
"python get_genotypes_and_plink.py $clinical_traits --all_strains $all_strains \\\n"+
"\t$sex_chrom_opt --plink_basename $plink_basename\n"+
"python make_pheno_file.py --clinical $clinical_traits --target $trait_name \\\n"+
"\t--tfam $tfam --output $pheno_file_name $normalize_opt\n"+
"python loco_kinship.py --plink $plink_basename --outdir $loco_outdir \\\n"+
"\t--pylmm $pylmmKinship\n"+
"python run_pylmm_loco.py --loco_dir $loco_outdir --pheno_file $pheno_file_name \\\n"+
"\t--pylmm $pylmmGWAS --outfile $gwas_outfile\n"


def parseargs():    # handle user arguments
	parser = argparse.ArgumentParser(description="Run a simple Pylmm GWAS" +
				" leave-one-chromosome-out analysis.")
	parser.add_argument('--clinical_traits', required = True,
		help = 'Clinitical trait tsv file.')
	parser.add_argument('--trait_name', required = True,
		help = 'EXACT name of the trait to study.')
	parser.add_argument('--output_dir', required = True,
		help = 'Which directory to output results to.')
	parser.add_argument('--include_sex_chromosomes', action='store_true',
		help = 'Whether to include sex chromosomes. Default: exclude.')
	args = parser.parse_args()
	return args


def write_qsub_script(args):
	"""
	Writes a temporary shell script to run the analysis as a job on hoffman2.
	After submitting the script, it will be deleted.
	Argument: args are the user arguments parsed from argparse.
	"""
	with(open('TEMP-QSUB-SCRIPT', 'w')) as outfile:
		outfile.write(SET_BLOCK_1)  # predifined parts of script
		outfile.write("clinical_traits='" + args.clinical_traits + "'\n")
		outfile.write("trait_name='" + args.trait_name + "'\n")
		outfile.write("output_directory='" + args.output_dir + "'\n")
		outfile.write("include_sex_chromosomes='" +
			args.include_sex_chromosomes + "'\n")
		outfile.write(SET_BLOCK_2)  # predifined parts of script


def main():
	args = parseargs()
	if not args.output_dir.endswith('/'):
		args.output_dir += '/'
	print('Launching a job on hoffman2 according to your specifications.')
	print('Your results will be stored in ' + args.output_dir)
	print('Plink results will be stored at ' + args.output_dir + 'plink12*')
	print('Kinship matrices will be stored at ' +args.output_dir + 'loco/*.kin')
	print('Pylmm GWAS results will be stored at ' + args.output_dir +
		'pylmm_gwas_results.txt')

	write_qsub_script(args)  # write temporary qsub script to run analysis
	# qsub the script and, after submitted, delete it
	subprocess.Popen(['qsub', 'TEMP-QSUB-SCRIPT']).wait()
	subprocess.Popen(['rm', 'TEMP-QSUB-SCRIPT']).wait()


if __name__ == '__main__':
	main()
