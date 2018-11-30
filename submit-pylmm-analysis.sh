#!/bin/sh
#$ -S /bin/bash
#$ -N pylmm_analysis_qsub_job
#$ -cwd
#$ -o stdout-pylmm-analysis.out
#$ -l h_data=10G,h_rt=24:00:00


# This script is a wrapper for all pylmm analysis scripts, i.e.
#  	get_genotypes_and_plink.py, make_pheno_file.py, loco_kinship.py, and
#  	run_pylmm_loco.py.
# There are a few options you should set for every new
#  	analysis you perform, and some that you can usually leave the same.
# For more information on these options, see the wiki.


. /u/local/Modules/default/init/modules.sh
module load R
module load plink


# Important options: specify location of clinical traits file, EXACT trait name,
#  	output directory, and whether to include sex chromosomes in the analysis.
# You will probably want to check/modify these for every analysis you run.

clinical_traits='/full/path/to/my_clinical_traits_file.tsv'
trait_name='my_trait_name'  # i.e. Fat_mass
output_directory='/full/path/to/my/desired/output/directory/'
include_sex_chromosomes=1  # 1 to include sex chromosomes, 0 to exclude


# Other options that you may want to set, but often the default is fine.

all_strains='/u/home/n/nlapier2/mousedata/all_strains.tped'
plink_basename=$output_directory+'plink12'
pheno_file_name=$output_directory+'pheno_file_pylmm.txt'
no_normalization=1  # 1 to normalize phenotypes, 0 for no normalization
loco_outdir=$output_directory+'loco/'
pylmm_loc='/u/home/n/nlapier2/mousedata/pylmm/'
gwas_outfile=$loco_outdir+'pylmm_gwas_results.txt'


# Generally do not change the part below -- it simply runs the analysis scripts
#  	based on the arguments you specified above. Only change if you know what
#  	you are doing. One reason would be not wanting to run all analysis steps.

if [ $include_sex_chromosomes = 1 ]
then
	sex_chrom_opt='--no_sex_chromosomes'
else
	sex_chrom_opt=''
fi
if [ $no_normalization = 1 ]
then
	normalize_opt='--no_normalization'
else
	normalize_opt=''
fi
tfam=$plink_basename+'.tfam'
pylmmKinship=$pylmm_loc+'pylmmKinship.py'
pylmmGWAS=$pylmm_loc+'pylmmGWAS.py'

python get_genotypes_and_plink.py $clinical_traits --all_strains $all_strains \
	$sex_chrom_opt --plink_basename $plink_basename
python make_pheno_file.py --clinical $clinical_traits --target $trait_name \
	--tfam $tfam --output $pheno_file_name $normalize_opt
python loco_kinship.py --plink $plink_basename --outdir $loco_outdir \
	--pylmm $pylmmKinship
python run_pylmm_loco.py --loco_dir $loco_outdir --pheno_file $pheno_file_name \
	--pylmm $pylmmGWAS --outfile $gwas_outfile
#
