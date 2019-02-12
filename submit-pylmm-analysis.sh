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

package_dir='/location/of/mouseGWASAnalysisPackage/'
clinical_traits='/full/path/to/my_clinical_traits_file.tsv'
trait_name='my_trait_name'  # i.e. Fat_mass
output_directory='/full/path/to/my/desired/output/directory/'
include_sex_chromosomes=0  # 1 to include sex chromosomes, 0 to exclude


# Other options that you may want to set, but often the default is fine.
all_strains='/u/home/n/nlapier2/mousedata/all_strains.tped'
plink_basename=$output_directory'plink12'
pheno_file_name=$output_directory'pheno_file_pylmm.txt'  
#quantile is quantile transformation, standardize is standardization, boxcox is box-cox transformation and none is none
transform='standardize'
plot=0 # 1 to plot transformations, 0 to forgo; if you plot, you cannot run the full pipeline after make_pheno_file.py
loco_outdir=$output_directory'loco/'
pylmm_loc='/u/home/n/nlapier2/mousedata/pylmm/'
gwas_outfile=$loco_outdir'pylmm_gwas_results.txt'

# Generally do not change the part below -- it simply runs the analysis scripts
#  	based on the arguments you specified above. Only change if you know what
#  	you are doing. One reason would be not wanting to run all analysis steps.

#if $plot = 1, then transformations are ignored
#if $plot = 1, the figure will be written to ${pheno_file_name_prefix}.png
if [ $plot = 1 ]
then
	plot_opt='--plot'
else
	plot_opt=''
fi

if [ $include_sex_chromosomes = 0 ]
then
	sex_chrom_opt='--no_sex_chromosomes'
else
	sex_chrom_opt=''
fi

tfam=$plink_basename'.tfam'
pylmmKinship=$pylmm_loc'scripts/pylmmKinship.py'
pylmmGWAS=$pylmm_loc'scripts/pylmmGWAS.py'
get_geno_loc=$package_dir'get_genotypes_and_plink.py'
make_pheno_loc=$package_dir'make_pheno_file.py'
loco_kinship_loc=$package_dir'loco_kinship.py'
run_pylmm_loc=$package_dir'run_pylmm_loco.py'

python $get_geno_loc --clinical $clinical_traits --all_strains $all_strains \
	$sex_chrom_opt --plink_basename $plink_basename
python $make_pheno_loc --clinical $clinical_traits \
	--target $trait_name --tfam $tfam --output $pheno_file_name --transform $transform $plot_opt
#these scripts will not run if you plot in the figure before
python $loco_kinship_loc --plink $plink_basename \
	--outdir $loco_outdir --pylmm $pylmmKinship
python $run_pylmm_loc --loco_dir $loco_outdir \
	--pheno_file $pheno_file_name --pylmm $pylmmGWAS --outfile $gwas_outfile
#
