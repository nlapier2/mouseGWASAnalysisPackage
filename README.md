# mouseGWASAnalysisPackage

This is a package for leave-one-chromosome-out (LOCO) GWAS analysis of mouse clinical traits via pylmm on hoffman2. The only dependencies are numpy and scipy, which are pylmm dependencies. These are probably already installed by default, but if you have problems you may want to run "pip install --user numpy scipy".

### Quick Usage

```
git clone https://github.com/nlapier2/mouseGWASAnalysisPackage.git
cd mouseGWASAnalysisPackage
python simple-analysis.py --clinical /path/to/my/clinical_traits_file.tsv --trait_name Fat_mass --output_dir /path/to/directory/to/output/results/to/
```

This is the simplest way to run things, but it is recommended to run submit-pylmm-analysis.sh instead of simple-analysis.py if you are comfortable with editing files via command line (more details below).


### Pipeline Overview and Examples

In general, the pipeline is to download a tab-delimited clinical traits file from the HMDP database, run preprocess.py on it, then edit and qsub submit-pylmm-analysis.sh. The latter script is really a wrapper that runs four analysis scripts: get_genotypes_and_plink.py, make_pheno_file.py, loco_kinship.py, and run_pylmm_loco.py. The user should modify the wrapper script with the specific variables they want (their preprocessed clinical traits file, trait name, and desired output directory, then qsub it. 

A concrete example is shown below, analyzing the Fat_mass trait in the chow2 study with sex chromosomes included, with results sent to /my/output/directory/. When editing the submit-pylmm-analysis.sh script with vim below, assume we change clinical_traits to /my/out/dir/preprocessed_clinical_traits_chow2.tsv, trait_name to Fat_mass, output_directory to /my/out/dir/, and include_sex_chromosomes to 1.

```
python preprocessing.py --clinical /path/to/my/clinical_traits_chow2.tsv --outname /my/out/dir/preprocessed_clinical_traits_chow2.tsv
vim submit-pylmm-analysis.sh  # change variables as described above
qsub submit-pylmm-analysis.sh
```

The simple-analysis.py script is an even simpler wrapper that allows the user to define a few options via the command line; the script then writes and submits the bash script for the user and tells them where the outputs will be. Below is the previous example except using the simple-analysis.py wrapper.

```
python preprocessing.py --clinical /path/to/my/clinical_traits_chow2.tsv --outname /my/out/dir/preprocessed_clinical_traits_chow2.tsv
python simple-analysis.py --clinical /my/out/dir/preprocessed_clinical_traits_chow2.tsv --trait_name Fat_mass --output_dir /my/out/dir/
# Check back on /my/out/dir/ in a few hours; if you see /my/out/dir/pylmm_gwas_results.txt, the job has finished running.
``` 


Finally, users can choose to run each of the four analysis scripts individually. Here's that same example again:

```
python preprocessing.py --clinical /path/to/my/clinical_traits_chow2.tsv --outname /my/out/dir/preprocessed_clinical_traits_chow2.tsv
python get_genotypes_and_plink.py --clinical /path/to/my/clinical_traits_chow2.tsv --plink_basename /my/out/dir/plink12
python make_pheno_file.py --clinical /path/to/my/clinical_traits_chow2.tsv --target Fat_mass --tfam /my/out/dir/plink12.tfam --output /my/out/dir/pheno_file_pylmm.txt
python loco_kinship.py --plink /my/out/dir/plink12 --outdir /my/out/dir/loco/
python run_pylmm_loco.py --loco_dir /my/out/dir/loco/ --pheno_file /my/out/dir/pheno_file_pylmm.txt --outfile /my/out/dir/pylmm_gwas_results.txt
```

The following section describes each script in more detail. For documentation of the actual code, please see the source files. Note that the arguments for each script can be seen on the command line by running "python some_script.py -h".


### Script Details

#### preprocessing.py

Overview: This script preprocesses a given clinical traits file, and should always be run upon downloading such a file from the SQL server. This file must be tab-delimited (see the wiki for more information). The main idea is that phenotype and strain names will be standardized between all studies. This is done by reading in a manually curated mapping file (which we call pheno_map here) that maps names seen in the clinical trait files to the names we want to have. Some other minor tweaks are made. The output is a new clinical_traits file with all these tweaks and substitutions. This script is not considered "analysis" and is thus not included in the submit-pylmm-analysis.sh or simple-analysis.py wrapper scripts.

Example:
```
python preprocessing.py --clinical /path/to/my/clinical_traits_chow2.tsv --outname /my/out/dir/preprocessed_clinical_traits_chow2.tsv
```

Arguments:
* (REQUIRED) --clinical: Clinitical trait tsv file.
* (DEPRECATED) --qtls: Clinitical QTLs tsv file. Do not use unless you know what you are doing.
* --outname: Output file name.
* --pheno_map: Specify phenotype map file to read. Default to known location.

#### submit-pylmm-analysis.sh

Overview: This script is a wrapper for all pylmm analysis scripts, i.e. get_genotypes_and_plink.py, make_pheno_file.py, loco_kinship.py, and run_pylmm_loco.py. There are a few options you should set for every new analysis you perform, namely clinical_traits, trait_name, output_directory, and include_sex_chromosomes. There are many other options that you can usually leave to the default, but you may occasionally want to change.

Example:
```
vim submit-pylmm-analysis.sh  # change variables as appropriate (see example in pipeline overview section)
qsub submit-pylmm-analysis.sh
```

Arguments (in this case, not command line arguments, but variables to modify in the script):
* (REQUIRED) clinical_traits: /full/path/to/my_clinical_traits_file.tsv
* (REQUIRED) trait_name: EXACT name of the trait to study as it appears in the (preprocessed) clinical traits file. Example: "Fat_mass".
* (REQUIRED) output_directory: /full/path/to/my/desired/output/directory/
* (REQUIRED) include_sex_chromosomes: 1 to include sex chromosomes, 0 to exclude
* all_strains: Location of all_strains.tped file, which matches mouse strains to their genotypes. Unless you know what you are doing, leave this to the default.
* plink_basename: Base name (i.e. do not specify file extension) of the output plink files from get_genotypes_and_plink.py
* pheno_file_name: Where to write the pylmm-formatted phenotype file.
* no_normalization: to normalize phenotypes, 0 for no normalization
* loco_outdir: Where to write the LOCO kinship matrices and other LOCO files.
* pylmm_loc: Location of the pylmm top level directory. Unless you know what you are doing, leave this to the default.
* gwas_outfile: Where to write the final pylmm LOCO GWAS results.

#### simple-analysis.py

Overview: This script is made to be a very simple wrapper to run the GWAS analysis steps in get_genotypes_and_plink.py, make_pheno_file.py, loco_kinship.py, and run_pylmm_loco.py. This script may be useful for users who are uncomfortable with editing files from the command line, which is necessary for running the submit-pylmm-analysis.sh script. However, users who are comfortable editing files should use the shell script for more flexibility.

Example:
```
python simple-analysis.py --clinical /my/out/dir/preprocessed_clinical_traits_chow2.tsv --trait_name Fat_mass --output_dir /my/out/dir/
```

Arguments:
* (REQUIRED) --clinical: Clinitical trait tsv file.
* (REQUIRED) --trait_name: EXACT name of the trait to study.
* (REQUIRED) --output_dir: Which directory to output results to.
* --include_sex_chromosomes: Whether to include sex chromosomes. Default: exclude.

#### get_genotypes_and_plink.py

Overview: This script takes the mice from a specified clinical traits file and matches them to their genotype information using a file called all_strains.tped. This information is then output in plink format, and then plink is used to recode alleles into the numbers 1 and 2, for GWAS analysis purposes. This is the first step of analysis, and importantly involves the decision of whether to include sex chromosomes or not.

Example: 
```
python get_genotypes_and_plink.py --clinical /path/to/my/clinical_traits_chow2.tsv --plink_basename /my/out/dir/plink12
```

Arguments:
* (REQUIRED) --clinical: Clinical trait tsv file.
* --all_strains: Location of all_strains.tped file. Default is a known location.
* --no_sex_chromosomes: Ignore sex chromosomes.
* --plink_basename: Base name of plink output.


#### make_pheno_file.py

Overview: This script makes a phenotype file in pylmm format, which pylmm then uses to read in the phenotype values for each mouse for GWAS analysis. Essentially, we are simply extracting relevant information from the clinical traits file. One of the plink files called a tfam is required input, as it specifies which mice are present and genotyped for this study. This script also normalizes phenotypes by subtracting their mean and dividing by their standard deviation. This can be turned off, but generally shouldn't be. Eventually this script will allow the user to specify covariates to regress out, but we currently do not have that option implemented.

Example:
```
python make_pheno_file.py --clinical /path/to/my/clinical_traits_chow2.tsv --target Fat_mass --tfam /my/out/dir/plink12.tfam --output /my/out/dir/pheno_file_pylmm.txt
```

Arguments:
* (REQUIRED) --clinical: Name of existing pheno file.
* (REQUIRED) --target: EXACT name of target trait to study.
* (REQUIRED) --tfam: Plink tfam file; exclude mice not in it.
* output: Name of output file.
* no_normalization: Do not normalize target phenotype (flag only; no text).

#### loco_kinship.py

Overview: This script generates leave one chromosome out (LOCO) kinship matrices via pylmm; pylmm by itself doesn't have the LOCO option. The kinship matrix measures the genetic relatedness between mice in the study and is a critical component of GWAS analysis. In general, the kinship matrix is computed by comparing all SNPs between all mice. In LOCO analysis, a kinship matrix is generated for each chromosome, where that kinship matrix only looks at the SNPs NOT in that chromosome. This helps avoid decreased power due to linkage disequilibrium effects. We also generate files containing the chromosome excluded from each kinship matrix; these are used for the pylmmGWAS analysis.

Example:
```
python loco_kinship.py --plink /my/out/dir/plink12 --outdir /my/out/dir/loco/
```

Arguments:
* (REQUIRED) --plink: Base name of plink files.
* --outdir: Name of output directory for LOCO and kinship files.
* --pylmm: Executable path to pylmmKinship.py. Default is known location.

#### run_pylmm_loco.py

Overview: This script runs pylmm LOCO analysis; that is, pylmmGWAS is run on the SNPs in in each chromosome using as input the kinship matrix with that chromosome left out. This helps avoid decreased power due to linkage disequilibrium effects. Results for each chromosome are then aggregated into a single file.

Example:
```
python run_pylmm_loco.py --loco_dir /my/out/dir/loco/ --pheno_file /my/out/dir/pheno_file_pylmm.txt --outfile /my/out/dir/pylmm_gwas_results.txt
```

Arguments:
* (REQUIRED) --loco_dir: Directory of LOCO files to run pylmmGWAS.py on.
* (REQUIRED) --pheno_file: Phenotype file in pylmm format.
* --outfile: Name of final pylmm results file.
* --pylmm: Executable path to pylmmGWAS.py. Default is known location.

#### compare-pvals.py

Not part of the analysis pipeline. Do not run this unless you know what you are doing. Also note that this script also requires matplotlib.

Overview: Compares p-values between SQL QTL file and pylmm results. Shows corerlations and significant hits, and makes some plots.

Example: 
```
python compare-pvals.py --pylmm /my/out/dir/pylmm_gwas_results.txt --qtls /path/to/my/clinical_QTL_chow2.tsv --trait_name Fat_mass
```

Arguments:
* (REQUIRED) --pylmm: pylmm pheno file.
* (REQUIRED) --qtls: Clinical QTLs from SQL DB.
* (REQUIRED) --trait_name: Name of trait being studied, needed for SQL QTL file.
* --fastlmm: use this boolean flag if the qtl file is actually raw fastlmm output.
* --out_basename: Plot output base names.

.
