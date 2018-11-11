# mouseGWASAnalysisPackage

Package for LOCO GWAS analysis of mouse clinical traits via pylmm

### Prerequisites

For the scripts in this package to work, you must have loaded the R and plink modules on hoffman2 and have an **executable** version of pylmm downloaded and set up (you can ensure this, for instance, with "chmod -R 755 pylmm/"). Also, make sure you unzip the all\_strains.tped.zip file.

### Usage

The scripts are all pretty simple to use. The python scripts are meant to be run by the user. Run "{scriptName}.py -h" to get a simple set of instructions.

### Notes on data analysis and scripts

The data is available at /u/home/n/nlapier2/scratch/mouse_data_lusis/studies

Some notes on the data analysis:
* Chow1 has NULL values for all of the traits we were going to study for all mice, so this study was excluded from further analysis
* For all mice in all studies, if one of the traits we were going to study is null, they all were, with one exception, a mouse in the hfhs_male study that had a NULL for GRAN % but non-null for the other traits. This mouse was excluded from analysis. For each study, mice with NULL values for the traits we're looking at were excluded.
* Each study (ath_female, ath_male, chow2, hfhs_female, hfhs_male, iron_female, iron_male) has its own directory under the studies directory, which contains the clinical traits file, clinical QTL file, intermediate plink analysis files, and loco_withsex and loco_nosex directories that contain all files for LOCO analysis, including LOCO tped/tfam/map files, one chromosome tped/tfam/map files, and the kinship matrices generated for each chromosome left out. The one chromosome files are used for the later GWAS analysis where one chromosome is tested against the kinship matrix generated with that chromosome held out.
* Kinships (and all other files) have been generated both with and without sex chromosomes being included.

The scripts are all pretty easy to use. I'd recommend running "python script.py -h" to get a quick overview of the options, but after seeing them it should be very clear what to do. Feel free to ask me if you have any questions though. I can add documentation if needed. Also, per the README, you'll need the R and plink modules loaded to run the trait setup script and an executable version of pylmm to run the kinship and gwas scripts.

Finished (I have already tested and run these, so you don't need to):
* trait_setup_and_plink.py: a script that takes a tab-separated clinical trait file from SQL and a specified trait, parses out the FIDs and IIDs and runs a slightly modified version of Calvin's get genotypes script on them, and then recodes this output to plink 12 format. Options to exclude sex chromosomes, name-change or delete intermediate files, specify location of the get genotypes script (assumes current directory). Get genotypes script must be in same folder as the unzipped all_genotypes.tped file.
* loco_kinship.py: Script that takes as input plink files, generates LOCO tped/tfam/map files and corresponding one chromosome files, then generates kinship matrices for each LOCO file via pylmm. Three ordered required arguments: plink files basename, output directory, path to the (executable) pylmmKinship.py script.

Possibly finished:
* run_pylmm_loco.py: a wrapper script for pylmmGWAS.py that takes as input your loco directory, a phenotype file, and path to the (executable) pylmmGWAS.py script, runs pylmm on each chromosome against the corresponding LOCO kinship matrix, and then aggregates the results to loco_dir/aggregated-gwas-results.txt. This script was cancelled by hoffman the first time I ran it in an interactive shell, but I think it was a time-out. When I ran it again last night it produced some unexpected stdout but the output seems to be as expected. Doing the LOCO analysis raised my previous correlations for chow1 fat mass form 84% to 91% and removed the consistent over-inflation. May be worth further testing this script to make sure it works. Could also be made a bit nicer.
* compare-pvals.py: Takes as input the rsid2jax.txt snp name mapping file included in the github, your pylmm pvalues file, and the clinical QTLs file from the SQL database, and outputs comparisons and charts. The script gives you the pearson and spearman correlations between the sets of pvalues, and also the number of uncorrected and naively bonferroni corrected significant hits for each. It also outputs two charts: one which plots the two pvalue sets in decreasing order of -log10(pval), which lets you check for general correlation and under/over inflation of one versus the other. The other is a zoomed in look at the top 1000 most significant hits for the SQL results, and the corresponding pylmm p-values. Obviously, since we don't have p-values for iron, this won't work for those. This script could be updated to make some other interesting graphs or do some more robust or interesting analysis, but otherwise it works.

Not finished:
* make_pheno_file.py: This is the script that generates the phenotypes value in pylmm format, so that pylmmGWAS.py can be run. This script already handles basic non-regressed phenotypes well: it will take in your clinical traits file, read in and normalize the phenotype values, and exclude NULL mice or those not in the tfam file. In progress was a method that cleans the other traits to regress on: removing those that have too many NULLs and inserting the mean for limited NULLs in other traits. By default, it regresses against everything else available; the option to choose needs to be added. It also does not perform the actual regression yet; this needs to be done.
.
