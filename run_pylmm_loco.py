"""
This script runs pylmm LOCO analysis; that is, pylmmGWAS is run on the SNPs in
  in each chromosome using as input the kinship matrix with that chromosome
  left out. This helps avoid decreased power due to linkage disequilibrium
  effects. Results for each chromosome are then aggregated into a single file.
"""


import argparse, glob, os, subprocess, sys


def parseargs():    # handle user arguments
  parser = argparse.ArgumentParser(description="Given plink files," +
        " generate LOCO tpeds/maps and kinship matrices.")
  parser.add_argument('--loco_dir', required = True,
    help = 'Directory of LOCO files to run pylmmGWAS.py on.')
  parser.add_argument('--pheno_file', required = True,
    help = 'Phenotype file in pylmm format.')
  parser.add_argument('--outfile', default='AUTO',
    help = 'Name of final pylmm results file.')
  parser.add_argument('--pylmm',
    default = '/u/home/n/nlapier2/mousedata/pylmm/scripts/pylmmGWAS.py',
    help = 'Executable path to pylmmGWAS.py. Default is known location.') 
  args = parser.parse_args()
  return args


def run_pylmm_loco(args):
  """
  Runs pylmm on each of the one chromosome (one_chr) files, using the kinship
      matrix generated by leaving that one chromosome out (LOCO).
  Argument: args are the user arguments parsed with argparse.
  """
  num_chroms = 21  # max for mouse, we also allow fewer i.e. exluding X & Y
  for i in range(num_chroms):
    # the kinship with this chromsome (i) left out
    kin_name = args.loco_dir + 'loco_chr' + str(i+1) + '.kin'
    if not os.path.exists(kin_name):
      break  # subset of chromosomes available
    # the one_chr file for thid chromosome (i)
    one_chr_name = args.loco_dir + 'one_chr' + str(i+1)
    # execute pylmm and put results in same directory
    subprocess.Popen(['python', args.pylmm, '--tfile', one_chr_name,
              '--kfile', kin_name, '--phenofile', args.pheno_file,
              args.loco_dir+'gwas-chrom'+str(i+1)+'.txt']).wait()
    # clean one_chr files (no longer needed)
        # to_del = glob.glob(one_chr_name + '.*')
        # for fname in to_del:
            # subprocess.Popen(['rm', fname]).wait()


def aggregate_chrom_files(args):
  """
  Combine the pylmm results for each chromosome into one file.
  Argument: args are the user arguments parsed with argparse.
  """
  num_chroms = 21  # max for mouse, we also allow fewer i.e. exluding X & Y
  with(open(args.outfile, 'a')) as outfile:
    for i in range(num_chroms):
      gwas_outfile = args.loco_dir + 'gwas-chrom' + str(i+1) + '.txt'
      if not os.path.exists(gwas_outfile):
        break  # only subset of chromosomes available
      if i == 0:  # keep header from first file
        subprocess.Popen(['cat', gwas_outfile], stdout=outfile).wait()
      else:
        subprocess.Popen(['tail', '-n', '+2', gwas_outfile],
           stdout=outfile).wait()
      subprocess.Popen(['rm', gwas_outfile]).wait()


def main():
  args = parseargs()  # handle user arguments
  if not args.loco_dir.endswith('/'):
    args.loco_dir += '/'
  if args.outfile == 'AUTO':
    args.outfile = args.loco_dir + 'pylmm_gwas_results.txt'
  open(args.outfile, 'w').close()  # clear file and check if writeable

  run_pylmm_loco(args)  # run pylmm analysis on each chromosome
  aggregate_chrom_files(args)  # aggregate all resulting files into one file


if __name__ == '__main__':
  main()
#
