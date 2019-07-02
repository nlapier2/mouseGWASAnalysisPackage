"""
The equivalent of get_genotypes_and_plink.py for Jonathan Flint's heterogeneous
	stock mice data.
"""


import argparse, glob, subprocess, sys


def parseargs():    # handle user arguments
  parser = argparse.ArgumentParser(description="Given phenotype file," +
        " get into proper format for pylmm.")
  parser.add_argument('--clinical', required = True,
    help = 'Clinical trait tsv file.')
  parser.add_argument('--geno_file', required = True,
    help = 'Original genotypes file.')
  parser.add_argument('--map_file', required = True,
    help = 'Original SNP map file.')
  parser.add_argument('--geno', default='0.1',
    help = 'Filter variants missing above this rate. Generally, do not modify.')
  parser.add_argument('--maf', default='0.05',
    help = 'Minor allele freq. filter, default 0.05. Generally, do not modify.')
  parser.add_argument('--include_sex_chromosomes', action='store_true',
    help = 'Include sex chromosomes.')
  parser.add_argument('--plink_basename', default='AUTO',
    help = 'Base name of plink output.')
  args = parser.parse_args()
  return args


def get_fid_iid(args):
  """
  Extract fid (family ID == mouse Strain) & iid (individual ID = mouse_number)
      from clinical traits file.
  Argument: args are the user arguments parsed from argparse.
  Returns: fid_iid, the strain / mouse number pairs from clinical traits file.
  """
  fid_iid = []  # fid == mouse strain, iid == individual mouse number
  with(open(args.clinical, 'r')) as infile:
    header = infile.readline().strip().split('\t')
    # determine columns with fid (strain) and iid (mouse number) info
    straincol, numcol = header.index('Strain'), header.index('mouse_number')
    for line in infile:
      splits = line.strip().split('\t')
      fid, iid = splits[straincol].replace(' ', '_'), splits[numcol]
      fid_iid.append([fid, iid])
  return fid_iid


def read_map(mapfile):
  """
  Read args.map_file to get the mapping from SNP name to chromosome & location.
  Argument: map.txt file from the downloaded Palmer data.
  Returns: dict mapping SNP name to chromosome & location, list of SNPs in order
  """
  snp_to_info, ordered_snps = {}, []
  with(open(mapfile, 'r')) as infile:
    infile.readline()  # skip header
    for line in infile:
      snp_name, chrom, pos = line.split()[:3]
      snp_to_info[snp_name] = [chrom, pos]
      ordered_snps.append(snp_name)
  return snp_to_info, ordered_snps


def make_tped(args, fid_iid, snp_to_info, ordered_snps, outpath):
  """
  Go through pruned_dosages (genotype dosage) file for each chromosome, removing
    columns for mice not in args.clinical, converting dosages to genotypes based
    on a heuristic, and formatting in plink 12 tped format.
  For more information on tped files, see plink documentation online.
  Arguments:
  -- args are the user arguments parsed from argparse.
  -- fid_iid is the list of FIDs and IIDs read from the args.clinical file
  -- snp_to_info maps SNP name to chromosome and position
  -- outpath is the path for writing output.
  Returns: the FIDs and IIDs for mice in both geno and pheno file originally
  """
  final_fid_iid = []  # includes mice in both geno and pheno file originally
  snp_to_dosages = {}  # maps SNP to the dosage values per mouse
  with(open(args.geno_file, 'r')) as infile:
    with(open(outpath + 'TEMP.tped', 'w')) as outfile:
      snp_list = infile.readline().strip().split()[2:]
      for snp in snp_list:
        snp_to_dosages[snp] = {}

      for line in infile:
        splits = line.strip().split()
        fid, discard = splits[0], splits[1]  # fid == iid here
        if discard == 'yes':
          continue
        if [fid, fid] not in fid_iid:
          continue
        final_fid_iid.append([fid, fid])  # now we know it's in both files
        for i in range(2, len(splits)):
          dosage, snp = float(splits[i]), snp_list[i - 2]
          # round obvious dosages, coding ambiguous ones as unknown
          if dosage > 1.9:
            snp_to_dosages[snp][fid] = '2 2'
          elif dosage < 0.1:
            snp_to_dosages[snp][fid] = '1 1'
          elif dosage > 0.9 and dosage < 1.1:
            snp_to_dosages[snp][fid] = '2 1'
          else:  # too ambiguous -- give '0 0' code which is "unknown" in plink
            snp_to_dosages[snp][fid] = '0 0'

      # now write the SNPs out using snp_to_info & ordered mice in final_fid_iid
      for snp in ordered_snps:
        if snp not in snp_to_dosages:
          continue
        chrom, pos = snp_to_info[snp]
        ordered_dosages = [snp_to_dosages[snp][fid[0]] for fid in final_fid_iid]
        outfile.write('\t'.join([chrom, snp, '0', pos] + ordered_dosages) + '\n')
  return final_fid_iid


def write_tfam(fid_iid, outpath):
  """
  Given the mouse strains (fid) and individual numbers (iid), write out a
      plink-format tfam file. In this case we don't specify the other tfam
    fields, so this is mostly just writing out the mice in the study.
  For more information on tfam files, see plink documentation online.
  Argument: fid_iid specifies mouse strains and numbers in this study.
  """
  other_cols = '\t0\t0\t0\t-9\n'  # father/mother/sex/pheno unspecified
  with(open(outpath+'TEMP.tfam', 'w')) as outfile:
    for pair in fid_iid:  # pair is [fid, iid] of a mouse
      outfile.write(pair[0] + '\t' + pair[1] + other_cols)


def run_plink_recode(args, outpath):
  """
  Run plink recode for final formatting and geno/maf filtering.
  Arguments:
  -- args are the user arguments parsed from argparse.
  -- outpath is the path for writing output.
    """
  num_chroms = 19 + (args.include_sex_chromosomes==1)  # 20 chroms if X included
  # now recode (reformat) and maf/geno filter the tfam & tped files using plink
  subprocess.Popen(['plink', '--mouse', '--recode', '12', '--transpose',
            '--maf', args.maf, '--geno', args.geno, '--tfile', outpath+'TEMP',
            '--out', args.plink_basename]).wait()
  # clean temporary files
  tempfiles = glob.glob(outpath + 'TEMP*')
  subprocess.Popen(['rm'] + tempfiles).wait()


def main():
  args = parseargs()
  if args.plink_basename == 'AUTO':
    # default basename is 'plink12' in same directory as args.clinical
    clinical_rpt_dir = ''
    splits = args.clinical.split('/')
    if len(splits) > 1:  # clinical report not in current dir
      clinical_dir = '/'.join(splits[:-1]) + '/'
    args.plink_basename = clinical_dir + 'plink12'
  if '/' not in args.plink_basename:
    outpath = './'
  else:
    outpath = '/'.join(args.plink_basename.split('/')[:-1]) + '/'

  fid_iid = get_fid_iid(args)  # get fid & iid values for mice in clinical file
  snp_to_info, ordered_snps = read_map(args.map_file)  # read snp info from map
  # here we make the tped file and remove mice from fid_iid not in that file
  fid_iid = make_tped(args, fid_iid, snp_to_info, ordered_snps, outpath)
  write_tfam(fid_iid, outpath)  # write temporary tfam file using fid_iid
  run_plink_recode(args, outpath)  # final plink formatting and filtering


if __name__ == '__main__':
  main()
#
