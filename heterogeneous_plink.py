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
  parser.add_argument('--genotypes', required = True,
    help = 'Folder containing downloaded heterogeneous stock mice genotypes.')
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


def parse_alleles(args):
  """
  The .alleles files in the genotypes folder contains the ordered SNP names
    for the ped files. We extract this information here.
  For more information on ped files, see plink documentation online.
  Arguments:
  -- args are the user arguments parsed from argparse.
  """
  marker_names = {}
  num_chroms = 19 + (args.include_sex_chromosomes==1)  # 20 chroms if X included
  for i in range(1, num_chroms+1):
    chrom_name = 'chr' + str(i)
    marker_names[chrom_name] = []
    with(open(args.genotypes + chrom_name + '.Build37.alleles', 'r')) as infile:
      for line in infile:
        # look for lines w/ marker names, parse them, add to list for this chrom
        if not line.startswith('marker '):
          continue
        marker_names[chrom_name].append(line.split(' ')[1])
  return marker_names


def reformat_map(args, marker_names, outpath):
  """
  Takes the mapfile.txt in the heterogeneous mice genotype folder, which is
    already nearly in plink map format, and puts it in the exact map format.
  For more information on map files, see plink documentation online.
  Arguments:
  -- args are the user arguments parsed from argparse.
  -- marker_names gives the ordered SNP names parsed from .alleles files
  -- outpath is the path for writing output.
  """
  num_chroms = 19 + (args.include_sex_chromosomes==1)  # 20 chroms if X included
  chrom_data = {}  # holds mapfile data, later written in order of marker_names
  for chrom_name in marker_names:
    chrom_data[chrom_name] = {}

  with(open(args.genotypes+'mapfile.txt', 'r')) as mapfile:
    with(open(outpath+'TEMP.map', 'w')) as outfile:
      mapfile.readline()  # skip header
      for line in mapfile:
        x_chromosome = []  # store x chromosome to write at end if desired
        splits = line.strip().split('\t')  # [SNP, chromosome, bp]
        if splits[1] == 'X':
          if not args.include_sex_chromosomes:
            continue
          splits[1] = '20'  # other files refer to this as chromosome 20
        # two layer dict storing line to print for a given snp in a given chrom.
        chrom_data['chr'+splits[1]][splits[0]] = '\t'.join(
            [splits[1], splits[0], '0', splits[2]]) + '\n'

      for i in range(1, num_chroms + 1):
        # we also need a mapfile for each chrom corresponding to chrom ped files
        with(open(outpath+'TEMP' + str(i) + '.map', 'w')) as map_chrom:
          chrom_name = 'chr' + str(i)
          for snp in marker_names[chrom_name]:
            if snp in chrom_data[chrom_name]:
              outfile.write(chrom_data[chrom_name][snp])
              map_chrom.write(chrom_data[chrom_name][snp])


def reformat_peds(args, outpath):
  """
  Removes mice from ped files that aren't in args.clinical file, and adds
    tab-delimiting between biallelic genotypes (keeps space within each pair).
  Also makes sure missing genotypes are formatted properly.
  For more information on ped files, see plink documentation online.
  Arguments:
  -- args are the user arguments parsed from argparse.
  -- outpath is the path for writing output.
  """
  # first parse clinical file to get the (ordered) list of mice
  mice_in_clinical = {}  # dict for instant name-based access
  with(open(args.clinical, 'r')) as clinical_file:
    header = clinical_file.readline().strip().split('\t')
    mousecol = header.index('mouse_number')
    for line in clinical_file:
      splits = line.strip().split('\t')
      mice_in_clinical[splits[mousecol]] = ''

  # now write the peds without the mice that are absent from the clinical file
  num_chroms = 19 + (args.include_sex_chromosomes==1)  # 20 chroms if X included
  for i in range(1, num_chroms + 1):
    chrom_name = 'chr' + str(i)
    with(open(args.genotypes + chrom_name + '.Build37.data', 'r')) as infile:
      with(open(outpath+'TEMP'+str(i)+'.ped', 'w')) as outfile:
        for line in infile:
          splits = line.strip().split(' ')
          if splits[1] not in mice_in_clinical:  # mouse not in args.clinical
            continue
          # convert missing genotype from 'NA' to plink default of '0'
          splits = splits[:6] + [i if i != 'NA' else '0' for i in splits[6:]]
          # fields 6 and onwards contain the alleles; we group them into
          #   biallelic genotypes then write tab-delim mouse fields & genotypes
          genos = [splits[i] + ' ' + splits[i+1]
                    for i in range(6, len(splits)) if i%2==0]
          for j in range(len(genos)):
            if '0' in genos[j]:
              genos[j] = '0 0'  # ensure no half-missing genotypes
          outfile.write('\t'.join(splits[:6]) + '\t' + '\t'.join(genos) + '\n')


def write_tfam(args, outpath):
  """
  Uses a ped file (which has already been subset and ordered according to
    args.clinical) to get the information needed to write the tfam file.
  For more information on tfam files, see plink documentation online.
  Arguments:
  -- args are the user arguments parsed from argparse.
  -- outpath is the path for writing output.
  """
  with(open(outpath + 'TEMP1.ped', 'r')) as infile:
    with(open(outpath+'TEMP.tfam', 'w')) as outfile:
      for line in infile:
        splits = line.split(' ')
        # tfam fields are just first 6 ped fields, so simply write them in order
        outfile.write(' '.join(splits[:6]) + '\n')


def transpose_and_combine(args, outpath):
  """
  Transpose the plink ped-formatted chromosome files into tped files, then
  	combine them into a single tped file.
  Then do the final recoding and maf/geno filtering.
  For more information on tped files, see plink documentation online.
  Arguments:
  -- args are the user arguments parsed from argparse.
  -- outpath is the path for writing output.
  """
  num_chroms = 19 + (args.include_sex_chromosomes==1)  # 20 chroms if X included
  individual_tpednames = []
  for i in range(1, num_chroms + 1):
    chrom_name = 'chr' + str(i)
    base_name = outpath + 'TEMP' + str(i)
    individual_tpednames.append(base_name + '.tped')
    subprocess.Popen(['plink', '--mouse', '--recode',
        '--transpose', '--file', base_name, '--out', base_name]).wait()

  # now that these are tpeds they can simply be catted together
  with(open(outpath + 'TEMP.tped', 'w')) as outfile:
    subprocess.Popen(['cat'] + individual_tpednames, stdout = outfile)

  # now recode (reformat) and maf/geno filter the tfam & tped files using plink
  subprocess.Popen(['plink', '--mouse', '--recode', '12', '--transpose',
            '--maf', args.maf, '--geno', args.geno, '--tfile', outpath+'TEMP',
            '--out', args.plink_basename]).wait()
  # clean temporary files
  tempfiles = glob.glob(outpath + 'TEMP*')
  subprocess.Popen(['rm'] + tempfiles).wait()


def main():
  args = parseargs()
  if not args.genotypes.endswith('/'):
    args.genotypes += '/'
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

  marker_names = parse_alleles(args)  # parse alleles files in genotypes folder
  reformat_map(args, marker_names, outpath)  # puts mapfile in plink .map format
  reformat_peds(args, outpath)  # remove entries not in clinical, tab-delimit
  write_tfam(args, outpath)  # write tfam file using ped file info
  transpose_and_combine(args, outpath)  # transpose & combine chrom. ped files


if __name__ == '__main__':
  main()
#
