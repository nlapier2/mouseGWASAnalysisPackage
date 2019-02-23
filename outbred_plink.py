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
  parser.add_argument('--dosages', required = True,
    help = 'Folder containing extracted tsv-format outbred mice dosages.')
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


def find_cols_to_keep(args):
  """
  We only want to keep columns (mice) in the pruned_dosages files that are in
    args.clinical. This method finds those columns.
  Arguments:
  -- args are the user arguments parsed from argparse.
  """
  # in clincal file, extract mouse names from mouse_number column
  clinical_mice = []
  with(open(args.clinical, 'r')) as clinical_file:
    mousecol= clinical_file.readline().strip().split('\t').index('mouse_number')
    for line in clinical_file:
      clinical_mice.append(line.strip().split('\t')[mousecol])

  # extract mouse names from nameList and find which ones are in clinical_mice
  cols_to_keep = []
  with(open(args.dosages + 'chr1_nameList.tsv', 'r')) as nameList:
    nameList.readline()  # first line not needed
    all_iids = nameList.readline().strip().split('\t')  # contains all iids
    all_iids = ['Q_CFW-SW/' + iid.split('_recal')[0].split('_')[-1]
      for iid in all_iids]  # format iids in same way as clinical file
    cols_to_keep = {i: 0 for i in range(len(all_iids))
      if all_iids[i] in clinical_mice}
  return cols_to_keep


def write_tfam(args, outpath):
  """
  Write tfam file using the mouse names from args.clinical.
  For more information on tfam files, see plink documentation online.
  Arguments:
  -- args are the user arguments parsed from argparse.
  -- outpath is the path for writing output.
  """
  other_cols = '\t0\t0\t0\t-9\n'  # father/mother/sex/pheno unspecified
  with(open(args.clinical, 'r')) as infile:
    with(open(outpath+'TEMP.tfam', 'w')) as outfile:
      mousecol = infile.readline().strip().split('\t').index('mouse_number')
      for line in infile:
        mouse_number = line.strip().split('\t')[mousecol]
        # no family id, set to same as mouse_number
        outfile.write(mouse_number + '\t' + mouse_number + other_cols)


def pruned_dosages_to_tped(args, outpath, cols_to_keep):
  """
  Go through pruned_dosages (genotype dosage) file for each chromosome, removing
    columns for mice not in args.clinical, converting dosages to genotypes based
    on a heuristic, and formatting in plink 12 tped format.
  For more information on tped files, see plink documentation online.
  Arguments:
  -- args are the user arguments parsed from argparse.
  -- outpath is the path for writing output.
  -- cols_to_keep is the list of columns to keep from the pruned_dosages file
  """
  num_chroms = 19 + (args.include_sex_chromosomes==1)  # 20 chroms if X included
  for i in range(1, num_chroms + 1):
    chrom_name = 'chr' + str(i)
    if i == 20:
        chrom_name = 'chrX'
    with(open(args.dosages + chrom_name + '_pruned_dosages.tsv','r')) as dfile:
      with(open(args.dosages + chrom_name + '_pruned_pos.tsv','r')) as posfile:
        dfile.readline() ; posfile.readline()  # skip headers
        with(open(outpath+'TEMP'+str(i)+'.tped', 'w')) as outfile:
          for dose_line in dfile:
            # extract corresponding pruned_pos line and format for tped
            pos_line = posfile.readline().strip().split('\t')
            if pos_line[0] == '"chrX"':
              pos_line[0] = '"chr20"'
            pos_line[0] = pos_line[0].strip('"')
            tped_fields = [pos_line[0][3:], pos_line[1] + pos_line[0], '0',
              pos_line[1]]  # chr1 6350458 "T" "A" --> 1 6350458chr1 0 6350458

            # reformat allelic dosages to plink 12 format and write output
            dose_line = dose_line.strip().split('\t')
            dose_line = [dose_line[i] for i in range(len(dose_line))
              if i in cols_to_keep]  # remove cols for mice not in clinical file
            for i in range(len(dose_line)):
              dose = float(dose_line[i])
              if dose > 0.95:  # likely both minor alleles (plink 12 code is 2 2)
                dose_line[i] = '2 2'
              elif dose < 0.05:  # likely both major alleles (1 1 in plink 12)
                dose_line[i] = '1 1'
              elif dose > 0.45 and dose < 0.55:  # heterozygous --> 2 1
                dose_line[i] = '2 1'
              else:  # ambiguous dosage, we just code as missing (0 0 in plink 12)
                dose_line[i] = '0 0'
            all_fields = tped_fields + dose_line
            outfile.write('\t'.join(all_fields) + '\n')


def combine_and_filter(args, outpath):
  """
  Combine individual chromosome tped files, then do maf/geno filtering.
  For more information on tped files, see plink documentation online.
  Arguments:
  -- args are the user arguments parsed from argparse.
  -- outpath is the path for writing output.
    """
  num_chroms = 19 + (args.include_sex_chromosomes==1)  # 20 chroms if X included
  individual_tpednames = [outpath + 'TEMP' + str(i) + '.tped'
    for i in range(1, num_chroms + 1)]  # all tped file names
  with(open(outpath + 'TEMP.tped', 'w')) as outfile:  # cat together
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
  if not args.dosages.endswith('/'):
    args.dosages += '/'
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

  # keep only the columns (mice) in pruned_dosages that are in args.clinical
  cols_to_keep = find_cols_to_keep(args)
  write_tfam(args, outpath)  # write plink tfam file based on clinical file
  pruned_dosages_to_tped(args, outpath, cols_to_keep)  # reformat to tped files
  combine_and_filter(args, outpath)  # combine chrom tpeds and maf+geno filter


if __name__ == '__main__':
  main()
#
