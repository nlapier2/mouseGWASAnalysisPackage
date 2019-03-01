"""
This script preprocesses a given clinical traits file, and should always be run
    upon downloading such a file from the SQL server. This file must be
  tab-delimited (see the wiki for more information).
The main idea is that phenotype and strain names will be standardized between
  all studies. This is done by reading in a manually curated mapping file
  (which we call pheno_map here) that maps names seen in the clinical trait
  files to the names we want to have.
Some other minor tweaks are made. The output is a new clinical_traits file with
    all these tweaks and substitutions.
This script is not considered "analysis" and is thus not included in the
    submit-pylmm-analysis.sh or simple-analysis.py wrapper scripts.
"""


import argparse, glob, subprocess, sys


def parseargs():    # handle user arguments
  parser = argparse.ArgumentParser(description="Prepare clincal trait" +
        " files for further processing and pylmm analysis.")
  parser.add_argument('--clinical', required = True,
    help = 'Clinitical trait tsv file.')
  parser.add_argument('--heterogeneous_genotypes', default = 'NONE/',
    help = 'genotypes folder for heterogeneous stock mice (if using them).')
  parser.add_argument('--outbred_dosages', default = 'NONE/',
    help = 'dosages folder for outbred mice (if using them).')
  parser.add_argument('--outname', default='AUTO', help = 'Output file name.')
  parser.add_argument('--pheno_map',
    default = '/u/home/n/nlapier2/mousedata/mouseGWASAnalysisPackage/pheno_map.txt',
    help = 'Specify phenotype map file to read. Default to known location.')
  args = parser.parse_args()
  return args


def read_pheno_map(mapfile):
  """
  Read in a file mapping original phenotype names to new names.
  This allows us to keep phenotype names consistent across traits.
  Argument: the args.pheno_map file specified by the user (or default "AUTO")
  Returns: dict mapping old phenotype names to new phenotype names
  """
  pheno_map = {}
  # Parse mapfile contents into a dict. The file contains
  #    tab separated original/new pairs, so splits[0] is the original and
  #    splits[1] is the new name. Dict thus maps original name to new name.
  with(open(mapfile, 'r')) as f:
    for line in f:
      if line.startswith('#'):
        continue
      splits = line.strip().split('\t')
      if len(splits) < 2:
        continue
      pheno_map[splits[0]] = splits[1]
  return pheno_map


def write_ordered_clinical(args, lines_to_write, straincol):
  """
  Write clinical file with mice in same order as ped files when using
    heterogeneous stock mice data, or nameList files when using outbred mice.
  Also makes heterogeneous FIDs consistent between ped files and pheno file.
  Arguments:
  -- args: the user arguments parsed by parseargs
  -- lines_to_write: lines from clinical file to write
  -- straincol: Column with strain (FID) field; make consistent
  """
  # extract ordered mouse names from either heterogeneous or outbred format
  ordered_mice = []  # ordered mice identified by iid (unique for this data)
  if args.heterogeneous_genotypes != 'NONE/':
    with(open(args.heterogeneous_genotypes + 'chr1.Build37.data', 'r')) as ped:
      for line in ped:
        splits = line.split(' ')
        ordered_mice.append([splits[0], splits[1]])
  elif args.outbred_dosages != 'NONE/':
    with(open(args.outbred_dosages + 'chr1_nameList.tsv', 'r')) as nameList:
      nameList.readline()  # first line is not needed
      all_iids = nameList.readline().strip().split('\t')  # contains all iids
      ordered_mice = ['Q_CFW-SW/' + iid.split('_recal')[0].split('_')[-1]
        for iid in all_iids]  # format iids in same way as clinical file
      ordered_mice = [[i, i] for i in ordered_mice]  # create fid/iid pairs

  # write clinical lines in order specified by ordered_mice
  with(open(args.outname, 'a')) as outfile:
    for fid, iid in ordered_mice:
      if iid in lines_to_write:
        lines_to_write[iid][straincol] = fid  # strain set to fid
        outfile.write('\t'.join(lines_to_write[iid]) + '\n')


def preprocess_traits(args, pheno_map):
  """
  Preprocess the clinical traits file for easier downstream usage.
  Mainly, replace the pheno_map strings and 'NULL' with 'NA' (for R).
  Arguments:
  -- args: the user arguments parsed by parseargs
  -- pheno_map: the phenotype map read in by the read_pheno_map method.
  """
  # this code block processes the header line with fid/iid and phenotype names
  lines_to_write = {}  # used to buffer lines for reordering if heterogeneous
  with(open(args.clinical, 'r')) as traitfile:
    header = traitfile.readline().strip().split('\t')
    # remove leading non-ascii characters, replace col names using pheno_map
    header[0] = ''.join([ch for ch in header[0] if ord(ch) < 128])
    if args.outbred_dosages != 'NONE/':
      header = ['Strain'] + header
    header = [i if i not in pheno_map else pheno_map[i] for i in header]
    # Determine columns that define mouse number and strain. If default
    #    phenotype map file is used, columns guaranteed to have these names.
    mousecol, straincol = header.index('mouse_number'), header.index('Strain')

    # this code block formats & replaces phenotype names according to pheno_map
    with(open(args.outname, 'w')) as outfile:
      outfile.write('\t'.join(header) + '\n')
      for line in traitfile:
        splits = line.strip().split('\t')
        if len(splits) < 2:  # trailing line at end of file
          break
        # no mouse_number in outbred; make same as strain name
        if args.outbred_dosages != 'NONE/':
          splits = [splits[mousecol-1]] + splits
        splits[straincol] = splits[straincol].replace(' ', '_')
        # replace strain name using pheno_map, if necessary
        if splits[straincol] in pheno_map:
          splits[straincol] = pheno_map[splits[straincol]]
        for i in range(len(splits)):  # replace NULL with NA, for R
          if splits[i] == 'NULL':
            splits[i] = 'NA'

        # this code block deals with writing the output
        # if using heterogeneous or outbred mice, store the lines to later write
		#    in the same order as they appear in the ped or nameList files
        if args.heterogeneous_genotypes != 'NONE/' or args.outbred_dosages != 'NONE/':
          lines_to_write[splits[mousecol]] = splits  # '\t'.join(splits) + '\n'
        else:  # otherwise just write the lines out
          outfile.write('\t'.join(splits) + '\n')
  if args.heterogeneous_genotypes != 'NONE/' or args.outbred_dosages != 'NONE/':
    write_ordered_clinical(args, lines_to_write, straincol)


def main():
  args = parseargs()
  if args.outname == 'AUTO':
    args.outname = 'preprocessed_clinical_traits.tsv'
  if not args.heterogeneous_genotypes.endswith('/'):
    args.heterogeneous_genotypes += '/'
  if not args.outbred_dosages.endswith('/'):
    args.outbred_dosages += '/'

  # Reads in a phenotype mapping file, uses this to replace phenotype names
  #    in the clincal traits file.
  pheno_map = read_pheno_map(args.pheno_map)
  preprocess_traits(args, pheno_map)


if __name__ == '__main__':
  main()
#
