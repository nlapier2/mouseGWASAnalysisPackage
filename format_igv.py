import argparse, glob, subprocess, sys


def parseargs():    # handle user arguments
  parser = argparse.ArgumentParser(description="Take PyLMM results" +
        " and create a new file with those results in IGV format.")
  parser.add_argument('--pylmm', required = True, help = 'PyLMM results file.')
  parser.add_argument('--output', required = True, help = 'Output file name.')
  parser.add_argument('--all_strains',
    default = '/u/home/n/nlapier2/mousedata/all_strains.tped',
    help= 'Location of all_strains.tped file. Default is a known location.')
  args = parser.parse_args()
  return args


def parse_pylmm(pylmm):  # map snp names to p-vals w/ pylmm results
  snp2info = {}
  with(open(pylmm, 'r')) as infile:
    infile.readline()  # skip header
    for line in infile:
      splits = line.strip().split('\t')  # [snp_id, beta, beta_sd, f_stat, p_value]
      snp2info[splits[0]] = [splits[-1]]
  return snp2info


def parse_tped(all_strains, snp2info):  # add chr and bp to snp2info via tped
  with(open(all_strains, 'r')) as infile:
    infile.readline()  # skip header
    for line in infile:
      splits = line.strip().split('\t')  # [snp_chr, snp_id, centimorgans, snp_bp_mm10, ...]
      if splits[1] in snp2info:
        snp2info[splits[1]].extend([splits[0], splits[3]])
  return snp2info


def write_igv(output, snp2info):  # write out snp2info in igv format
  # note snp2info now contains snp_id: [p_value, snp_chr, snp_bp_mm10]
  with(open(output, 'w')) as outfile:
    outfile.write('CHR\tBP\tSNP\tP\n')  # header line
    for snp_id in snp2info:
      p_val, snp_chr, snp_bp_mm10 = snp2info[snp_id]
      outfile.write('\t'.join([snp_chr, snp_bp_mm10, snp_id, p_val]) + '\n')


def main():
  args = parseargs()
  snp2info = parse_pylmm(args.pylmm)  # map snp names to p-vals w/ pylmm results
  snp2info = parse_tped(args.all_strains, snp2info)  # add chr and bp from tped
  write_igv(args.output, snp2info)  # write out in igv format


if __name__ == '__main__':
  main()
#
