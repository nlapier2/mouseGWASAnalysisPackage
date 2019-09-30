import argparse, glob, subprocess, sys


def parseargs():    # handle user arguments
  parser = argparse.ArgumentParser(description="Take PyLMM results" +
        " and create a new file with those results in IGV format.")
  parser.add_argument('--pylmm', required = True, help = 'PyLMM results file.')
  parser.add_argument('--output', required = True, help = 'Output file name.')
  parser.add_argument('--no_translate', action = 'store_true',
    help = 'Use if SNPs are already rsid.')
  parser.add_argument('--tped',
    default = '/u/home/n/nlapier2/mousedata/all_strains.tped',
    help= 'Location of tped file. Default is known all_strains.tped filepath.')
  parser.add_argument('--plink12_formatted', action = 'store_true',
    help = 'Use this flag if --tped is plink12 format.')
  parser.add_argument('--rsid2jax',
    default = '/u/home/n/nlapier2/mousedata/rsid2jax.txt',
    help= 'Location of rsid2jax.txt file. Default is a known location.')
  args = parser.parse_args()
  return args


def parse_pylmm(pylmm, jax2rsid, no_translate):  # map snp names to p-vals w/ pylmm results
  snp2info = {}
  with(open(pylmm, 'r')) as infile:
    infile.readline()  # skip header
    for line in infile:
      splits = line.strip().split('\t')  # [snp_id, beta, beta_sd, f_stat, p_value]
      if not no_translate:
        if splits[0] not in jax2rsid:
          continue
        snp2info[jax2rsid[splits[0]]] = [splits[-1]]
      else:
        snp2info[splits[0]] = [splits[-1]]
  return snp2info


# add chr and bp to snp2info via tped
def parse_tped(tped, jax2rsid, snp2info, plink12_formatted, no_translate):
  if plink12_formatted:
    delim = ' '
  else:
    delim = '\t'
  with(open(tped, 'r')) as infile:
    if not plink12_formatted:
      infile.readline()  # skip header
    for line in infile:
      splits = line.strip().split(delim)  # [snp_chr, snp_id, centimorgans, snp_bp_mm10, ...]
      if not no_translate:  # all_strains has to be converted to rsid
        if splits[1] not in jax2rsid:
          continue
        splits[1] = jax2rsid[splits[1]]
      if splits[1] in snp2info:
        snp2info[splits[1]].extend([splits[0], splits[3]])
  return snp2info


# Reads a file that translates between JAX SNP locations and rsids
def parse_rsid2jax(rsid2jax):
  jax2rsid = {}
  with(open(rsid2jax, 'r')) as infile:
    for line in infile:
      rsid, jax = line.strip().split()
      jax2rsid[jax] = rsid
  return jax2rsid


'''def write_igv(output, snp2info):  # write out snp2info in igv format
  # note snp2info now contains snp_id: [p_value, snp_chr, snp_bp_mm10]
  with(open(output, 'w')) as outfile:
    outfile.write('CHR\tBP\tSNP\tP\n')  # header line
    for snp_id in snp2info:
      p_val, snp_chr, snp_bp_mm10 = snp2info[snp_id]
      outfile.write('\t'.join([snp_chr, snp_bp_mm10, snp_id, p_val]) + '\n')'''


# Organizes SNPs by chromosome, then orders them by position in the chromosome
def get_ordered_snps(snp2info, jax2rsid):
  ordered_snps = {}
  for snp_id in snp2info:
    p_val, snp_chr, snp_bp_mm10 = snp2info[snp_id]
    rsid_snp_id = snp_id #jax2rsid[snp_id]
    if snp_chr not in ordered_snps:
      ordered_snps[snp_chr] = [[rsid_snp_id, snp_chr, int(snp_bp_mm10), p_val]]
    else:
      ordered_snps[snp_chr].append([rsid_snp_id, snp_chr, int(snp_bp_mm10), p_val])
  for chr in ordered_snps:  # order by location in snp (snp_bp_mm10)
    ordered_snps[chr].sort(key=lambda x: x[2])
  return ordered_snps


def write_igv(output, ordered_snps):  # write out ordered snp info in igv format
  with(open(output, 'w')) as outfile:
    outfile.write('SNP\tCHR\tBP\tP\n')  # header line
    for chrom in range(1, len(ordered_snps) + 1):
      if chrom == 20 and 'X' in ordered_snps:
        chrom = 'X'
      elif chrom == 21 and 'Y' in ordered_snps:
        chrom = 'Y'
      chrom_snps = ordered_snps[str(chrom)]
      for snp_info in chrom_snps:
        rsid_snp_id, snp_chr, snp_bp_mm10, p_val = snp_info
        if snp_chr == '20':
          snp_chr = 'X'  # required for IGV
        outfile.write('\t'.join([rsid_snp_id, snp_chr, str(snp_bp_mm10), p_val]) + '\n')


def main():
  args = parseargs()
  jax2rsid = parse_rsid2jax(args.rsid2jax)  # allows conversion from jax to rsid
  snp2info = parse_pylmm(args.pylmm, jax2rsid, args.no_translate)  # map snp names to p-vals
  snp2info = parse_tped(args.tped, jax2rsid, snp2info, args.plink12_formatted, args.no_translate)
  ordered_snps = get_ordered_snps(snp2info, jax2rsid)  # order by chrom then loc
  write_igv(args.output, ordered_snps)  # write out in igv format


if __name__ == '__main__':
  main()
#
