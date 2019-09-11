import argparse, glob, subprocess, sys


def parseargs():    # handle user arguments
  parser = argparse.ArgumentParser(description="Take PyLMM results" +
        " and create a new file with those results in IGV format.")
  parser.add_argument('--pylmm', required = True, help = 'PyLMM results file.')
  parser.add_argument('--output', required = True, help = 'Output file name.')
  parser.add_argument('--all_strains',
    default = '/u/home/n/nlapier2/mousedata/all_strains.tped',
    help= 'Location of all_strains.tped file. Default is a known location.')
  parser.add_argument('--rsid2jax',
    default = '/u/home/n/nlapier2/mousedata/rsid2jax.txt',
    help= 'Location of rsid2jax.txt file. Default is a known location.')
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
    rsid_snp_id = jax2rsid[snp_id]
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
      chrom_snps = ordered_snps[str(chrom)]
      for snp_info in chrom_snps:
        rsid_snp_id, snp_chr, snp_bp_mm10, p_val = snp_info
        if snp_chr == '20':
          snp_chr = 'X'  # required for IGV
        outfile.write('\t'.join([rsid_snp_id, snp_chr, str(snp_bp_mm10), p_val]) + '\n')


def main():
  args = parseargs()
  snp2info = parse_pylmm(args.pylmm)  # map snp names to p-vals w/ pylmm results
  snp2info = parse_tped(args.all_strains, snp2info)  # add chr and bp from tped
  jax2rsid = parse_rsid2jax(args.rsid2jax)  # allows conversion from jax to rsid
  ordered_snps = get_ordered_snps(snp2info, jax2rsid)  # order by chrom then loc
  write_igv(args.output, ordered_snps)  # write out in igv format


if __name__ == '__main__':
  main()
#
