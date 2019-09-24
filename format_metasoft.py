import argparse, glob, subprocess, sys


def parseargs():    # handle user arguments
  parser = argparse.ArgumentParser(description="Take multiple PyLMM results" +
        " and create a new file with those results in Metasoft format.")
  parser.add_argument('--pylmm', nargs = '+', required = True, help = 'PyLMM results files.')
  parser.add_argument('--output', required = True, help = 'Output file name.')
  parser.add_argument('--rsid2jax',
    default = '/u/home/n/nlapier2/mousedata/rsid2jax.txt',
    help= 'Location of rsid2jax.txt file. Default is a known location.')
  args = parser.parse_args()
  return args


# Read pylmm files, extracting beta and beta_sd, keep only SNPs present in all
def parse_pylmm(pylmm_files):
  # first read in the beta and beta_sd from each pylmm results file
  snp2info = {}
  for pylmm in pylmm_files:
    with(open(pylmm, 'r')) as infile:
      infile.readline()  # skip header
      for line in infile:
        splits = line.strip().split('\t')  # [snp_id, beta, beta_sd, f_stat, p_value]
        if splits[0] not in snp2info:
          snp2info[splits[0]] = [[splits[1], splits[2]]]
        else:
          snp2info[splits[0]].append([splits[1], splits[2]])

  # now intersect to make sure only SNPs in all results files are kept
  del_list = []
  for snp in snp2info:
    if len(snp2info[snp]) < len(pylmm_files):
      del_list.append(snp)
  for snp in del_list:
    del snp2info[snp]
  return snp2info


# Reads a file that translates between JAX SNP locations and rsids
def parse_rsid2jax(rsid2jax):
  jax2rsid = {}
  with(open(rsid2jax, 'r')) as infile:
    for line in infile:
      rsid, jax = line.strip().split()
      jax2rsid[jax] = rsid
  return jax2rsid


# write snp info in Metasoft input format
def write_metasoft(output, snp2info, jax2rsid):
  with(open(output, 'w')) as outfile:
    for snp in snp2info:
      if snp not in jax2rsid:
        continue
      rsid_snp_name = jax2rsid[snp]
      outfile.write(rsid_snp_name)
      for beta_pair in snp2info[snp]:
        outfile.write(' ' + str(beta_pair[0]) + ' ' + str(beta_pair[1]))
      outfile.write('\n')


def main():
  args = parseargs()
  snp2info = parse_pylmm(args.pylmm)  # map snp names to p-vals w/ pylmm results
  jax2rsid = parse_rsid2jax(args.rsid2jax)  # allows conversion from jax to rsid
  write_metasoft(args.output, snp2info, jax2rsid)  # write out in igv format


if __name__ == '__main__':
  main()
#
