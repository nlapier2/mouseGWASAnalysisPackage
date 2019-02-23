# Script for extracting genotypes from Rdata files for the outbred mice

# parse argument; should just be one argument to the outbred dosages directory
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop(paste("Requires a user argument to the outbred mice dosages directory",
               "i.e. extract_outbred_rdata.R dosages/"))
}
geno_dir <- args[1]
if (!(endsWith(geno_dir, '/'))) {
  geno_dir <- paste(geno_dir, '/', sep = '')
}

# iterate through the Rdata files in the dosages directory, extract and write
#   the three tables within (pruned_dosages, nameList, pruned_pos)
for (i in 1:20){
  str_i <- toString(i)
  if (str_i == '20') {
    str_i <- 'X'
  }
  chr_file <- paste(geno_dir, "chr", str_i,
    ".prunedgen.final.maf001.0.98.RData", sep = '')
  print(paste("Extracting", chr_file))
  chr_data <- load(chr_file)
  write.table(pruned_dosages, row.names = FALSE, sep = "\t",
    file = paste(geno_dir, "chr", str_i, "_pruned_dosages.tsv", sep = ''))
  write.table(nameList, row.names = FALSE, sep = "\t",
    file = paste(geno_dir, "chr", str_i, "_nameList.tsv", sep = ''))
  write.table(pruned_pos, row.names = FALSE, sep = "\t",
    file = paste(geno_dir, "chr", str_i, "_pruned_pos.tsv", sep = ''))
  rm(pruned_dosages, nameList, pruned_pos, chr_data)
}
#
