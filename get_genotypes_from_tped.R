library(tidyverse)
library(magrittr)

args = commandArgs(trailingOnly=TRUE)

geno = data.table::fread('all_strains.tped', sep = '\t', header = T)
src_strains = colnames(geno)[-c(1:4)]

make_geno = function(strain_fn, dst_stem = str_remove(strain_fn, '.strains.txt')){
  dst_tped = str_glue('{dst_stem}.tped')
  dst_tfam = str_glue('{dst_stem}.tfam')

  targets = data.table::fread(strain_fn, sep = '\t', header = F) %>%
    set_names(c('FID', 'IID'))

  absent = setdiff(targets$FID, src_strains)
  if (length(absent)) {
    message(paste(c('Strains not in source file:', absent), collapse = '\n'))
    targets %<>% filter(FID %in% src_strains)
  }

  # write .tped file
  geno[, c(colnames(geno)[1:4], targets$FID), with = F] %>%
    data.table::fwrite(dst_tped, sep = '\t', quote = F, row.names = F, col.names = F)

  # write .tfam file
  data_frame(FID = targets$FID %>% str_replace_all('[ /]', '.'),
             IID = targets$IID %>% str_replace_all('[ /]', '.'),
             fa = 0,
             ma = 0,
             sex = 0,
             pheno = -9) %>%
    write.table(dst_tfam, sep = '\t', quote = F, row.names = F, col.names = F)
}

# example
#make_geno('bile_acid_plasma_lipidome.strains.txt')
#make_geno('fid_iid_chow1.txt')
make_geno(args[1])
#
