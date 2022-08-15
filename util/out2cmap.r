#!/data/wangjiaxuan/biosoft/miniconda3/bin/Rscript

#输出到cmap中--------
require(tidyverse)

bug_all = read_tsv("3.Result_Sum/all.sample_buglist.tsv" )
bug_all %>% 
  filter(feature != "#SampleID") %>%
  select(-feature,feature) %>%
  mutate(feature = gsub("\\|",";",feature)) %>%
  add_column(Taxmony = paste0("Sp",seq(1:nrow(.))),.before = 1) -> bug_all_fix

bug_all_fix[is.na(bug_all_fix)] = "0"
write_tsv(bug_all_fix,"4.Out2CAMP/bugs_CMAP_input.txt")


# 基因功能分析
gene_cpm = read_tsv("3.Result_Sum/all.sample_Functionfamilie_cpms.tsv")

# 输出gene表格
colnames(gene_cpm)[1] = "Function"
gene_cpm %>%  
  select(-Function,Function) %>%
  #mutate(Function = gsub("\\|",";",feature)) %>%
  add_column(Taxmony = paste0("Fc",seq(1:nrow(.))),.before = 1) -> gene_cpm_fix
write_tsv(bug_all_fix,"4.Out2CAMP/bugs_CMAP_input.txt")

# #  path_ab
# path_ab = read_tsv("all.sample_pathabundance.tsv")

# # path_cov
# path_cov = read_tsv("all.sample_pathcoverage.tsv")
