#add three types of tags to vcf
library(purrr)
library(dplyr)
library(jsonlite)
library(stringr)
library(tidyr)
library(data.table)
library(GenomicAlignments)
library(Rsamtools)
library(optigrab)
runbcftools <- function(BCFoptions = "") system(paste('bcftools',BCFoptions))
runtabix <- function(TABIXoptions = "") system(paste('tabix',TABIXoptions))
runbgzip <- function(BGZIPoptions = "") system(paste('bgzip',BGZIPoptions))

replace_with_third_highest <- function(pairs_list) {
  lapply(pairs_list, function(x) {
    if (length(x) <= 2){
      return(0)
    } else {
      return(sort(as.numeric(x), decreasing = TRUE)[3])
    }
  })
}

# does not handle multistrs but could if needed be changed
calc_NTSR_lc <- function(sample, multi_str){
  json_file <- paste0(sample, ".json")
  json <- read_json(json_file)$LocusResults
  json <- json[!names(json) %in% multi_str]
  count_of_spanning_reads <- map_chr(json, ~ .x$Variants[[1]]$CountsOfSpanningReads)
  pairs_list <- str_extract_all(count_of_spanning_reads, "\\(\\d+|\\d+\\)")
  pairs_list <- lapply(pairs_list, function(x) gsub("\\(|\\)", "", x))
  pairs_list <- lapply(pairs_list, function(x) if(length(x)==0){0} else {x[seq(2, length(x), by=2)]})
  new_list <- replace_with_third_highest(pairs_list)
  
  coverage <- map_chr(json, ~ .x$Coverage)
  av_cov <- sum(as.numeric(coverage))/length(coverage)
  norm_locus_cov <- round(as.numeric(coverage)/av_cov, 3)
  df <- data.frame(locus = unlist(names(count_of_spanning_reads)), NTSR = unlist(new_list), NLC = norm_locus_cov)

  return(df)
}

# Calculate AAP for each locus
calc_aap_vect <- function(temp_file){
  bam_df <- read.table(temp_file, col.names=c("locus", "cigar"))
  loci <- unique(bam_df$locus) # Get loci names
  total_cigars <- tapply(bam_df$cigar, bam_df$locus, paste, sep='', collapse = '') # collapse all cigar strings together for each locus
  aaps <- sapply(total_cigars, function(x) round(cigarOpTable(x)[1]/sum(cigarOpTable(x)), 4)) # find percentage of bases matching the reference
  df <- data.frame(aap=unname(aaps), locus=names(aaps))
  return(df)
}

calc_MADSP <- function(sample){
  vcf <- paste0(sample, ".vcf.gz")
  runbcftools(paste0("query -f '%VARID\t[%ADSP]\n' ", vcf, " > ", vcf, ".txt"))
  file_path <- paste0(vcf, ".txt")
  data <- read.table(file_path, header=FALSE, sep="\t")
  new_names <- c('locus', 'ADSP')
  colnames(data) <- new_names
  data$MADSP <- sapply(strsplit(as.character(data$ADSP), split="/"), function(x) min(as.integer(x)))
  data$MADSP[is.na(data$MADSP)] <- 0
  data <- data %>%
    select(c(locus,MADSP))
  return(data)
}

add_all_tags <- function(sample, aap, df_extra_tags, MADSP){
  # Merge the first two data frames
  temp_merge <- merge(aap, df_extra_tags, by='locus')
  
  # Merge the result with the third data frame
  df_tags <- merge(temp_merge, MADSP, by='locus')
  
  df <- merge(read.table(opt_get('variants'), col.names= c('locus', 'CHROM', 'POS')), 
              df_tags, by='locus', all.x = T) %>% 
    select(CHROM, POS, aap, NTSR, NLC, MADSP) %>%
    replace(is.na(.), '.')
  # write bcf annotation file
  write.table(df, 'temp.annot.tab', sep='\t', quote=F, row.names = F, col.names = F)
  # Sort zip and index annotation file
  system('sort -V -k1,1 -k2,2 temp.annot.tab > temp.sorted_annot.tab')
  runbgzip('temp.sorted_annot.tab')
  runtabix('-s1 -b2 -e2 temp.sorted_annot.tab.gz')
  # Run bcftools to annotate the vcf with values
  runbcftools(paste0('annotate -a temp.sorted_annot.tab.gz -h annot.hdr -c CHROM,POS,FMT/AAP,FMT/NTSR,FMT/NLC,FMT/MADSP ', sample, '.vcf.gz -Oz -o ', sample, '.annot.vcf.gz'))
  runtabix(paste0(sample, '.annot.vcf.gz'))
  system("rm temp*")
}

# Folder contains sorted/index bamlets in genotypes subfolder with .vcf.gz and index
sample <- opt_get('sample')
multi_str <- opt_get('multi-strs')
# Write header lines for 3 annotations and store in inputs directory
write('##FORMAT=<ID=AAP,Number=1,Type=Float,Description="Average Absolute Purity">', 'annot.hdr')
write('##FORMAT=<ID=NTSR,Number=1,Type=Float,Description="No. Third Species Reads">', 'annot.hdr', append = T)
write('##FORMAT=<ID=NLC,Number=1,Type=Float,Description="Normalised Locus Coverage">', 'annot.hdr', append = T)
write('##FORMAT=<ID=MADSP,Number=1,Type=Integer,Description="Minimum value of ADSP">', 'annot.hdr', append = T)
# This script takes the reads that overlap a locus with only one repeat 
# and combines the CIGAR strings for all the reads overlapping that locus, saving the output in the temp directory
system(paste('sh aap.sh', sample, multi_str))

AAP <- calc_aap_vect('temp_data.txt')
extra_tags <- calc_NTSR_lc(sample, multi_str)
MADSP <- calc_MADSP(sample)
add_all_tags(sample, AAP, extra_tags, MADSP)



