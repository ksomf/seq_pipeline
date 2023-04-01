library(tidyverse)
library(deq)


## TEST
#metadata <- read_tsv('metadata.tsv')
#metadata_condition1 <- filter( metadata, condition=='MAVS' )
#metadata_condition2 <- filter( metadata, condition=='d103-467' )
#metadata_condition1_samples_ip    <- filter( metadata_condition1, method=='IP' )   $sample_id
#metadata_condition1_samples_input <- filter( metadata_condition1, method=='Input' )$sample_id
#metadata_condition2_samples_ip    <- filter( metadata_condition2, method=='IP' )   $sample_id
#metadata_condition2_samples_input <- filter( metadata_condition2, method=='Input' )$sample_id
#
#condition1_ip_filenames    <- paste0( '03_aligned/', metadata_condition1_samples_ip   , '_star_aligned.bam' )
#condition1_input_filenames <- paste0( '03_aligned/', metadata_condition1_samples_input, '_star_aligned.bam' )
#condition2_ip_filenames    <- paste0( '03_aligned/', metadata_condition2_samples_ip   , '_star_aligned.bam' )
#condition2_input_filenames <- paste0( '03_aligned/', metadata_condition2_samples_input, '_star_aligned.bam' )
#peak_filenames             <- paste0( '04_peakcalling/', c(metadata_condition1_samples_ip, metadata_condition2_samples_ip), '_normal_peaks.narrowPeak' )
#gtf_filename               <- '../reference/hg38.gtf'
#
#output_filename <- '../04_peakcalling/deq/MAVS_vs_d103-467.txt'

# SNAKEMAKE
condition1_ip_filenames    <- snakemake@input[['condition1_ips']]
condition1_input_filenames <- snakemake@input[['condition1_inputs']]
condition2_ip_filenames    <- snakemake@input[['condition2_ips']]
condition2_input_filenames <- snakemake@input[['condition2_inputs']]
peak_filenames             <- snakemake@input[['peaks']]
gtf_filename               <- snakemake@input[['gtf']]

output_filename <- snakemake@output[['peaks']]

deq( input.bams         = condition2_input_filenames
   , ip.bams            = condition2_input_filenames
   , treated.input.bams = condition1_input_filenames
   , treated.ip.bams    = condition1_input_filenames
   , peak.files         = peak_filenames
   , gtf                = gtf_filename
   , paired.end         = TRUE
   , outfi              = output_filename            )
