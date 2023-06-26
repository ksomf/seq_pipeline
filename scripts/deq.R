library(yaml)
library(locfit)
library(tidyverse)
library(deq)

## TEST
#metadata <- read_tsv('metadata.tsv')
#metadata_condition1 <- filter( metadata, condition=='{condition_1}' )
#metadata_condition2 <- filter( metadata, condition=='{condition_2}' )
#metadata_condition1_samples_ip    <- filter( metadata_condition1, method=='IP' )   $sample_id
#metadata_condition1_samples_input <- filter( metadata_condition1, method=='Input' )$sample_id
#metadata_condition2_samples_ip    <- filter( metadata_condition2, method=='IP' )   $sample_id
#metadata_condition2_samples_input <- filter( metadata_condition2, method=='Input' )$sample_id
#
#condition1_ip_filenames    <- paste0( '03_aligned/', metadata_condition1_samples_ip   , '.star_aligned.bam' )
#condition1_input_filenames <- paste0( '03_aligned/', metadata_condition1_samples_input, '.star_aligned.bam' )
#condition2_ip_filenames    <- paste0( '03_aligned/', metadata_condition2_samples_ip   , '.star_aligned.bam' )
#condition2_input_filenames <- paste0( '03_aligned/', metadata_condition2_samples_input, '.star_aligned.bam' )
#peak_filename              <- '04_peakcalling/deq/{condition_1}_vs_{condition_2}_combined_peaks.bed'
#gtf_filename               <- '../reference/ucsc/hg38.gtf'
#
#read_length                <- 75

#output_filename <- '04_peakcalling/deq/{condition_1}_vs_{condition_2}.txt'

# SNAKEMAKE
condition1_ip_filenames    <- snakemake@input[['condition1_ips']]
condition1_input_filenames <- snakemake@input[['condition1_inputs']]
condition2_ip_filenames    <- snakemake@input[['condition2_ips']]
condition2_input_filenames <- snakemake@input[['condition2_inputs']]
peak_filename              <- snakemake@input[['peaks']]
gtf_filename               <- snakemake@input[['gtf']]
bam_infos                  <- snakemake@input[['bam_infos']]

read_length                <- snakemake@params[['read_length']]

threads                   <- snakemake@threads

avg_fragment_length        <- bam_infos %>%
	map_dfr( read_tsv ) %>%
	mutate( reads = `1st fragments`) %>%
	mutate( avg_insert_length   = `insert size average`) %>%
	mutate( avg_fragment_length = 2*read_length + avg_insert_length) %>%
	mutate( total_reads = sum(reads) ) %>%
	mutate( inte = (reads/total_reads) * avg_fragment_length ) %>%
	pull( inte ) %>%
	sum()

output_peaks_filename  <- snakemake@output[['peaks']]
output_counts_filename <- snakemake@output[['counts']]

deq_results <- deq( input.bams         = condition2_input_filenames
                  , ip.bams            = condition2_ip_filenames
                  , treated.input.bams = condition1_input_filenames
                  , treated.ip.bams    = condition1_ip_filenames
                  , peak.files         = peak_filename
                  , gtf                = gtf_filename
                  , paired.end         = TRUE
                  , readlen            = read_length
                  , fraglen            = avg_fragment_length
                  , nthreads           = threads            )

write_tsv( as.data.frame(deq_results$counts) , output_counts_filename )
write_tsv( as.data.frame(deq_results$results), output_peaks_filename )
