library(GenomicRanges)
library(ggplot2)
library(ggthemes)
library(tidyverse)

pepr_peaks_filename <- '04_peakcalling/pepr/MAVSvsd103-467__PePr_chip1_peaks.bed'
gff_filename        <- '../reference/hg38.gff'
output_dir          <- '04_peakcalling/analysis/plots/'

#pepr_peaks_filename <- snakemake@input[['pepr_peaks']]
#gff_filename        <- snakemake@input[['gff']]
#output_dir          <- snakemake@output[['plot_dir']]

broadpeak_colnames <- c( 'chrom', 'start', 'end', 'name', 'score', 'strand', 'enrichment', 'p_value', 'q_value' )
pepr_peaks <- read_tsv( pepr_peaks_filename, col_names=broadpeak_colnames, show_col_types=FALSE ) %>%
	mutate( strand=ifelse(strand == '.', '*', strand) )
pepr_peaks <- GRanges( seqnames=pepr_peaks$chrom
                     , ranges=IRanges(pepr_peaks$start, end=pepr_peaks$end, names=pepr_peaks$name)
                     , strand=pepr_peaks$strand
                     , p_value=pepr_peaks$p_value
                     , method='pepr' )

gff_columns <- c( 'chrom', 'label', 'level', 'start', 'end', 'score', 'strand', 'something', 'metadata' )
gff <- read_tsv( gff_filename, col_names=gff_columns, comment='#', show_col_types=FALSE ) %>%
	filter( str_starts( metadata, 'ID=gene:' ) ) %>%
	mutate( metadata=str_match( x, 'ID=gene:([^;]+);' )[,2] ) #TODO continue here make this line work
print(gff)
