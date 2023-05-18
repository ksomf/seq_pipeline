# Taken from DEQ, using some of it's helper functions

library(yaml)
library(locfit)
library(tidyverse)
library(deq)

diffbind_file   <- snakemake@input[['diffbind_peaks']]
gtf_filename    <- snakemake@input[['gtf']]
gene2name_file  <- snakemake@input[['gene_conversion']]
output_file     <- snakemake@output[['diffbind_peaks']]

# Get ready gtf
txdb       <- GenomicFeatures::makeTxDbFromGFF( gtf_filename, format='gtf' )
gtf        <- rtracklayer::import(gtf_filename)
gene_names <- unique(as.data.frame(gtf)[,c("gene_id","gene_name","strand")])
anno.order <- c("utr3","utr5","exon","intron","utr3*","utr5*")

# Read in tsv
peaks_df          <- read_tsv(diffbind_file)
peaks_df$old_name <- peaks_df$name
peaks             <- GenomicRanges::makeGRangesFromDataFrame( peaks_df, keep.extra.columns=TRUE )
peaks             <- GenomeInfoDb::keepStandardChromosomes( peaks, pruning.mode='coarse' )
names(peaks)      <- c(paste0( 'peak', 1:length(peaks) ))
peaks.anno        <- deq:::annotateFeatures( peaks, txdb=txdb, genenames=gene_names, utr5.extend=2000, utr3.extend=2000 )
names(peaks.anno) <- c( 'peaks', 'overlaps' )
peaks.anno$peaks  <- deq:::prioritizeAnnotations( peaks.anno$peaks, feature.order=anno.order, genenames=gene_names )

# Extract name and type from the anno
gene_id2gene_name <- read_tsv(gene2name_file)
annotation <- peaks.anno$peaks %>%
	as.data.frame() %>%
	select( old_name, annot, main.gene.id ) %>%
	rename( name=old_name, gene_id=main.gene.id ) %>%
	left_join( gene_id2gene_name, by='gene_id' )

res <- peaks_df %>%
	select(-old_name) %>%
	left_join( annotation, by='name' )
write_tsv( res, output_file )
