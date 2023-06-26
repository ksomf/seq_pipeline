library(Biostrings)
library(GenomicRanges)
library(tidyverse)
library(Rsamtools)
library(ggseqlogo)

edit_file       <- snakemake@input[['edit_sites']]
genome_file     <- snakemake@input[['genome']]
gtf_file        <- snakemake@input[['gtf']]

output_file    <- snakemake@output[['logo']]

#metadata <- read_tsv('metadata.tsv')
#bam_files       <- paste0('03_aligned/', filter( metadata, condition=='FL' )$sample_id, '.star_aligned.bam')
#bam_index_files <- paste0('03_aligned/', filter( metadata, condition=='FL' )$sample_id, '.star_aligned.bam.bai')
#edit_file       <- '04_stamp/complex_normal_condition_{name}.tsv'
#genome_file     <- '../reference/ucsc/hg38.fasta'
#gtf_file        <- '../reference/ucsc/hg38.gtf'
#output_file     <- '04_stamp/motif_logo.svg'

genome <- readDNAStringSet(genome_file)
gtf    <- rtracklayer::import(gtf_file)

motif_width <- 20

edit_sites <- read_tsv( edit_file, show_col_types=F ) %>% 
	dplyr::select( chrom, start, end, strand, edit_location, gene_id, gene_name, gene_biotype ) %>% 
	mutate( start=start+1 )
edit_sites_gr <- makeGRangesFromDataFrame( edit_sites, keep.extra.columns=T )

gtf_exons <- subset( gtf, type == 'exon' )
overlaps <- findOverlaps( edit_sites_gr, gtf_exons, type='within' )
exon_hits <- subset( gtf_exons, overlaps@to )
GenomicRanges::elementMetadata(exon_hits)[,'seq'] <- BSgenome::getSeq( genome, exon_hits )
exon_hits

df_exon_hits <- as.data.frame(exon_hits) %>% 
	dplyr::select( seqnames, start, end, strand, gene_id, transcript_id, exon_number, exon_id, seq )
glimpse(df_exon_hits)

df_edit_sites <- edit_sites %>% 
	dplyr::select( start, gene_name, gene_biotype ) %>% 
	rename(position=start)

df_combined <- map2_dfr( overlaps@from, overlaps@to, function(x,y){
	bind_cols( df_edit_sites[x,], df_exon_hits[y,])
}) %>% 
	mutate( seq_offset=position-start+1, seq_width=str_length(seq) ) %>% 
	mutate( slice_start=pmax(seq_offset-motif_width, 1), slice_end=pmin(seq_offset+motif_width, end-start) ) %>% 
	mutate(seq_pattern=substr( seq, slice_start, slice_end )) %>% 
	filter(str_length(seq_pattern) == 2*motif_width+1) %>% 
	filter(str_detect(seq_pattern, 'N', negate=TRUE )) %>% 
	filter(substr( seq_pattern, motif_width+1, motif_width+1 ) == 'C')
ggseqlogo(df_combined$seq_pattern)
ggsave( output_file, width=14, height=9 )

#pileup_params <- PileupParam( distinguish_nucleotides=TRUE, distinguish_strands=FALSE ) 
#p <- pileup( bam_files[1], index=bam_index_files[1], scanBamParam=scan_params, pileupParam=pileup_params ) %>% 
#	group_by(seqnames, which_label, pos) %>% 
#	slice_max( desc(count), n=1, with_ties=F ) %>% 
#	dplyr::select(-count) %>% 
#	ungroup() %>% 
#	group_by(seqnames, which_label) %>% 
#	group_modify(~data.frame(sequence=.x %>% arrange(pos) %>% pull(nucleotide) %>% str_flatten()))
#p
