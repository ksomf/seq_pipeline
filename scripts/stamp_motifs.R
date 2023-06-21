library(Biostrings)
library(GenomicRanges)
library(tidyverse)
library(Rsamtools)
library(ggseqlogo)

bam_files       <- snakemake@input[['bam']]
bam_index_files <- snakemake@input[['bai']]
edit_file       <- snakemake@input[['edit_sites']]
genome_file     <- snakemake@input[['genome']]
gtf_file        <- snakemake@input[['gtf']]

metadata <- read_tsv('metadata.tsv')
bam_files       <- paste0('03_aligned/', filter( metadata, condition=='FL' )$sample_id, '.star_aligned.bam')
bam_index_files <- paste0('03_aligned/', filter( metadata, condition=='FL' )$sample_id, '.star_aligned.bam.bai')
edit_file       <- '04_stamp/complex_normal_condition_fl_vs_all.tsv'
genome_file     <- '../reference/ucsc/hg38.fasta'
gtf_file        <- '../reference/ucsc/hg38.gtf'

genome <- readDNAStringSet(genome_file)
gtf    <- rtracklayer::import(gtf_file)

motif_width <- 4


edit_sites <- read_tsv( edit_file, show_col_types=F ) %>% 
	select( chrom, start, end, strand, edit_location, gene_id, gene_name, gene_biotype ) %>% 
	mutate( start=start+1 )
edit_sites_gr <- makeGRangesFromDataFrame( edit_sites, keep.extra.columns=T )
	

GenomicRanges::elementMetadata(edit_sites_gr)[,'seq'] <- BSgenome::getSeq( genome, edit_sites_gr + motif_width ) %>% as.character()

ggseqlogo(edit_sites_gr$seq)

scan_params   <- ScanBamParam( which=edit_sites_gr+motif_width, what=scanBamWhat() )
pileup_params <- PileupParam( distinguish_nucleotides=TRUE, distinguish_strands=FALSE ) 
pileup( bam_files[1], index=bam_index_files[1], scanBamParam=scan_params, pileupParam=pileup_params ) %>% 
	group_by(seqnames, which_label, pos) %>% 
	slice_max(desc(count), n=1) %>% 
	select(-count) %>% 
	ungroup() %>% 
	group_by(seqnames, which_label) %>% 
	group_modify(~data.frame(sequence=.x %>% arrange(pos) %>% pull(nucleotide))
	
	
