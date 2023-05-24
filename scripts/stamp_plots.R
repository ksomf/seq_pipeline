library(GenomicRanges)
library(Rsamtools)
library(tidyverse)


#simple_bullseye_results  <- snakemake@input[['simple_bullseye']]
#complex_bullseye_results <- snakemake@input[['complex_gbullseye']]
#simple_multi_results     <- snakemake@input[['multi_simple_bullseye']]
#complex_multi_results    <- snakemake@input[['multi_complex_gbullseye']]
#gtf_filename             <- snakemake@input[['gtf']]
#gtf_database             <- snakemake@params[['gtf_database']]
#read_files               <- snakemake@input[['bams']]
#read_conditions          <- snakemake@params[['conditions']]
#condition_order          <- snakemake@params[['condition_order']]
#ouptup_dir               <- snakemake@output[['plot_dir']]

metadata <- read_tsv('metadata.tsv')
sample2condition <- metadata %>% 
	select( sample_id, condition ) %>% 
	deframe()
simple_bullseye_results  <- list.files( '04_stamp/', full.names=T ) %>% .[str_starts( ., '04_stamp//simple_normal' )  ]  %>% .[str_ends( ., '.tsv' )] %>% .[str_ends( ., 'edited_genes.tsv', negate=T )]
complex_bullseye_results <- list.files( '04_stamp/', full.names=T ) %>% .[str_starts( ., '04_stamp//complex_normal' ) ]  %>% .[str_ends( ., '.tsv' )] %>% .[str_ends( ., 'edited_genes.tsv', negate=T )]
simple_multi_results     <- list.files( '04_stamp/', full.names=T ) %>% .[str_starts( ., '04_stamp//simple_relaxed' )  ] %>% .[str_ends( ., '.tsv' )] %>% .[str_ends( ., 'edited_genes.tsv', negate=F )] 
complex_multi_results    <- list.files( '04_stamp/', full.names=T ) %>% .[str_starts( ., '04_stamp//complex_relaxed' ) ] %>% .[str_ends( ., '.tsv' )] %>% .[str_ends( ., 'edited_genes.tsv', negate=F )]
gtf_filename             <- '../reference/ucsc/hg38.gtf'
gtf_database             <- 'ucsc'
read_files               <- list.files( '04_stamp/', full.names=T ) %>% .[str_ends( ., '.matrix.gz' )]
read_conditions          <- bam_files %>% lapply(function(x) sample2condition[[str_split(x, '\\.')[[1]][[1]]]]) %>% unlist()
condition_order          <- c( 'FL', 'C-term', 'CTRL' )
ouptup_dir               <- '04_stamp/plots'

message('Reading GTF')
gtf <- rtracklayer::import(gtf_filename)
if( gtf_database == 'ucsc' ){
	mart <- biomaRt::useEnsembl('ensembl', 'hsapiens_gene_ensembl')
	missing_metadata <- biomaRt::getBM( attributes=c('ensembl_gene_id', 'external_gene_name', 'gene_biotype'), mart=mart, verbose=TRUE ) %>% 
		rename(gene_name='external_gene_name', gene_id='ensembl_gene_id')
	
	gtf_genes <- gtf %>% 
		subset( type=='transcript' ) %>% 
		as.data.frame() %>% 
		group_by(gene_id) %>% 
		summarise( seqnames=seqnames[[1]], start=min(start), end=max(end), strand=strand[[1]] ) %>% 
		left_join( missing_metadata, by='gene_id') %>% 
		mutate( label=ifelse(gene_name != '', gene_name, gene_id) ) %>% 
		mutate( width=end - start, type='gene' )
	
	transcripts2gene  <- gtf %>% 
		subset( type=='transcript') %>% 
		as.data.frame() %>% 
		select(gene_id, transcript_id) %>% 
		distinct() %>% 
		group_by(gene_id) %>% 
		mutate( num_siblings=n(), cannonical=seq_len(n()) == 1 )  %>% 
		left_join(gtf_genes, select(any_of('gene_id', 'gene_biotype')), by='gene_id') %>% 
		select(-seqnames, -start, -end, -width, -strand)
	
	gtf_exons <- gtf %>%
		subset( type=='exon' ) %>% 
		as.data.frame() %>% 
		select(any_of(c('seqnames', 'start', 'end', 'width', 'strand', 'type', 'gene_id', 'transcript_id', 'constitutive', 'exon_id'))) %>% 
		left_join( transcripts2gene, by=c('gene_id', 'transcript_id') ) %>% 
		mutate( type='exon' )
	
	gtf_genes <- mutate( gtf_genes, num_siblings=1, cannonical=TRUE )
	
	common_gtf_columns <- c('seqnames', 'start', 'end', 'width', 'strand', 'type', 'gene_id', 'label', 'gene_biotype', 'num_siblings', 'cannonical')
	gtf <- rbind( select(gtf_genes, any_of(common_gtf_columns)), select(gtf_exons, any_of(common_gtf_columns)) )
}else{
	gtf_genes <- gtf %>%
		subset( type=='gene' ) %>% 
		as.data.frame() %>% 
		select(any_of(c('seqnames', 'start', 'end', 'width', 'strand', 'type', 'gene_id', 'Name', 'gene_biotype', 'biotype', 'description', 'version'))) %>% 
		mutate( num_siblings=1, cannonical=TRUE ) %>% 
		rename( any_of(c(gene_biotype='biotype')) ) %>% 
		mutate( label = ifelse(is.na(Name), gene_id, Name))
	
	transcripts2gene  <- gtf %>% 
		subset( type=='transcript') %>% 
		as.data.frame() %>% 
		select(gene_id, transcript_id) %>% 
		distinct() %>% 
		group_by(gene_id) %>% 
		mutate( num_siblings=n(), cannonical=seq_len(n()) == 1 ) %>% 
		left_join(select(any_of('gtf_genes', 'gene_id', 'gene_biotype', 'gene_biotype')), by='gene_id')
	gtf_exons <- gtf %>%
		subset( type=='exon' ) %>% 
		as.data.frame() %>% 
		select(any_of(c('seqnames', 'start', 'end', 'width', 'strand', 'type', 'gene_id', 'transcript_id', 'constitutive', 'exon_id'))) %>% 
		left_join( transcripts2gene, by=c('gene_id', 'transcript_id') )
	common_gtf_columns <- c('seqnames', 'start', 'end', 'width', 'strand', 'type', 'gene_id', 'gene_biotype', 'num_siblings', 'cannonical')
	gtf <- rbind( select(gtf_genes, any_of(common_gtf_columns)), select(gtf_exons, any_of(common_gtf_columns)) )
}

message('Reading Bullseye Data')
site_reads <- read_files %>% 
	set_names(read_conditions) %>% 
	map_dfr(function(x) read_tsv( gzfile(x), col_names=c('chrom', 'site', 'A', 'T', 'C', 'G', '_', 'N'), show_col_types=F ), .id='condition')

bullseye_genes <- c( simple_bullseye_results, complex_bullseye_results, simple_multi_results, complex_multi_results ) %>% 
	map_dfr( read_tsv, show_col_types=F ) %>% 
	pull(gene_id) %>% 
	unique()

display_gtf <- filter( gtf, gene_id %in% bullseye_genes )

gene_list <- unique(display_gtf$gene_id)

for( gene in gene_list ){
	gene_gtf <- filter( display_gtf, gene_id == gene )
	gene_chrom <- unique(gene_gtf$seqnames)
	gene_min_site <- min(gene_gtf$start, gene_gtf$end)
	gene_max_site <- max(gene_gtf$start, gene_gtf$end)
	gene_sites <- filter( site_reads, chrom==gene_chrom & between( site, gene_min_site, gene_max_site ) )
}
