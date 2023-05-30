library(GenomicRanges)
library(Rsamtools)
library(tidyverse)

source('workflow/scripts/tracks.R')

#simple_bullseye_results  <- snakemake@input[['simple_bullseye']]
#complex_bullseye_results <- snakemake@input[['complex_gbullseye']]
#simple_multi_results     <- snakemake@input[['multi_simple_bullseye']]
#complex_multi_results    <- snakemake@input[['multi_complex_gbullseye']]
#genome_filename          <- snakemake@input[['genome']]
#gtf_filename             <- snakemake@input[['gtf']]
#gtf_database             <- snakemake@params[['gtf_database']]
#read_files               <- snakemake@input[['bams']]
#read_conditions          <- snakemake@params[['conditions']]
#condition_order          <- snakemake@params[['condition_order']]
#output_dir               <- snakemake@output[['plot_dir']]

metadata <- read_tsv('metadata.tsv')
sample2condition <- metadata %>% 
	select( sample_id, condition ) %>% 
	deframe()
simple_bullseye_results  <- list.files( '04_stamp/', full.names=T ) %>% .[str_starts( ., '04_stamp//simple_normal' )  ]  %>% .[str_ends( ., '.tsv' )] %>% .[str_ends( ., 'edited_genes.tsv', negate=T )]
complex_bullseye_results <- list.files( '04_stamp/', full.names=T ) %>% .[str_starts( ., '04_stamp//complex_normal' ) ]  %>% .[str_ends( ., '.tsv' )] %>% .[str_ends( ., 'edited_genes.tsv', negate=T )]
simple_multi_results     <- list.files( '04_stamp/', full.names=T ) %>% .[str_starts( ., '04_stamp//simple_relaxed' )  ] %>% .[str_ends( ., '.tsv' )] %>% .[str_ends( ., 'edited_genes.tsv', negate=F )] 
complex_multi_results    <- list.files( '04_stamp/', full.names=T ) %>% .[str_starts( ., '04_stamp//complex_relaxed' ) ] %>% .[str_ends( ., '.tsv' )] %>% .[str_ends( ., 'edited_genes.tsv', negate=F )]
genome_filename          <- '../reference/ucsc/hg38.fasta'
gtf_filename             <- '../reference/ucsc/hg38.gtf'
gtf_database             <- 'ucsc'
read_files               <- list.files( '04_stamp/', full.names=T ) %>% .[str_ends( ., '.matrix.gz' )]
read_filenames           <- list.files( '04_stamp/', full.names=F ) %>% .[str_ends( ., '.matrix.gz' )]
read_conditions          <- read_filenames %>% lapply(function(x) sample2condition[[str_split(x, '\\.')[[1]][[1]]]]) %>% unlist()
condition_order          <- c( 'FL', 'C-term', 'CTRL' )
output_dir               <- '04_stamp/plots'

message('Reading Genome')
genome <- Rsamtools::FaFile(genome_filename)
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
	set_names(read_files) %>% 
	map_dfr(function(x) read_tsv( gzfile(x), col_names=c('chrom', 'site', 'A', 'T', 'C', 'G', '_', 'N'), show_col_types=F ), .id='file') %>% 
	left_join( set_names(read_files, read_conditions) %>% enframe( name='condition', value='file'), by='file' )

bullseye_result_files <- c( simple_bullseye_results, complex_bullseye_results, simple_multi_results, complex_multi_results )

bullseye_results <- bullseye_result_files %>% 
	set_names(bullseye_result_files) %>% 
	map_dfr( read_tsv, show_col_types=F, .id='file' )

bullseye_genes <- bullseye_results %>% 
	pull(gene_id) %>% 
	unique()

display_gtf <- filter( gtf, gene_id %in% bullseye_genes )

gene_list <- unique(display_gtf$gene_id)

gene <- gene_list[[1]]
for( gene in gene_list ){
	gene_gtf <- filter( display_gtf, gene_id == gene )
	gene_chrom <- unique(gene_gtf$seqnames)
	gene_min_site <- min(gene_gtf$start, gene_gtf$end)
	gene_max_site <- max(gene_gtf$start, gene_gtf$end)
	
	gene_bullseye <- filter( bullseye_results, chrom==gene_chrom & between( start, gene_min_site, gene_max_site ) )
	
	gene_gr <- GRanges( gene_chrom, IRanges( start=gene_min_site, end=gene_max_site ) )
	
	sequence <- getSeq( genome, gene_gr )[[1]] %>% 
		as.character() 
	
	sequence_df <- str_split(sequence, '')[[1]] %>%
		enframe( name='site', value='reference' ) %>% 
		mutate( site=site + gene_min_site - 1 )
	
	gene_sites <- filter( site_reads, chrom==gene_chrom & between( site, gene_min_site, gene_max_site ) )
	gene_sites <- left_join( gene_sites, sequence_df, by='site')
	
	edit_sites <- gene_sites %>% 
		filter(reference == 'C') %>% 
		group_by( condition, site ) %>% 
		summarise( N_tot = sum(N), C_tot=sum(`C`), T_tot=sum(`T`), .groups='drop') %>% 
		mutate( prop_edited = T_tot / N_tot ) %>% 
		filter(T_tot != 0) %>% 
		arrange(desc(prop_edited)) %>% 
		mutate(prop=prop_edited)

	min_site <- min(gene_sites$site)
	max_site <- max(gene_sites$site)
	
	stamp_pileups <- gene_sites %>% 
		select( file, condition, chrom, site, N ) %>% 
		group_by( file, condition, chrom ) %>% 
		group_modify(function(df, g){ #Add missing zeros
			zero_sites <- setdiff( min_site:max_site, df$site )
			rbind( df, data.frame( site=zero_sites, N=0 ) )
		}) %>% 
		ungroup() %>% 
		arrange(site) %>% 
		group_by( condition, chrom, site ) %>% 
		mutate( mean = mean(N), std = sd(N) ) %>% 
		select( condition, chrom, site, mean, std ) %>% 
		ungroup()
	stamp_norm <- ceiling(max(stamp_pileups$mean, stamp_pileups$mean + 2*stamp_pileups$std, na.rm=T))
	
	source('workflow/scripts/tracks.R')
	colour_scheme1 <- khroma::colour('muted')(9)
	colour_map <- c(pileup=colour_scheme1[[2]], prop=colour_scheme1[[1]])
	track <- tracks_create() %>% 
		tracks_annotation(gene_gtf) 
	for( condition_i in condition_order ){
		sites_condition_1 <- filter(stamp_pileups, condition==condition_i)
		edits_condition_1 <- filter(edit_sites   , condition==condition_i) %>% select(site, prop)
		track <- tracks_pileup2( track, sites_condition_1, norm=stamp_norm, prop=edits_condition_1, axis_label=condition_i )
	}	
	p <- tracks_plot(track)
	p <- p + 
		scale_fill_manual  (values=colour_map) +
		scale_colour_manual(values=colour_map)
	
	gene_label <- ifelse(is.null(gene_gtf$label), gene_gtf$gene_id, gene_gtf$label)[[1]]
	save_name <- paste0( output_dir, '/', gene_chrom, '.', gene_min_site, '-', gene_max_site, '_', gene_label, '.svg' )
	dir.create( output_dir, showWarnings=FALSE, recursive=TRUE )
	ggsave( save_name, plot=p, width=24, height=8 )
	write_tsv( gene_bullseye, str_replace( save_name, '.svg', '.tsv' ) )
}

