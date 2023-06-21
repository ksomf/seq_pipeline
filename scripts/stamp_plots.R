library(GenomicRanges)
library(Rsamtools)
library(tidyverse)

source('workflow/scripts/tracks.R')

simple_bullseye_results  <- snakemake@input[['simple_bullseye']]
complex_bullseye_results <- snakemake@input[['complex_bullseye']]
simple_multi_results     <- snakemake@input[['multi_simple_bullseye']]
complex_multi_results    <- snakemake@input[['multi_complex_bullseye']]
genome_filename          <- snakemake@input[['genome']]
gtf_filename             <- snakemake@input[['gtf']]
gtf_database             <- snakemake@params[['gtf_database']]
read_files               <- snakemake@input[['counts']]
read_conditions          <- snakemake@params[['conditions']]
condition_order          <- snakemake@params[['condition_order']]
output_dir               <- snakemake@output[['plot_dir']]
extra_bam_files          <- snakemake@params[['extra_bam_files']]

#metadata <- read_tsv('metadata.tsv')
#sample2condition <- metadata %>% 
#	select( sample_id, condition ) %>% 
#	deframe()
#simple_bullseye_results  <- list.files( '04_stamp/', full.names=T ) %>% .[str_starts( ., '04_stamp//simple_normal' )  ]  %>% .[str_ends( ., '.tsv' )] %>% .[str_ends( ., 'edited_genes.tsv', negate=T )]
#complex_bullseye_results <- list.files( '04_stamp/', full.names=T ) %>% .[str_starts( ., '04_stamp//complex_normal' ) ]  %>% .[str_ends( ., '.tsv' )] %>% .[str_ends( ., 'edited_genes.tsv', negate=T )]
#simple_multi_results     <- list.files( '04_stamp/', full.names=T ) %>% .[str_starts( ., '04_stamp//simple_relaxed' )  ] %>% .[str_ends( ., '.tsv' )] %>% .[str_ends( ., 'edited_genes.tsv', negate=F )] 
#complex_multi_results    <- list.files( '04_stamp/', full.names=T ) %>% .[str_starts( ., '04_stamp//complex_relaxed' ) ] %>% .[str_ends( ., '.tsv' )] %>% .[str_ends( ., 'edited_genes.tsv', negate=F )]
#genome_filename          <- '../reference/ucsc/hg38.fasta'
#gtf_filename             <- '../reference/ucsc/hg38.gtf'
#gtf_database             <- 'ucsc'
#read_files               <- list.files( '04_stamp/', full.names=T ) %>% .[str_ends( ., '.matrix.gz' )]
#read_filenames           <- list.files( '04_stamp/', full.names=F ) %>% .[str_ends( ., '.matrix.gz' )]
#read_conditions          <- read_filenames %>% lapply(function(x) sample2condition[[str_split(x, '\\.')[[1]][[1]]]]) %>% unlist()
#condition_order          <- c( 'FL', 'C-term', 'CTRL' )
#output_dir               <- '04_stamp/plots'
#ripseq_location <- '../nandan_mavs_ripseq/'
#extra_bam_files <- read_tsv(paste0(ripseq_location,'metadata.tsv')) %>% 
#	filter(method == 'IP') %>% 
#	nest(.by='condition') %>% 
#	mutate( files = lapply( data, function(df){
#		 list( ip    = paste0(ripseq_location, '/03_aligned/', df$sample_id             , '.star_aligned.bam')
#		     , input = paste0(ripseq_location, '/03_aligned/', df$matching_input_control, '.star_aligned.bam') )
#	} ) ) %>% 
#	select(-data) %>% 
#	deframe()

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
bullseye_result_files <- c( simple_bullseye_results, complex_bullseye_results, simple_multi_results, complex_multi_results )
bullseye_results <- bullseye_result_files %>% 
	set_names(bullseye_result_files) %>% 
	map_dfr( read_tsv, show_col_types=F, .id='file' )
bullseye_genes <- bullseye_results %>% 
	pull(gene_id) %>% 
	unique()

display_gtf <- filter( gtf, gene_id %in% bullseye_genes )
gene_widths <- display_gtf %>%
	group_by(gene_id) %>% 
	summarise( start=min(start), end=max(end) ) %>% 
	mutate(width=end-start) %>% 
	arrange(width)
gene_list <- gene_widths$gene_id

sites_of_interest <- display_gtf %>%
	makeGRangesFromDataFrame() %>% 
	GenomicRanges::reduce() %>% 
	as.data.frame() %>% 
	rename(chrom=seqnames) %>% 
	group_by(chrom) %>% 
	group_modify(function(df, g){
		data.frame(sites=map2( df$start, df$end, ~.x:.y ) %>% unlist())
	}) %>% 
	summarise(sites=list(sites)) %>% 
	deframe()

site_reads <- read_files %>% 
	set_names(read_files) %>% 
	map_dfr(function(x){
		print(x)
		read_tsv( gzfile(x), col_names=c('chrom', 'site', 'A', 'T', 'C', 'G', '_', 'N'), show_col_types=F ) %>% 
			filter(chrom %in% names(sites_of_interest)) %>% 
			group_by(chrom) %>% 
			group_modify(function(df, g){
				#print(g$chrom)
				#print(df)
				#print(sites_of_interest[[g$chrom]])
				filter(df, site %in% sites_of_interest[[g$chrom]])
				})
		}, .id='file') %>% 
	ungroup() %>% 
	left_join( set_names(read_files, read_conditions) %>% enframe( name='condition', value='file'), by='file' )


gene <- gene_list[[1]]
for( i in 1:length(gene_list) ){
	gene <- gene_list[[i]]
	gene_gtf <- filter( display_gtf, gene_id == gene )
	gene_chrom <- unique(gene_gtf$seqnames)
	gene_min_site <- min(gene_gtf$start, gene_gtf$end)
	gene_max_site <- max(gene_gtf$start, gene_gtf$end)
	gene_width    <- gene_max_site - gene_min_site
	print(paste0(gene, '(', gene_width, '): ', i, '/', length(gene_list)))
	
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
	edit_norm <- max( edit_sites$prop, na.rm=T )

	min_site <- min(gene_sites$site)
	max_site <- max(gene_sites$site)
	
	stamp_pileups <- gene_sites %>% 
		select( file, condition, chrom, site, N ) %>% 
		group_by( file, condition, chrom ) %>% 
		group_modify(function(df, g){ #Add missing zeros
			zero_sites <- setdiff( min_site:max_site, df$site )
			if(length(zero_sites) != 0){
				df <- rbind( df, data.frame( site=zero_sites, N=0 ) )
			}
			df
		}) %>% 
		ungroup() %>% 
		arrange(site) %>% 
		group_by( condition, chrom, site ) %>% 
		mutate( min=min(N), max=max(N), median=median(N), q1=quantile(N, .25), q3=quantile(N, .75), mean=mean(N), std=sd(N) ) %>% 
		select( condition, chrom, site, min, max, median, q1, q3 ) %>% 
		ungroup()
	stamp_norm <- c(pileup=ceiling(max(stamp_pileups$max, na.rm=T)))
	
	scan_params   <- ScanBamParam( which=gene_gr, what=scanBamWhat() )
	pileup_params <- PileupParam( distinguish_nucleotides=FALSE, distinguish_strands=FALSE )
	
	extra_pileups <- extra_bam_files %>% 
		enframe( name='condition', value='data' ) %>% 
		mutate( data=lapply( data, data.frame ) ) %>% 
		mutate( bams=lapply( data, function(df){
			df %>% 
				mutate(across( everything(), function(col) lapply( col, function(f) pileup( f, index=str_replace(f, '.bam', '.bam.bai'), scanBamParam=scan_params, pileupParam=pileup_params ) ) ))
		}) ) %>% 
		select(-data) %>% 
		unnest(bams) %>% 
		pivot_longer(cols=c('ip', 'input'), names_to='type', values_to='pileup') %>% 
		mutate( pileup=lapply( pileup, function(df){
			zero_sites <- setdiff( min_site:max_site, df$pos )
			df <- df %>% 
				select( -seqnames, -which_label ) %>% 
				rename( site='pos', N='count' ) 
			if(!is.null(nrow(zero_sites))){
				df <- df %>% 
					rbind(data.frame( site=zero_sites, N=0 )) %>% 
					filter(between( site, min_site, max_site ))
			}
			df
		})) %>% 
		unnest(pileup) %>% 
		group_by( condition, type, site ) %>% 
		summarise( min=min(N), max=max(N), median=median(N), q1=quantile(N, .25), q3=quantile(N, .75), mean=mean(N), std=sd(N), .groups='drop' )
	extra_types <- unique(extra_pileups$type)
	extra_norm <-  extra_pileups %>% 
		group_by(type) %>% 
		summarise(norm=ceiling(max(max, na.rm=T))) %>% 
		deframe()
	
	#source('workflow/scripts/tracks.R')
	colour_scheme1 <- khroma::colour('muted')(9)
	colour_map <- c(pileup=colour_scheme1[[2]], prop=colour_scheme1[[1]], input=colour_scheme1[[3]], ip=colour_scheme1[[4]])
	track <- tracks_create()
	
	for( condition_j in unique(extra_pileups$condition) ){
		sites_condition_j <- filter( extra_pileups, condition==condition_j )
		for( t in extra_types ){
			track <- tracks_pileup_shade( obj=track, df=filter( sites_condition_j, type==t ), separation_variable='type', norm=extra_norm[t], colour_is_group=T, axis_label=condition_j)
		}
	}
	
	track <- tracks_annotation( track, gene_gtf ) 
	for( condition_i in condition_order ){
		sites_condition_i <- filter(stamp_pileups, condition==condition_i)
		edits_condition_i <- filter(edit_sites   , condition==condition_i) %>% select(site, prop)
		track <- tracks_pileup_shade( track, sites_condition_i, separation_variable='pileup', norm=stamp_norm, prop=edits_condition_i, prop_norm=edit_norm, axis_label=condition_i )
	}	
	p <- tracks_plot(track)
	p <- p + 
		scale_fill_manual  (values=colour_map) +
		scale_colour_manual(values=colour_map)
	p
		
	gene_label <- ifelse(is.null(gene_gtf$label), gene_gtf$gene_id, gene_gtf$label)[[1]]
	save_name <- paste0( output_dir, '/', gene_chrom, '.', gene_min_site, '-', gene_max_site, '_', gene_label, '.svg' )
	dir.create( output_dir, showWarnings=FALSE, recursive=TRUE )
	ggsave( save_name, plot=p, width=24, height=8 )
	ggsave( str_replace(save_name, '.svg', 'png'), plot=p, width=24, height=8 )
	write_tsv( gene_bullseye, str_replace( save_name, '.svg', '.tsv' ) )
}

