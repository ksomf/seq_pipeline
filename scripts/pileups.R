library(GenomicRanges)
library(Rsamtools)
library(ggplot2)
library(ggthemes)
library(furrr)
library(rtracklayer)
library(tidyverse)

LERP <- function(y1, y2, x){y1 + (y2-y1)*x}

condition_peaks <- '04_peakcalling/analysis/condition_peaks.tsv'
diffbind_peaks  <- '04_peakcalling/analysis/diffbind_peaks.tsv'
bam_files       <- list.files('03_aligned/', full.names=TRUE ) %>% .[str_ends(.,'star_aligned.bam')]
bam_index_files <- list.files('03_aligned/', full.names=TRUE ) %>% .[str_ends(.,'star_aligned.bam.bai')]
library_sizes   <- '04_peakcalling/analysis/library_sizes.tsv'
gtf_filename    <- '../reference/ucsc/hg38.gtf'

sample_ids            <- list.files('03_aligned/', full.names=FALSE) %>% .[str_ends(.,'star_aligned.bam')] %>% str_remove('.star_aligned.bam')
treatment_conditions  <- 'MAVS'
control_condition     <- 'd103-467'
metadata_filename     <- 'metadata.tsv'
gtf_database          <- 'ucsc'

output_dir            <- '04_peakcalling/analysis/plots/'
output_signal         <- '04_peakcalling/analysis/summary.txt'

threads               <- 32

#condition_peaks <- snakemake@input[['condition_peaks']]
#diffbind_peaks  <- snakemake@input[['diffbind_peaks' ]]
#bam_files       <- snakemake@input[['bam_files'      ]]
#bam_index_files <- snakemake@input[['bam_index_files']]
#library_sizes   <- snakemake@input[['library_sizes'  ]]
#gtf_filename    <- snakemake@input[['gtf'            ]]
#
#treatment_conditions  <- snakemake@params[['treatment_conditions']]
#control_condition     <- snakemake@params[['control_condition'   ]]
#metadata_filename     <- snakemake@params[['metadata'            ]]
#sample_ids            <- snakemake@params[['sample_ids'          ]]
#gtf_database          <- snakemake@params[['database'            ]]
#
#output_dir              <- snakemake@output[['plot_dir']]
#output_signal           <- snakemake@output[['summary_file']]
#
#threads                 <- snakemake@threads

write_tsv(data.frame(), output_signal)

plan(multisession, workers=threads)
#TODO normalise bam files against each other using RPKM

print('Reading Metadata')
ordered_conditions <- c(treatment_conditions, control_condition)
metadata <- read_tsv( metadata_filename, show_col_types=FALSE )
sample2condition <- select(metadata, sample_id, condition) %>% deframe()
sample2input     <- select(metadata, sample_id, method   ) %>% mutate(method=method=='Input') %>% deframe()
condition2count  <- table(metadata$condition) / 2 #NOTE(KIM): Assuming equal input and ip

sample2libsize_factor <- read_tsv(library_sizes) %>% deframe()
sample2libsize_factor <- sample2libsize_factor / mean(sample2libsize_factor)

conds      = unique(sample2condition)
num_conditions = length(conds)
condition_offset = set_names((seq_along(conds) - 1)*-3, conds)

print('Reading GTF')
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
	
print('Reading Diffbind Results')
de_peaks        <- read_tsv(diffbind_peaks, show_col_types=FALSE) %>% mutate(siblings=1) %>% rename(any_of(c(chrom='seqnames')))
de_peaks_gr     <- makeGRangesFromDataFrame( de_peaks, keep.extra.columns=TRUE )
de_search_peaks <- de_peaks

print('Reading vs Input Peakcaller Results')
cond_peaks <- read_tsv(condition_peaks, show_col_types=FALSE) %>% rename(any_of(c(chrom='seqnames')))
cond_search_peaks <- cond_peaks %>%
	group_by(method) %>%
	group_modify(function(df, g){
		df_treatment <- subset( df, condition %in% treatment_conditions ) %>% makeGRangesFromDataFrame(keep.extra.columns=TRUE) 
		df_control   <- subset( df, condition %in% control_condition    ) %>% makeGRangesFromDataFrame(keep.extra.columns=TRUE)
		non_overlapping_ranges <- findOverlaps( df_treatment, df_control, select='first' ) %>% is.na()
		df_treatment %>% 
			subset(non_overlapping_ranges) %>% 
			subset(significant) %>% 
			as.data.frame()
	}) %>% 
	rename(any_of(c(chrom='seqnames')))

search_peaks <- rbind( de_search_peaks, select(cond_search_peaks, -width) ) %>% 
	group_by(method) %>% 
	mutate(method_rank=seq_len(n()))
search_peaks <- search_peaks %>% 
	subset( method == 'deq' ) %>% 
	subset(significant)
peaks        <- rbind( de_peaks       ,        cond_peaks                 )

print('Loading simple track tools')
method_offset <- function(xs){
	res <- c()
	if( length(xs) > 0 ){
		x2d <- xs %>% 
			unique() %>% 
			set_names((seq_along(.)-1)/length(.),.)
		res <- map_dbl(xs, ~x2d[[.x]])
}
res
}

tracks_create <- function(padding=0.1){
	list(width=0, tracks=list(), track_guides=c(), padding=padding)
}

tracks_shade_bar <- function( obj=tracks_create(), df, group_column, direction, range ){
	tab_height <- 0.1
	groups     <- unique(df[,group_column])
	group2num  <- set_names(seq_along(groups), groups)
	num_groups <- length(groups)
	shade_width <- tab_height*num_groups
	
	shade_index <- length(obj$tracks) + 1
	base_y <- -obj$width - shade_width
	  
	df_plot <- df %>% 
		mutate( group=df[[group_column]] ) %>% 
		mutate( tab_top=map_dbl(group, ~base_y+tab_height*(group2num[[.x]]  ))
		      , tab_bot=map_dbl(group, ~base_y+tab_height*(group2num[[.x]]-1)) ) %>% 
		mutate( xmin=start, xmax=end, alpha=1.0/siblings ) %>%
		mutate( alpha=ifelse(significant, alpha, alpha*0.25) ) %>% 
		select( xmin, xmax, tab_top, tab_bot, group, alpha )
	
	obj$width <- obj$width + shade_width + obj$padding
	obj$tracks <- append( obj$tracks, list(list( type='shade', df=df_plot, y_bot=min(df_plot$tab_bot), y_top=max(df_plot$tab_top), direction=direction, start_index=shade_index, range=range )) )
	obj
}

tracks_pileup <- function( obj=tracks_create(), df, condition_column ){
	base_y <- -obj$width - 1
	plot_df <- df %>% 
		mutate( xmin=pos-0.5, xmax=pos+0.5, ymin=base_y, ymax=base_y+normcount, condition=df[,condition_column], alpha=alpha ) %>% 
		select( xmin, xmax, ymin, ymax, condition, alpha )
	obj$width <- obj$width + 1 + obj$padding
	obj$tracks <- append( obj$tracks, list(list( type='pileup', df=plot_df, y_bot=base_y-.1, y_top=base_y-.1+1 )) )
	obj$track_guides <- c( obj$track_guides, base_y )
	obj
}

tracks_annotation <- function( obj=tracks_create(), df ){
	gene_height        <- 0.0125
	exon_height        <- 0.05
	tss_marker_height  <- 0.20
	tss_marker_width   <- 0.05
	label_height       <- 0.1
	
	box_height <- gene_height+2*exon_height + tss_marker_height + label_height + 0.1
	
	annotation_sep <- .01 *(max(df$end) - min(df$start))
	
	plot_df <- df %>% 
		subset(type == 'gene') %>% 
		mutate( tss1_xmin = ifelse( strand=='-', end -   annotation_sep, start                   )
		       , tss1_xmax = ifelse( strand=='-', end                   , start +   annotation_sep)
		       , tss2_xmin = ifelse( strand=='-', end - 3*annotation_sep, start                   )
		       , tss2_xmax = ifelse( strand=='-', end                   , start + 3*annotation_sep)
		       , width     = end - start
		       , y_index   = 0 ) %>% 
		arrange(desc(width))
	
	for( i in 1:nrow(plot_df) ){
		element_to_place     <- plot_df[i, ]
		placed_elements      <- slice_head( plot_df, n=i-1 )
		overlapping_elements <- subset( placed_elements, (start <= element_to_place$end) & (end >= element_to_place$start) )
		plot_df$y_index[i] <- min(setdiff(0:1+max(overlapping_elements$y_index,0), unique(overlapping_elements$y_index)))
	}
	
	obj$width <- obj$width + (1+max(plot_df$y_index))*box_height + obj$padding
	y_base <- -(obj$width - obj$padding)
	
	plot_df <- plot_df %>% 
		mutate( y_bot = y_base + y_index*box_height + label_height ) %>% 
		mutate( y_top = y_bot + gene_height ) %>% 
		mutate( tss1_ymin = y_top 
		      , tss1_ymax = y_top + tss_marker_height 
		      , tss2_ymin = y_top + tss_marker_height - tss_marker_width
		      , tss2_ymax = y_top + tss_marker_height
		      , label_y   = y_bot - exon_height - label_height/2 )
	
	
	df_exons <- df %>% 
		subset(type == 'exon') %>% 
		select(-type) %>% 
		left_join(select(plot_df, gene_id, y_bot, y_top), by='gene_id')
	
	df_upper_exons <- df_exons %>% 
		mutate( xmin=start, xmax=end, ymin=y_top, ymax=y_top+exon_height, alpha=1.0/num_siblings, colour=gene_biotype ) %>% 
		select( xmin, xmax, ymin, ymax, alpha, colour )
	
	df_lower_exons <- df_exons %>%
		subset( cannonical == TRUE ) %>% 
		mutate( xmin=start, xmax=end, ymin=y_bot, ymax=y_bot-exon_height, alpha=1.0            , colour=gene_biotype ) %>% 
		select( xmin, xmax, ymin, ymax, alpha, colour )
	  
	obj$tracks <- append( obj$tracks, list(list( type='annotation', df=plot_df, df_exons=rbind(df_upper_exons, df_lower_exons) )) )
	obj
}

tracks_plot <- function( obj ){
	# obj <- track
	tracks <- obj$tracks
	
	bars <- keep(tracks, map_chr(tracks, ~.$type) == 'shade' ) %>% 
		map_dfr(function(d){
			df <- d$df
			if( d$direction == 'down' ){
				df$shade_top <- df$tab_bot
				if( d$range > 0 ){
					df$shade_bot <- tracks[[d$start_index + d$range]]$y_bot
				}else{
					df$shade_bot <- -obj$width
				}
			}else if( d$direction == 'up' ){
				df$shade_bot <- df$tab_top
				if( d$range > 0 ){
					df$shade_top <- tracks[[d$start_index - d$range]]$y_top
				}else{
					df$shade_top <- 0
				}
			}
			df
		})
	
	pileups <- keep(tracks, map_chr(tracks, ~.$type) == 'pileup' ) %>% 
		map(~.$df) %>% 
		bind_rows()
		
	annotations <- keep(tracks, map_chr(tracks, ~.$type) == 'annotation' ) %>% 
		map(~.$df) %>% 
		bind_rows()
	
	exon_annotations <- keep(tracks, map_chr(tracks, ~.$type) == 'annotation' ) %>% 
		map(~.$df_exons) %>% 
		bind_rows()
	
	xmin <- min(pileups$xmin)
	xmax <- max(pileups$xmax)
	guides <- data.frame(x=xmin, xend=xmax, y=obj$track_guides, yend=obj$track_guides)
	
	p <- ggplot()
	p <- p +
		geom_rect(aes(xmin=xmin, xmax=xmax, ymin=shade_bot, ymax=shade_top, fill=group, alpha=0.1*alpha), show.legend=FALSE, data=bars) +
		geom_rect(aes(xmin=xmin, xmax=xmax, ymin=tab_bot  , ymax=tab_top  , fill=group, alpha=1.0*alpha), show.legend=TRUE , data=bars)  
	if( nrow(annotations) != 0 ){
		p <- p +
			geom_rect(aes(xmin=start    , xmax=end      , ymin=y_bot    , ymax=y_top    , fill=gene_biotype       ), data=annotations) +
			geom_rect(aes(xmin=tss1_xmin, xmax=tss1_xmax, ymin=tss1_ymin, ymax=tss1_ymax, fill=gene_biotype       ), data=annotations) +
			geom_rect(aes(xmin=tss2_xmin, xmax=tss2_xmax, ymin=tss2_ymin, ymax=tss2_ymax, fill=gene_biotype       ), data=annotations) +
			geom_text(aes(x=(start+end)/2, y=label_y, label=label), data=annotations)
	}
	if( nrow(exon_annotations) != 0 ){
		p <- p +
			geom_rect(aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=colour, alpha=alpha), data=exon_annotations) 
	}
	p <- p + 
		geom_segment(aes(x  =xmin , xend=xmax, y   =y   , yend=y                                ), data=guides) +
		geom_rect   (aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=condition, alpha=alpha), data=pileups)
	p +
		scale_colour_tableau() +
		scale_fill_tableau() +
		theme_void()
}

print('Working out ranges to plot')
search_width <- 10000
search_peaks <- makeGRangesFromDataFrame(search_peaks, keep.extra.columns=TRUE)
data_ranges <- GenomicRanges::reduce( search_peaks + search_width, ignore.strand=TRUE )

print('Plotting ranges')
#future_map( seq_along(data_ranges), function(i){
for( i in seq_along(data_ranges) ){
	print(paste0(i, '/', length(data_ranges)))
	range_plot      <- data_ranges[i]
	range_chrom     <- seqnames(range_plot) %>% as.character()
	range_gtf       <- filter( gtf, (gtf$seqnames == range_chrom) & (gtf$start <= end(range_plot)) & (gtf$end >= start(range_plot)))
	range_gtf_genes <- filter( range_gtf, type=='gene' )
	
	# Expand search area for all gtf annotatins
	range_plot <- GRanges(seqnames=range_chrom, ranges=IRanges(start=min(range_gtf_genes$start, start(range_plot)), end=max(range_gtf_genes$end, end(range_plot))))
	
	scan_params     <- ScanBamParam( which=range_plot, what=scanBamWhat() )
	pileup_params   <- PileupParam( distinguish_nucleotides=FALSE, distinguish_strands=FALSE ) #PileupParam(distinguish_nucleotides=FALSE)
	
	range_peaks <- GenomicRanges::pintersect(search_peaks, range_plot, drop.nohit.ranges=TRUE)
	min_rank <- min(range_peaks$method_rank)
	
	if( length(subset(range_peaks, significant)) && nrow(range_gtf) ){
		sample_range_pileup <- map_dfr(seq_along(bam_files), function(i){
			sample_id       <- sample_ids[i]
			condition       <- sample2condition[[sample_id]]
			condition_count <- condition2count[[condition]]
			input           <- sample2input[[sample_id]]
			offset          <- condition_offset[[condition]] + !input
			libsize_factor  <- sample2libsize_factor[[sample_id]]
			pileup( bam_files[i], index=bam_index_files[i], scanBamParam=scan_params, pileupParam=pileup_params ) %>%
				select(-which_label) %>% 
				mutate( base=offset, sample=sample_id, condition=condition, input=input, libsize_factor=libsize_factor, alpha=1.0/condition_count ) %>% 
				mutate(count = count / libsize_factor)
		})
		if( length( sample_range_pileup )){
			sample_range_pileup <- sample_range_pileup %>% 
				mutate(normcount=count/max(count))
			
			track <- tracks_create() 
			
			range_de_peaks <- GenomicRanges::pintersect(de_peaks_gr, range_plot, drop.nohit.ranges=TRUE) %>% 
				as.data.frame() 
			if(nrow(range_de_peaks) != 0) {
				track <- tracks_shade_bar( track, range_de_peaks, group_column='method', direction='down', range=-1 )
			}
			
			if(nrow(range_gtf) != 0){
				track <- tracks_annotation( track, range_gtf )
			}
			
			range_cond_peaks <- GenomicRanges::pintersect(makeGRangesFromDataFrame(cond_peaks, keep.extra.columns=TRUE), range_plot, drop.nohit.ranges=TRUE) %>% 
				as.data.frame()
			for( c in ordered_conditions ){
				pc <- subset(range_cond_peaks   , condition==c)
				pi <- subset(sample_range_pileup, condition==c)
			  
				track <- track %>% 
					tracks_pileup(    subset(pi, input == FALSE), condition_column='condition' ) %>% 
					tracks_pileup(    subset(pi, input == TRUE ), condition_column='condition' )
				if( nrow(pc) > 0 ){
					track <- tracks_shade_bar( track, pc , group_column='method', direction='up', range=2 )
				}
			}
			
			p <- tracks_plot(track)
			dir.create( output_dir, showWarnings=TRUE, recursive=TRUE )
			save_name = paste0(output_dir, '/pileup_', min_rank, '_chr', as.character(seqnames(range_plot)), ':', start(range_plot), '-', end(range_plot), '.pdf')
			save_name = str_replace( save_name, 'chrchr', 'chr' )
			ggsave(save_name, plot=p)
			print(save_name)
			
			#pmap(list(range_gtf_genes$start, range_gtf_genes$end, ifelse(range_gtf_genes$name=='', range_gtf_genes$parent_gene, range_gtf_genes$name)), function(start, end, name){
			#  q <- p +
			#    xlim(start, end)
			#  ggsave(str_replace(save_name, '.pdf', paste0('_', name, '.pdf')), plot=q)
			#})
		}
	}
}
#})
#}
