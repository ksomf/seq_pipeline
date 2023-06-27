library(GenomicRanges)
library(Rsamtools)
library(ggplot2)
library(furrr)
library(rtracklayer)
library(tidyverse)

source('workflow/scripts/tracks.R')

theme_set(theme_void())

LERP <- function(y1, y2, x){y1 + (y2-y1)*x}

#condition_peaks <- '04_peakcalling/analysis/condition_peaks.tsv'
#diffbind_peaks  <- '04_peakcalling/analysis/diffbind_peaks.tsv'
#manual_peaks    <- '04_peakcalling/analysis/manual_peaks.tsv'
#bam_files       <- list.files('03_aligned/', full.names=TRUE ) %>% .[str_ends(.,'star_aligned.bam')]
#bam_index_files <- list.files('03_aligned/', full.names=TRUE ) %>% .[str_ends(.,'star_aligned.bam.bai')]
#library_sizes   <- '04_peakcalling/analysis/library_sizes.tsv'
#gtf_filename    <- '../reference/ucsc/hg38.gtf'
#
#sample_ids            <- list.files('03_aligned/', full.names=FALSE) %>% .[str_ends(.,'star_aligned.bam')] %>% str_remove('.star_aligned.bam')
#treatment_conditions  <- '{condition_1}'
#control_condition     <- '{condition_2}'
#metadata_filename     <- 'metadata.tsv'
#gtf_database          <- 'ucsc'
#
#output_dir            <- '04_peakcalling/analysis/plots/'
#output_signal         <- '04_peakcalling/analysis/summary.txt'
#
#threads               <- 32

condition_peaks <- snakemake@input[['condition_peaks']]
diffbind_peaks  <- snakemake@input[['diffbind_peaks' ]]
manual_peaks    <- snakemake@input[['manual_peaks'   ]]
bam_files       <- snakemake@input[['bam_files'      ]]
bam_index_files <- snakemake@input[['bam_index_files']]
library_sizes   <- snakemake@input[['library_sizes'  ]]
gtf_filename    <- snakemake@input[['gtf'            ]]

treatment_conditions  <- snakemake@params[['treatment_conditions']]
control_condition     <- snakemake@params[['control_condition'   ]]
metadata_filename     <- snakemake@params[['metadata'            ]]
sample_ids            <- snakemake@params[['sample_ids'          ]]
gtf_database          <- snakemake@params[['database'            ]]

output_dir              <- snakemake@output[['plot_dir']]
output_signal           <- snakemake@output[['summary_file']]

threads                 <- snakemake@threads

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

print('Reading Manual Peaks')
m_peaks         <- read_tsv(manual_peaks, show_col_types=FALSE) %>% mutate(siblings=1) %>% rename(any_of(c(chrom='seqnames')))
m_peaks_gr      <- makeGRangesFromDataFrame( m_peaks, keep.extra.columns=TRUE )

print('Reading vs Input Peakcaller Results')
cond_peaks <- read_tsv(condition_peaks, show_col_types=FALSE) %>% rename(any_of(c(chrom='seqnames')))
print(colnames(de_peaks))
print(colnames(cond_peaks))
print(colnames(m_peaks))
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

# search_peaks <- rbind( de_search_peaks, select(cond_search_peaks, -width) ) %>% 
# 	group_by(method) %>% 
# 	mutate(method_rank=seq_len(n()))
search_peak_deq <- de_search_peaks %>% 
	subset( method == 'deq' ) %>% 
	arrange( desc(significant), stat ) %>% 
	slice_head( n=30 )
search_peaks_rest <- de_search_peaks %>%
	subset( method != 'deq' ) %>%
	group_by(method) %>%
	slice_head( n=10 )
search_peaks <- rbind( select(search_peak_deq  , any_of(colnames(m_peaks)))
                     , select(search_peaks_rest, any_of(colnames(m_peaks)))
                     , m_peaks )

peaks <- rbind( select(de_peaks  , any_of(colnames(m_peaks)))
              , select(cond_peaks, any_of(colnames(m_peaks)))
              , m_peaks )

print('Working out ranges to plot')
search_width <- 10000
search_peaks <- makeGRangesFromDataFrame( search_peaks, keep.extra.columns=TRUE )
peaks        <- makeGRangesFromDataFrame( peaks       , keep.extra.columns=TRUE )
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
	
	range_peaks     <- GenomicRanges::pintersect(search_peaks, range_plot, drop.nohit.ranges=TRUE)
	displayed_peaks <- GenomicRanges::pintersect(peaks, range_plot, drop.nohit.ranges=TRUE)
	
	if( length(range_peaks) && nrow(range_gtf) ){
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
				mutate(count = count / libsize_factor) %>% 
				group_by( condition, input, pos ) %>% 
				summarise( min=min(count), max=max(count), median=median(count), q1=quantile(count, .25), q3=quantile(count, .75), mean=mean(count), std=sd(count), .groups='drop' )
		})
		if( length( sample_range_pileup )){
			norm <- sample_range_pileup %>% 
				group_by(input) %>% 
				slice_max(max, n=1, with_ties=F) %>% 
				ungroup() %>% 
				select(input, max) %>% 
				deframe()
			colour_scheme1 <- khroma::colour('muted')(9)
			colour_map <- c(idr=colour_scheme1[[2]], deq=colour_scheme1[[1]], input=colour_scheme1[[3]], ip=colour_scheme1[[4]])
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
				pi <- subset(sample_range_pileup, condition==c) %>% 
					mutate(input_var=ifelse(input, 'input', 'ip'))
				
				norm <- set_names(norm, c('ip', 'input'))
				
				pi_ip <- pi %>% 
					subset(input == FALSE) %>% 
					mutate(across(all_of(c('min', 'max', 'median', 'q1', 'q3', 'mean')), ~./norm[['ip']]))
			  
				pi_input <- pi %>% 
					subset(input == TRUE) %>% 
					mutate(across(all_of(c('min', 'max', 'median', 'q1', 'q3', 'mean')), ~./norm[['input']]))
				
				track <- track %>% 
					tracks_pileup_shade( pi_ip   , separation_variable='input_var', position_variable='pos', norm=norm['ip']   , colour='ip'   , axis_label=c ) %>% 
					tracks_pileup_shade( pi_input, separation_variable='input_var', position_variable='pos', norm=norm['input'], colour='input'  )
				if( nrow(pc) > 0 ){
					track <- tracks_shade_bar( track, pc , group_column='method', direction='up', range=2 )
				}
			}
			
			range_manual_peaks <- GenomicRanges::pintersect(m_peaks_gr, range_plot, drop.nohit.ranges=TRUE) %>% as.data.frame()
			if( nrow(range_manual_peaks) != 0 ){
				print('Adding manual hit')
				track <- tracks_shade_bar( track, range_manual_peaks, group_column='condition', direction='up', range=-1 )
			}
			
			p <- tracks_plot(track)
			p <- p +
				scale_fill_manual  (values=colour_map) +
				scale_colour_manual(values=colour_map)
			dir.create( output_dir, showWarnings=TRUE, recursive=TRUE )
			save_name = paste0(output_dir, '/pileup_chr', as.character(seqnames(range_plot)), '.', start(range_plot), '-', end(range_plot), '.svg')
			save_name = str_replace( save_name, 'chrchr', 'chr' )
			print(save_name)
			ggsave(save_name, plot=p, width=24, height=8)
			ggsave(str_replace(save_name, '.svg', '.png'), plot=p, width=24, height=8)
			as.data.frame(displayed_peaks) %>% write_tsv(str_replace( save_name, '.svg', '.tsv' ))
			
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
