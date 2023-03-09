library(GenomicRanges)
library(Rsamtools)
library(ggplot2)
library(ggthemes)
library(tidyverse)

LERP <- function(y1, y2, x){y1 + (y2-y1)*x}

pepr_peaks_filename <- '../04_peakcalling/pepr/MAVSvsd103-467__PePr_chip1_peaks.bed'
genrich_filenames   <- list.files('../04_peakcalling/genrich/', full.names=TRUE) %>% .[str_ends(.,'_peaks.narrowPeak')]
genrich_conditions  <- list.files('../04_peakcalling/genrich/', full.names=FALSE) %>% .[str_ends(.,'_peaks.narrowPeak')]%>% str_remove('_peaks.narrowPeak')
bam_files           <- list.files('../03_aligned/', full.names=TRUE ) %>% .[str_ends(.,'star_aligned.bam')]
bam_index_files     <- list.files('../03_aligned/', full.names=TRUE ) %>% .[str_ends(.,'star_aligned.bam.bai')]
sample_ids          <- list.files('../03_aligned/', full.names=FALSE) %>% .[str_ends(.,'star_aligned.bam')] %>% str_remove('.star_aligned.bam')
metadata_filename   <- '../metadata.tsv'
gff_filename        <- '../../reference/hg38.gff'
output_dir          <- '../04_peakcalling/analysis/plots/'
threads             <- 32
pepr_cuttoff        <- 1e-15

#pepr_peaks_filename <- snakemake@input[['pepr_peaks']]
#genrich_filenames   <- snakemake@input[['genrich_peaks']]
#genrich_conditions  <- snakemake@params[['genrich_conditions']]
#bam_files           <- snakemake@input[['bam_files']]
#bam_index_files     <- snakemake@input[['bam_index_files']]
#sample_ids          <- snakemake@params[['sample_ids']]
#metadata_filename   <-snakemake@params[['metadata']]
#gff_filename        <- snakemake@input[['gff']]
#output_dir          <- snakemake@output[['plot_dir']]
#threads             <- snakemake@threads
#pepr_cuttoff        <- snakemake@params[['pepr_cuttoff']]

# Read metadata
metadata <- read_tsv( metadata_filename, show_col_types=FALSE )
sample2condition <- select(metadata, sample_id, condition) %>% deframe()
sample2input     <- select(metadata, sample_id, method   ) %>% mutate(method=method=='Input') %>% deframe()
condition2count  <- table(metadata$condition) / 2 #NOTE(KIM): Assuming equal input and ip

conds      = unique(sample2condition)
num_conditions = length(conds)
condition_offset = set_names((seq_along(conds) - 1)*-3, conds)

# Read in Annotation
name_matcher <- function(m){
  res <- ''
  if( str_detect(m, 'Name=([^;]+);')){
    res <- str_match(m, 'Name=([^;]+);')[2]
  }
  res
}
gff_columns <- c( 'chrom', 'label', 'level', 'start', 'end', 'score', 'strand', 'something', 'metadata' )
gff <- read_tsv( gff_filename, col_names=gff_columns, comment='#', show_col_types=FALSE )
gff <- gff %>%
  filter( str_starts( metadata, 'ID=gene:' ) ) %>%
  mutate(name=map_chr(metadata, name_matcher)) %>% 
  mutate(metadata=map_chr(metadata, ~str_match( .x, 'ID=gene:([^;]+);' )[2] )) %>% 
  mutate(label=ifelse(name!='', paste0(metadata, ' - ', name), metadata))

# Read in DiffBind results
broadpeak_colnames <- c( 'chrom', 'start', 'end', 'name', 'score', 'strand', 'enrichment', 'p_value', 'q_value' )
pepr_peaks <- read_tsv( pepr_peaks_filename, col_names=broadpeak_colnames, show_col_types=FALSE ) %>%
	mutate( strand=ifelse(strand == '.', '*', strand) )
pepr_peaks <- GRanges( seqnames=pepr_peaks$chrom
                     , ranges=IRanges(pepr_peaks$start, end=pepr_peaks$end, names=pepr_peaks$name)
                     , strand=pepr_peaks$strand
                     , p_value=pepr_peaks$p_value
                     , plot=pepr_peaks$p_value < pepr_cuttoff
                     , method='pepr', condition='all' )
de_peaks <- c(pepr_peaks)

# Read single condition peaks
genrich_colnames <- c( 'chrom', 'start', 'end', 'name', 'scaled_auc', 'strand', 'auc', '-log10(pvalue)', 'qvalue', 'peak' )
genrich_peaks <- map2( genrich_filenames, genrich_conditions, function(filename, condition){
  x <- read_tsv( filename, col_names=genrich_colnames, show_col_types=FALSE ) %>% 
    mutate( strand=ifelse(strand == '.', '*', strand) )
  GRanges( seqnames=x$chrom
         , ranges=IRanges(x$start, end=x$end, names=x$name)
         , strand=x$strand
         , p_value=10**(-x$`-log10(pvalue)`)
         , method='genrich'
         , condition=condition)  
})


cond_peaks <- c(do.call('c',genrich_peaks))

# Work Out Ranges To plot
peaks <- c(de_peaks,cond_peaks)
search_width <- 10000
data_ranges <- GenomicRanges::reduce(peaks + search_width)

for( i in seq_along(data_ranges) ){
  range_plot     <- data_ranges[i]
  range_chrom    <- seqnames(range_plot) %>% as.character()
  
  scan_params   <- ScanBamParam( which=range_plot, what=scanBamWhat() )
  pileup_params <- PileupParam( distinguish_nucleotides=FALSE, distinguish_strands=FALSE ) #PileupParam(distinguish_nucleotides=FALSE)
  
  range_peaks <- GenomicRanges::pintersect(de_peaks, range_plot, drop.nohit.ranges=TRUE)
  
  if( length(subset(range_peaks, plot))){
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
    
    range_de_peaks <- GenomicRanges::pintersect(de_peaks, range_plot, drop.nohit.ranges=TRUE) %>% 
      as.data.frame() 
    if(nrow(range_de_peaks) != 0) {
      range_de_peaks <- range_de_peaks %>% 
        mutate(method_offset=method_offset(method)) %>% 
        mutate( shade_top = max(sample_range_pileup$base) + 1.5
              , shade_bot = min(sample_range_pileup$base) - 0.5
              , bar_bot   = max(sample_range_pileup$base) + 1.5 + method_offset
              , bar_top   = max(sample_range_pileup$base) + 1.5 + method_offset + .1 )
    }

    range_cond_peaks <- GenomicRanges::pintersect(cond_peaks, range_plot, drop.nohit.ranges=TRUE) %>% 
      as.data.frame()
    if( nrow(range_cond_peaks) ){ 
      range_cond_peaks <- range_cond_peaks %>% 
        left_join(enframe(condition_offset, name='condition', value='offset'), by='condition') %>% 
        group_by(method) %>% 
        mutate(display_offest=method_offset(method))
      range_cond_peaks$shade_top <- map_dbl(range_cond_peaks$condition, ~condition_offset[[.x]]+ 2  )
      range_cond_peaks$shade_bot <- map_dbl(range_cond_peaks$condition, ~condition_offset[[.x]]-.5  )
      range_cond_peaks$bar_bot   <- map_dbl(range_cond_peaks$condition, ~condition_offset[[.x]]-.6  )-range_cond_peaks$display_offest
      range_cond_peaks$bar_top   <- map_dbl(range_cond_peaks$condition, ~condition_offset[[.x]]-.5  )-range_cond_peaks$display_offest
    }
    
    sample_range_pileup <- map_dfr(seq_along(bam_files), function(i){
      sample_id       <- sample_ids[i]
      condition       <- sample2condition[[sample_id]]
      condition_count <- condition2count[[condition]]
      input           <- sample2input[[sample_id]]
      offset          <- condition_offset[[condition]] + !input
      pileup( bam_files[i], index=bam_index_files[i], scanBamParam=scan_params, pileupParam=pileup_params ) %>%
        select(-which_label) %>% 
        mutate( base=offset, sample=sample_id, condition=condition, input=input, alpha=1.0/condition_count )
    }) %>% 
      mutate(normcount=count/max(count))
    #%>% 
    #  mutate(count=ifelse( strand=='-', -count, count))
    
    
    range_gff <- gff %>% 
      filter( (gff$chrom == range_chrom) & (gff$start >= start(range_plot)) & (gff$end <= end(range_plot))) %>% 
      mutate( y=min(condition_offset)-1.5) %>% 
      mutate( tss_x =ifelse(strand=='-', LERP(start,end,.995), LERP(start,end,.005))
            , tss_x2=ifelse(strand=='-', LERP(start,end,.975 ), LERP(start,end,.025)))
    # Plot
    line_df = map_dfr(unique(sample_range_pileup$base), function(b){
      data.frame(xmin=start(range_plot), xmax=end(range_plot), y=b)
    })
    
    p <- ggplot()
    if( nrow(range_de_peaks) != 0 ){
      p <- p +
        geom_rect(aes(xmin=start, xmax=end, ymin=shade_bot, ymax=shade_top, fill=method), alpha=0.1, show.legend=FALSE, data=range_de_peaks) +
        geom_rect(aes(xmin=start, xmax=end, ymin=bar_bot  , ymax=bar_top  , fill=method), alpha=1.0, show.legend=TRUE , data=range_de_peaks)
    }
    if( nrow(range_cond_peaks) != 0 ){
       p <- p +
         geom_rect(aes(xmin=start, xmax=end, ymin=shade_bot, ymax=shade_top, fill=method), alpha=0.1, show.legend=FALSE, data=range_cond_peaks) +
         geom_rect(aes(xmin=start, xmax=end, ymin=bar_bot  , ymax=bar_top  , fill=method), alpha=1.0, show.legend=TRUE , data=range_cond_peaks)
    }
    p <- p +
      geom_segment(aes(x=xmin, xend=xmax, y=y, yend=y), data=line_df) +
      geom_rect(aes(xmin=pos-0.5, xmax=pos+0.5, ymin=base, ymax=base+normcount, fill=condition, alpha=alpha), data=sample_range_pileup) +
      geom_segment(aes(x=start, xend=end   , y=y   , yend=y   ), lineend='round', linewidth=2, data=range_gff) +
      geom_segment(aes(x=tss_x, xend=tss_x , y=y   , yend=y+.3), lineend='round', linewidth=1, data=range_gff) +
      geom_segment(aes(x=tss_x, xend=tss_x2, y=y+.3, yend=y+.3), lineend='round', linewidth=1, data=range_gff) +
      geom_text(aes(x=(start+end)/2, y=y-.15, label=label), data=range_gff) +
      theme_void()
    
    dir.create( output_dir, showWarnings=TRUE, recursive=TRUE )
    save_name = paste0(output_dir, '/pileup_chr', as.character(seqnames(range_plot)), ':', start(range_plot), '-', end(range_plot), '.pdf')
    ggsave(save_name, plot=p)
  }
}
