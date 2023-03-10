library(GenomicRanges)
library(Rsamtools)
library(ggplot2)
library(ggthemes)
library(tidyverse)

LERP <- function(y1, y2, x){y1 + (y2-y1)*x}

pepr_peaks_filename   <- '../04_peakcalling/pepr/MAVSvsd103-467__PePr_chip1_peaks.bed'
genrich_filenames     <- list.files('../04_peakcalling/genrich/', full.names=TRUE)   %>% .[str_ends(.,'_peaks.narrowPeak')]
genrich_conditions    <- list.files('../04_peakcalling/genrich/', full.names=FALSE)  %>% .[str_ends(.,'_peaks.narrowPeak')]%>% str_remove('_peaks.narrowPeak')
idr_filenames         <- list.files('../04_peakcalling/idr/'    , full.names=TRUE  ) %>% .[str_ends(.,'_true.tsv')]
idr_conditions        <- list.files('../04_peakcalling/idr/'    , full.names=FALSE ) %>% .[str_ends(.,'_true.tsv')] %>% map_chr( ~str_split(.x, '_', n=2)[[1]][1] )
idr_elements          <- list.files('../04_peakcalling/idr/'    , full.names=FALSE ) %>% .[str_ends(.,'_true.tsv')] %>% map_chr( ~str_split(.x, '_', n=2)[[1]][2] ) %>% str_remove('_true.tsv')
bam_files             <- list.files('../03_aligned/', full.names=TRUE ) %>% .[str_ends(.,'star_aligned.bam')]
bam_index_files       <- list.files('../03_aligned/', full.names=TRUE ) %>% .[str_ends(.,'star_aligned.bam.bai')]
sample_ids            <- list.files('../03_aligned/', full.names=FALSE) %>% .[str_ends(.,'star_aligned.bam')] %>% str_remove('.star_aligned.bam')
treatment_conditions  <- 'MAVS'
control_condition     <- 'd103-467'
metadata_filename     <- '../metadata.tsv'
gff_filename          <- '../../reference/hg38.gff'
output_dir            <- '../04_peakcalling/analysis/plots/'
threads               <- 32
pepr_cuttoff          <- 1e-15
genrich_cuttoff       <- 1e-5
idr_cuttoff           <- 0.01

#pepr_peaks_filename   <- snakemake@input[['pepr_peaks']]
#genrich_filenames     <- snakemake@input[['genrich_peaks']]
#genrich_conditions    <- snakemake@params[['genrich_conditions']]
#idr_filenames         <- snakemake@input[['idr_peaks']]
#idr_conditions        <- snakemake@params[['idr_conditions']]
#idr_elements          <- snakemake@params[['idr_elements']]
#bam_files             <- snakemake@input[['bam_files']]
#bam_index_files       <- snakemake@input[['bam_index_files']]
#sample_ids            <- snakemake@params[['sample_ids']]
#metadata_filename     <- snakemake@params[['metadata']]
#treatment_conditions  <- snakemake@params[['treatment_conditions']]
#control_condition     <- snakemake@params[['control_conditions']]
#gff_filename          <- snakemake@input[['gff']]
#output_dir            <- snakemake@output[['plot_dir']]
#threads               <- snakemake@threads
#pepr_cuttoff          <- snakemake@params[['pepr_cuttoff']]
#genrich_cuttoff       <- snakemake@params[['genrich_cuttoff']]
#idr_cuttoff           <- snakemake@params[['idr_cuttoff']]

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
pepr_peaks <- map2( pepr_peaks_filename, treatment_conditions, function(pepr_file, treatment_condition){
  x <- read_tsv( pepr_file, col_names=broadpeak_colnames, show_col_types=FALSE ) %>%
    mutate( strand=ifelse(strand == '.', '*', strand) )
  GRanges( seqnames=x$chrom
         , ranges=IRanges(x$start, end=x$end, names=paste( 'pepr', condition, x$name, sep='_' ))
         , strand=x$strand
         , p_value=x$p_value
         , sig=x$p_value < pepr_cuttoff
         , method='pepr'
         , condition=treatment_condition )  
}) %>% do.call('c',.)

de_peaks        <- c(pepr_peaks)
de_search_peaks <- de_peaks

# Read single condition peaks
genrich_colnames <- c( 'chrom', 'start', 'end', 'name', 'scaled_auc', 'strand', 'auc', '-log10(pvalue)', 'qvalue', 'peak' )
genrich_peaks <- map2( genrich_filenames, genrich_conditions, function(genrich_file, condition){
  x <- read_tsv( genrich_file, col_names=genrich_colnames, show_col_types=FALSE ) %>% 
    mutate( strand=ifelse(strand == '.', '*', strand) )
  GRanges( seqnames=x$chrom
         , ranges=IRanges(x$start, end=x$end, names=paste( 'genrich', condition, x$name, sep='_' ))
         , strand=x$strand
         , p_value=10**(-x$`-log10(pvalue)`)
         , sig=10**(-x$`-log10(pvalue)`) < genrich_cuttoff
         , method='genrich'
         , condition=condition )
}) %>% do.call('c',.)
genrich_treatment <- subset( genrich_peaks, condition %in% treatment_conditions )
genrich_control   <- subset( genrich_peaks, condition %in% control_condition    )
nonoverlapping_ranges <- findOverlaps( genrich_treatment, genrich_control, select='first' ) %>% is.na()
genrich_unique_treatment_ranges <- subset( genrich_treatment, nonoverlapping_ranges ) %>% 
  subset(sig)

idr_colnames <- c( 'chrom', 'start', 'end', 'name', 'scaled_idr', 'strand', 'enrichment', 'p_value', 'q_value', 'peak', 'local_idr', 'global_idr', 's1_start', 's1_end', 's1_enrichment', 's1_summit', 's2_start', 's2_end', 's2_enrichment', 's2_summit')
idr_peaks <- pmap( list(idr_filenames, idr_conditions, idr_elements), function(idr_filename, idr_condition, idr_element){
  x <- read_tsv( idr_filename, col_names=idr_colnames, show_col_types = FALSE) %>%
    mutate( strand=ifelse(strand == '.', '*', strand) ) %>% 
    mutate( idr=2**(scaled_idr/-125) )
  GRanges( seqnames=x$chrom
         , ranges=IRanges(x$start, end=x$end, names=paste( 'idr', idr_condition, idr_element, seq_along(x$name), sep='_' ))
         , strand=x$strand
         , p_value=x$idr
         , sig=x$idr < idr_cuttoff
         , method='idr'
         , condition=idr_condition )
}) %>% do.call('c',.)
idr_treatment <- subset( idr_peaks, condition %in% treatment_conditions )
idr_control   <- subset( idr_peaks, condition %in% control_condition    )
nonoverlapping_ranges <- findOverlaps( idr_treatment, idr_control, select='first' ) %>% is.na()
idr_unique_treatment_ranges <- subset( idr_treatment, nonoverlapping_ranges ) %>% 
  subset(sig)

cond_search_peaks <- c( genrich_unique_treatment_ranges, idr_unique_treatment_ranges )
cond_peaks <- c( genrich_peaks, idr_peaks )

search_peaks <- c( de_search_peaks, cond_search_peaks )
peaks        <- c( de_peaks       , cond_peaks        )

# Work Out Ranges To plot
for( treatment_condition in treatment_conditions ){
  treatment_de_peaks     <- subset( de_peaks    , condition %in% c( treatment_condition, control_condition ) )
  treatment_cond_peaks   <- subset( cond_peaks  , condition %in% c( treatment_condition, control_condition ) )
  treatment_peaks        <- subset( peaks       , condition %in% c( treatment_condition, control_condition ) )
  treatment_search_peaks <- subset( search_peaks, condition %in% c( treatment_condition, control_condition ) )
  
  search_width <- 10000
  data_ranges <- GenomicRanges::reduce(treatment_search_peaks + search_width)
  
  for( i in seq_along(data_ranges) ){
    range_plot     <- data_ranges[i]
    range_chrom    <- seqnames(range_plot) %>% as.character()
    
    scan_params   <- ScanBamParam( which=range_plot, what=scanBamWhat() )
    pileup_params <- PileupParam( distinguish_nucleotides=FALSE, distinguish_strands=FALSE ) #PileupParam(distinguish_nucleotides=FALSE)
    
    range_peaks <- GenomicRanges::pintersect(treatment_search_peaks, range_plot, drop.nohit.ranges=TRUE)
    
    gff_y_start=min(condition_offset)-1.5
    range_gff <- gff %>% 
      filter( (gff$chrom == range_chrom) & (gff$start >= start(range_plot)) & (gff$end <= end(range_plot))) %>% 
      mutate( tss_x =ifelse(strand=='-', LERP(start,end,.995), LERP(start,end,.005))
            , tss_x2=ifelse(strand=='-', LERP(start,end,.975 ), LERP(start,end,.025))) %>% 
      mutate( width=end-start, y_offset=0 ) %>% 
      arrange(desc(width))
    for( i in range(nrow(range_gff)) ){
      element_to_place     <- slice( range_gff, i )
      placed_elements      <- slice_head( range_gff, n=i-1 )
      overlapping_elements <- subset( placed_elements, (start <= element_to_place$end) & (end >= element_to_place$start) )
      range_gff$y_offset[i] <- min(setdiff(0:1+max(overlapping_elements$y_offset,0), unique(overlapping_elements$y_offset)))
    }
    range_gff$y = gff_y_start - range_gff$y_offset
    
    print(range_gff)
    
    if( length(subset(range_peaks, sig))){
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
      
      range_de_peaks <- GenomicRanges::pintersect(treatment_de_peaks, range_plot, drop.nohit.ranges=TRUE) %>% 
        as.data.frame() 
      if(nrow(range_de_peaks) != 0) {
        range_de_peaks <- range_de_peaks %>% 
          mutate( method_offset=method_offset(method) ) %>% 
          mutate( shade_top = max(sample_range_pileup$base) + 1.5
                , shade_bot = min(sample_range_pileup$base) - 0.5
                , bar_bot   = max(sample_range_pileup$base) + 1.5 + method_offset
                , bar_top   = max(sample_range_pileup$base) + 1.5 + method_offset + .1 )
      }
  
      range_cond_peaks <- GenomicRanges::pintersect(treatment_cond_peaks, range_plot, drop.nohit.ranges=TRUE) %>% 
        as.data.frame()
      if( nrow(range_cond_peaks) ){ 
        range_cond_peaks <- range_cond_peaks %>% 
          left_join(enframe(condition_offset, name='condition', value='offset'), by='condition') %>% 
          mutate( method_offset=method_offset(method) ) %>% 
          mutate( shade_top = map_dbl(condition, ~condition_offset[[.x]]+ 2  )
                , shade_bot = map_dbl(condition, ~condition_offset[[.x]]-.5  )-method_offset
                , bar_bot   = map_dbl(condition, ~condition_offset[[.x]]-.5  )-method_offset -.1
                , bar_top   = map_dbl(condition, ~condition_offset[[.x]]-.5  )-method_offset
                , alpha     = ifelse(method == 'idr', 1.0/max(table(idr_condition)), 1.0) )
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
           geom_rect(aes(xmin=start, xmax=end, ymin=shade_bot, ymax=shade_top, fill=method, alpha=alpha*0.1), show.legend=FALSE, data=range_cond_peaks) +
           geom_rect(aes(xmin=start, xmax=end, ymin=bar_bot  , ymax=bar_top  , fill=method, alpha=alpha    ), show.legend=TRUE , data=range_cond_peaks)
      }
      p <- p +
        geom_segment(aes(x=xmin, xend=xmax, y=y, yend=y), data=line_df) +
        geom_rect(aes(xmin=pos-0.5, xmax=pos+0.5, ymin=base, ymax=base+normcount, fill=condition, alpha=alpha), data=sample_range_pileup) +
        geom_segment(aes(x=start, xend=end   , y=y   , yend=y   , colour=level), lineend='round', linewidth=2, data=range_gff) +
        geom_segment(aes(x=tss_x, xend=tss_x , y=y   , yend=y+.3, colour=level), lineend='round', linewidth=1, data=range_gff) +
        geom_segment(aes(x=tss_x, xend=tss_x2, y=y+.3, yend=y+.3, colour=level), lineend='round', linewidth=1, data=range_gff) +
        geom_text(aes(x=(start+end)/2, y=y-.15, label=label, colour=level), data=range_gff) +
        theme_void()
      
      dir.create( output_dir, showWarnings=TRUE, recursive=TRUE )
      save_name = paste0(output_dir, '/pileup_', treatment_condition, '_chr', as.character(seqnames(range_plot)), ':', start(range_plot), '-', end(range_plot), '.pdf')
      ggsave(save_name, plot=p)
      
      pmap(list(range_gff$start, range_gff$end, ifelse(range_gff$name=='', range_gff$metadata, range_gff$name)), function(start, end, name){
        q <- p +
          xlim(start, end)
        ggsave(str_replace(save_name, '.pdf', paste0('_', name, '.pdf')), plot=q)
      })
    }
  }
}
