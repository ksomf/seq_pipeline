library(GenomicRanges)
library(Rsamtools)
library(ggplot2)
library(ggthemes)
library(furrr)
library(tidyverse)

LERP <- function(y1, y2, x){y1 + (y2-y1)*x}

#pepr_peaks_filename   <- '../04_peakcalling/pepr/MAVSvsd103-467__PePr_chip1_peaks.bed'
#genrich_filenames     <- list.files('../04_peakcalling/genrich/', full.names=TRUE)   %>% .[str_ends(.,'_peaks.narrowPeak')]
#genrich_conditions    <- list.files('../04_peakcalling/genrich/', full.names=FALSE)  %>% .[str_ends(.,'_peaks.narrowPeak')]%>% str_remove('_peaks.narrowPeak')
#idr_filenames         <- list.files('../04_peakcalling/idr/'    , full.names=TRUE  ) %>% .[str_ends(.,'_true.tsv')]
#idr_conditions        <- list.files('../04_peakcalling/idr/'    , full.names=FALSE ) %>% .[str_ends(.,'_true.tsv')] %>% map_chr( ~str_split(.x, '_', n=2)[[1]][1] )
#idr_elements          <- list.files('../04_peakcalling/idr/'    , full.names=FALSE ) %>% .[str_ends(.,'_true.tsv')] %>% map_chr( ~str_split(.x, '_', n=2)[[1]][2] ) %>% str_remove('_true.tsv')
#bam_files             <- list.files('../03_aligned/', full.names=TRUE ) %>% .[str_ends(.,'star_aligned.bam')]
#bam_index_files       <- list.files('../03_aligned/', full.names=TRUE ) %>% .[str_ends(.,'star_aligned.bam.bai')]
#sample_ids            <- list.files('../03_aligned/', full.names=FALSE) %>% .[str_ends(.,'star_aligned.bam')] %>% str_remove('.star_aligned.bam')
#treatment_conditions  <- 'MAVS'
#control_condition     <- 'd103-467'
#metadata_filename     <- '../metadata.tsv'
#gff_filename          <- '../../reference/hg38.gff'
#output_dir            <- '../04_peakcalling/analysis/plots/'
#threads               <- 32
#pepr_cuttoff          <- 1e-15
#genrich_cuttoff       <- 1e-5
#idr_cuttoff           <- 0.01

pepr_peaks_filename     <- snakemake@input[['pepr_peaks']]
thor_peaks_filename     <- snakemake@input[['thor_peaks']]
diffbind_peaks_filename <- snakemake@input[['diffbind_peaks']]
genrich_filenames       <- snakemake@input[['genrich_peaks']]
genrich_conditions      <- snakemake@params[['genrich_conditions']]
idr_filenames           <- snakemake@input[['idr_peaks']]
idr_conditions          <- snakemake@params[['idr_conditions']]
idr_elements            <- snakemake@params[['idr_elements']]
bam_files               <- snakemake@input[['bam_files']]
bam_index_files         <- snakemake@input[['bam_index_files']]
sample_ids              <- snakemake@params[['sample_ids']]
metadata_filename       <- snakemake@params[['metadata']]
treatment_conditions    <- snakemake@params[['treatment_conditions']]
control_condition       <- snakemake@params[['control_conditions']]
gff_filename            <- snakemake@input[['gff']]
output_dir              <- snakemake@output[['plot_dir']]
threads                 <- snakemake@threads
pepr_cuttoff            <- snakemake@params[['pepr_cuttoff']]
thor_cuttoff            <- snakemake@params[['thor_cuttoff']]
genrich_cuttoff         <- snakemake@params[['genrich_cuttoff']]
idr_cuttoff             <- snakemake@params[['idr_cuttoff']]
diffbind_cuttoff        <- snakemake@params[['diffbind_cuttoff']]

ordered_conditions    <- c(treatment_conditions, control_condition)

plan(multisession, workers=threads)
#TODO normalise bam files against each other using RPKM

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
gff_full <- read_tsv( gff_filename, col_names=gff_columns, comment='#', show_col_types=FALSE )
gff_genes <- gff_full %>%
  filter( str_starts( metadata, 'ID=gene:' ) ) %>%
  mutate(name=map_chr(metadata, name_matcher)) %>% 
  mutate(parent_gene=map_chr(metadata, ~str_match( .x, 'ID=gene:([^;]+);' )[2] )) %>%
  mutate(label=ifelse(name!='', paste0(parent_gene, ' - ', name), parent_gene)) %>% 
  select(-metadata)

transcript2gene <- gff_full %>%
  filter( str_starts( metadata, 'ID=transcript:' ) ) %>%
  mutate(parent_gene=map_chr(metadata, ~str_match( .x, 'Parent=gene:([^;]+);' )[2] )) %>% 
  mutate(metadata=map_chr(metadata, ~str_match( .x, 'ID=transcript:([^;]+);' )[2] )) %>%
  select(metadata, parent_gene) %>%
  distinct() %>% 
  rename(parent_transcript='metadata') %>% 
  group_by(parent_gene) %>% 
  mutate( num_siblings=n()
        , cannonical=seq_len(n()) == 1 )

gff_exons <- gff_full %>%
  filter( str_starts( metadata, 'Parent=transcript:' ) ) %>%
  mutate(parent_transcript=map_chr(metadata, ~str_match( .x, 'Parent=transcript:([^;]+);' )[2] )) %>% 
  left_join(transcript2gene, by='parent_transcript') %>% 
  select(-parent_transcript, -metadata)

gff <- bind_rows( list(gene=gff_genes, exons=gff_exons), .id='meta_level' )

# Read in DiffBind results
broadpeak_colnames <- c( 'chrom', 'start', 'end', 'name', 'score', 'strand', 'enrichment', 'p_value', 'q_value' )
pepr_peaks <- map2( pepr_peaks_filename, treatment_conditions, function(pepr_file, treatment_condition){
  x <- read_tsv( pepr_file, col_names=broadpeak_colnames, show_col_types=FALSE ) %>%
    mutate( strand=ifelse(strand == '.', '*', strand) )
  GRanges( seqnames=x$chrom
         , ranges=IRanges(x$start, end=x$end, names=paste( 'pepr', condition, x$name, sep='_' ))
         , strand=x$strand
         , p_value=x$p_value
         , alpha=1.0
         , sig=x$p_value < pepr_cuttoff
         , method='pepr'
         , condition=treatment_condition )  
}) %>% do.call('c',.)

thor_colnames <- c( 'chrom', 'start', 'end', 'name', 'score', 'strand', 'enrichment', 'p_value', 'q_value' )
thor_peaks <- map2( thor_peaks_filename, treatment_conditions, function(thor_file, treatment_condition){
  x <- read_tsv( thor_file, col_names=broadpeak_colnames, show_col_types=FALSE ) %>%
    mutate( strand=ifelse(strand == '.', '*', strand) )
  GRanges( seqnames=x$chrom
         , ranges=IRanges(x$start, end=x$end, names=paste( 'thor', condition, x$name, sep='_' ))
         , strand=x$strand
         , p_value=x$p_value
         , alpha=1.0
         , sig=x$p_value < thor_cuttoff
         , method='thor'
         , condition=treatment_condition )  
}) %>% do.call('c',.)

x <- read_tsv( diffbind_peaks_filename, show_col_types=FALSE ) %>%
    mutate( strand=ifelse(strand == '.', '*', strand) )
diffbind_peaks <- GRanges( seqnames=x$chrom
           , ranges=IRanges(x$start, end=x$end, names=paste( 'diffbind', condition, x$name, sep='_' ))
           , strand=x$strand
           , p_value=x$p_value
           , alpha=1.0
           , sig=x$p_value < diffbind_cuttoff
           , method='macs2-diffbind'
           , condition=treatment_condition )  

de_peaks        <- c(pepr_peaks, thor_peaks, diffbind_peaks)
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
         , alpha=1.0
         , method='genrich'
         , condition=condition )
}) %>% do.call('c',.)
genrich_treatment <- subset( genrich_peaks, condition %in% treatment_conditions )
genrich_control   <- subset( genrich_peaks, condition %in% control_condition    )
nonoverlapping_ranges <- findOverlaps( genrich_treatment, genrich_control, select='first' ) %>% is.na()
genrich_unique_treatment_ranges <- subset( genrich_treatment, nonoverlapping_ranges ) %>% 
  subset(sig)

idr_colnames <- c( 'chrom', 'start', 'end', 'name', 'scaled_idr', 'strand', 'enrichment', 'p_value', 'q_value', 'peak', 'local_idr', 'global_idr', 's1_start', 's1_end', 's1_enrichment', 's1_summit', 's2_start', 's2_end', 's2_enrichment', 's2_summit')
condition2idr_number <- table(idr_conditions)
idr_peaks <- pmap( list(idr_filenames, idr_conditions, idr_elements), function(idr_filename, idr_condition, idr_element){
  x <- read_tsv( idr_filename, col_names=idr_colnames, show_col_types = FALSE) %>%
    mutate( strand=ifelse(strand == '.', '*', strand) ) %>% 
    mutate( idr=2**(scaled_idr/-125) )
  GRanges( seqnames=x$chrom
         , ranges=IRanges(x$start, end=x$end, names=paste( 'idr', idr_condition, idr_element, seq_along(x$name), sep='_' ))
         , strand=x$strand
         , p_value=x$idr
         , sig=x$idr < idr_cuttoff
         , alpha=1.0/condition2idr_number[idr_condition]
         , method='macs2-idr'
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
    rename( group=group_column ) %>% 
    mutate( tab_top=map_dbl(group, ~base_y+tab_height*(group2num[[.x]]  ))
          , tab_bot=map_dbl(group, ~base_y+tab_height*(group2num[[.x]]-1)) ) %>% 
    mutate( xmin=start, xmax=end ) %>%
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
  gene_height        <- 0.05
  exon_height        <- 0.05
  tss_marker_height  <- 0.20
  tss_marker_width   <- 0.05
  label_height       <- 0.1
  
  box_height <- gene_height+2*exon_height + tss_marker_height + label_height + 0.1
  
  annotation_sep <- .01 *(max(df$end) - min(df$start))
  
  plot_df <- df %>% 
   subset(meta_level == 'gene')%>% 
    mutate( tss1_xmin = ifelse( strand=='-', end -   annotation_sep, start                   )
          , tss1_xmax = ifelse( strand=='-', end                   , start +   annotation_sep)
          , tss2_xmin = ifelse( strand=='-', end - 3*annotation_sep, start                   )
          , tss2_xmax = ifelse( strand=='-', end                   , start + 3*annotation_sep)
          , width     = end - start
          , y_index   = 0 ) %>% 
    arrange(desc(width))
  
  for( i in 1:nrow(plot_df) ){
    element_to_place     <- slice( plot_df, i )
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
    subset(meta_level == 'exons' & level == 'exon') %>% 
    select(-level) %>% 
    left_join(select(plot_df, parent_gene, level, y_bot, y_top), by='parent_gene')
  
  df_upper_exons <- df_exons %>% 
    mutate( xmin=start, xmax=end, ymin=y_top, ymax=y_top+exon_height, alpha=1.0/num_siblings, colour=level ) %>% 
    select( xmin, xmax, ymin, ymax, alpha, colour )
  
  df_lower_exons <- df_exons %>%
    subset( cannonical == TRUE ) %>% 
    mutate( xmin=start, xmax=end, ymin=y_bot, ymax=y_bot-exon_height, alpha=1.0            , colour=level ) %>% 
    select( xmin, xmax, ymin, ymax, alpha, colour )
    
  obj$tracks <- append( obj$tracks, list(list( type='annotation', df=plot_df, df_exons=rbind(df_upper_exons, df_lower_exons) )) )
  obj
}

tracks_plot <- function( obj ){
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
  p <- p +
    geom_rect(aes(xmin=start    , xmax=end      , ymin=y_bot    , ymax=y_top    , fill=level             ), data=annotations) +
    geom_rect(aes(xmin=tss1_xmin, xmax=tss1_xmax, ymin=tss1_ymin, ymax=tss1_ymax, fill=level             ), data=annotations) +
    geom_rect(aes(xmin=tss2_xmin, xmax=tss2_xmax, ymin=tss2_ymin, ymax=tss2_ymax, fill=level             ), data=annotations) +
    geom_rect(aes(xmin=xmin     , xmax=xmax     , ymin=ymin     , ymax=ymax     , fill=colour, alpha=alpha), data=exon_annotations) +
    geom_text(aes(x=(start+end)/2, y=label_y, label=label), data=annotations)
  p <- p + 
    geom_segment(aes(x  =xmin , xend=xmax, y   =y   , yend=y                                ), data=guides) +
    geom_rect   (aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=condition, alpha=alpha), data=pileups)
  p +
    scale_colour_tableau() +
    scale_fill_tableau() +
    theme_void()
}

# Work Out Ranges To plot
#for( treatment_condition in treatment_conditions ){
#  treatment_de_peaks     <- subset( de_peaks    , condition %in% c( treatment_condition, control_condition ) )
#  treatment_cond_peaks   <- subset( cond_peaks  , condition %in% c( treatment_condition, control_condition ) )
#  treatment_peaks        <- subset( peaks       , condition %in% c( treatment_condition, control_condition ) )
#  treatment_search_peaks <- subset( search_peaks, condition %in% c( treatment_condition, control_condition ) )
  
search_width <- 10000
data_ranges <- GenomicRanges::reduce(search_peaks + search_width)

#future_map( seq_along(data_ranges), function(i){
for( i in seq_along(data_ranges) ){
  range_plot      <- data_ranges[i]
  range_chrom     <- seqnames(range_plot) %>% as.character()
  range_gff       <- filter( gff, (gff$chrom == range_chrom) & (gff$start <= end(range_plot)) & (gff$end >= start(range_plot)))
  range_gff_genes <- filter( range_gff, meta_level=='gene' )
  
  # Expand search area for all gff annotatins
  range_plot <- GRanges(seqnames=range_chrom, ranges=IRanges(start=min(range_gff_genes$start, start(range_plot)), end=max(range_gff_genes$end, end(range_plot))))
  
  scan_params     <- ScanBamParam( which=range_plot, what=scanBamWhat() )
  pileup_params   <- PileupParam( distinguish_nucleotides=FALSE, distinguish_strands=FALSE ) #PileupParam(distinguish_nucleotides=FALSE)
  
  range_peaks <- GenomicRanges::pintersect(search_peaks, range_plot, drop.nohit.ranges=TRUE)
  
  if( length(subset(range_peaks, sig))){
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
    
    track <- tracks_create() 
    
    range_de_peaks <- GenomicRanges::pintersect(de_peaks, range_plot, drop.nohit.ranges=TRUE) %>% 
      as.data.frame() 
    if(nrow(range_de_peaks) != 0) {
      track <- tracks_shade_bar( track, range_de_peaks, group_column='method', direction='down', range=-1 )
    }
    
    track <- tracks_annotation( track, range_gff )
    
    range_cond_peaks <- GenomicRanges::pintersect(cond_peaks, range_plot, drop.nohit.ranges=TRUE) %>% 
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
    save_name = paste0(output_dir, '/pileup', '_chr', as.character(seqnames(range_plot)), ':', start(range_plot), '-', end(range_plot), '.pdf')
    ggsave(save_name, plot=p)
    print(save_name)
    
    pmap(list(range_gff_genes$start, range_gff_genes$end, ifelse(range_gff_genes$name=='', range_gff_genes$parent_gene, range_gff_genes$name)), function(start, end, name){
      q <- p +
        xlim(start, end)
      ggsave(str_replace(save_name, '.pdf', paste0('_', name, '.pdf')), plot=q)
    })
  }
}
#})
#}
