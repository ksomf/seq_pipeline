require(tidyverse)

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
		mutate( label = paste0( stat_type, ':', stat ) ) %>%
		select( xmin, xmax, tab_top, tab_bot, group, alpha, label )
	
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

tracks_pileup2 <- function( obj=tracks_create(), df, norm=1, axis_label=NA, colour='pileup', prop=tibble(site=0, prop=0, .rows=0), prop_colour='prop' ){
	base_y <- -obj$width - 1
	
	# Perform simplification
	df <- df %>% 
		select( site, mean, std ) %>% 
		mutate( mean = mean / norm, std = std / norm )
	plot_df <- df %>% 
		mutate( std=ifelse(is.na(std), -Inf, std)) %>% 
		mutate( prev_mean=c(NA, head(mean, n=-1))
		      , next_mean=c(tail(mean, n=-1), NA)
		      , prev_std =c(NA, head(std , n=-1))
		      , next_std =c(tail(std , n=-1), NA) ) %>% 
		mutate( same_std       = (std == prev_std) & (std == next_std)                     # Error shade size doesn't change
		      , redundant_mean = signif((prev_mean + next_mean) / 2) == signif(mean) ) %>% # Linear relationship
		mutate( redundant_site = redundant_mean & same_std ) %>% 
		filter(!redundant_site) %>% 
		mutate( std=ifelse(std==-Inf, NA, std)) %>% 
		mutate( y=base_y+mean ) %>% 
		mutate( x=site, ymin=pmax(base_y, y-2*std), ymax=y+2*std, colour=colour, alpha=ifelse(mean==0 & next_mean==0, 0, 1) ) %>% 
		select( x, y, ymin, ymax, colour, alpha )
	
	prop_df <- prop %>% 
		mutate( colour=prop_colour ) %>% 
		mutate( x=site, ymin=base_y, ymax=base_y+prop )
	
	obj$width <- obj$width + 1 + obj$padding
	obj$tracks <- append( obj$tracks, list(list( type='pileup', df=plot_df, prop_df=prop_df, y_bot=base_y-.1, y_top=base_y-.1+1, norm=norm, axis_label=axis_label )) )
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
	
	xmin <-  Inf
	xmax <- -Inf
	pileups <- keep(tracks, map_chr(tracks, ~.$type) == 'pileup' )
	jpileups <- bind_rows(map(pileups, ~.$df))
	if( nrow(jpileups) != 0 ){
		xmin <- min( xmin, jpileups$x )
		xmax <- max( xmax, jpileups$x )
	}
		
	annotations <- keep(tracks, map_chr(tracks, ~.$type) == 'annotation' ) %>% 
		map(~.$df) %>% 
		bind_rows()
	if( nrow(annotations) != 0 ){
		xmin <- min( xmin, annotations$start )
		xmax <- max( xmax, annotations$end   )
	}
	
	exon_annotations <- keep(tracks, map_chr(tracks, ~.$type) == 'annotation' ) %>% 
		map(~.$df_exons) %>% 
		bind_rows()
	
	if(!is.null(obj$track_guides)){
		guides <- data.frame(x=xmin, xend=xmax, y=obj$track_guides, yend=obj$track_guides)
	}else{
		guides <- data.frame()
	}
	
	
	p <- ggplot()
	
	if( nrow(bars) != 0 ){
		p <- p +
			geom_rect(aes(xmin=xmin, xmax=xmax, ymin=shade_bot, ymax=shade_top, fill=group, alpha=0.1*alpha), show.legend=FALSE, data=bars) +
			geom_rect(aes(xmin=xmin, xmax=xmax, ymin=tab_bot  , ymax=tab_top  , fill=group, alpha=1.0*alpha), show.legend=TRUE , data=bars)
			#geom_text(aes(x=(xmin+xmax)/2, y=(tab_bot+tab_top)/2, label=label), colour='white'              , show.legend=FALSE, data=bars)
	}
	
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
	
	if( nrow(guides) != 0 ){
		p <- p + 
			geom_segment(aes(x=xmin ,xend=xmax, y=y, yend=y), data=guides)
	}
	
	for( pileup in pileups ){
		label_df <- data.frame( x=xmax, y=pileup$y_top - .1, label=paste0(pileup$norm), colour=pileup$df$colour      %>% unique())
		p <- p + 
			geom_path(   aes(x=x, y=y, colour=colour, alpha=alpha)             , data=pileup$df ) +
			geom_ribbon( aes(x=x, ymin=ymin, ymax=ymax, alpha=1/3, fill=colour), data=pileup$df ) +
			geom_text(   aes(x=x, y=y, label=label, colour=colour), data=label_df  )
		if( nrow(pileup$prop_df) != 0 ){
			prop_df <- data.frame( x=xmin, y=pileup$y_top - .1, label='100%', colour=pileup$prop_df$colour %>% unique())
			p <- p +
				geom_rect( aes(xmin=x-20, xmax=x+20, ymin=ymin, ymax=ymax, fill=colour ), data=pileup$prop_df ) +
				geom_text( aes(x=x, y=y, label=label, colour=colour), data=prop_df )
		}
		if(!is.null(pileup$axis_label)){
			axis_df = data.frame( x=xmin, y=(pileup$y_top+pileup$y_bot)/2, label=pileup$axis_label )
			p <- p + 
				geom_text( aes(x=x, y=y, label=label), data=axis_df, hjust='right')
		}
	}
	
	p + theme_void()
}
