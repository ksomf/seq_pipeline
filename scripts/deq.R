library(yaml)
library(locfit)
library(tidyverse)
library(deq)

## TEST
#metadata <- read_tsv('metadata.tsv')
#metadata_condition1 <- filter( metadata, condition=='MAVS' )
#metadata_condition2 <- filter( metadata, condition=='d103-467' )
#metadata_condition1_samples_ip    <- filter( metadata_condition1, method=='IP' )   $sample_id
#metadata_condition1_samples_input <- filter( metadata_condition1, method=='Input' )$sample_id
#metadata_condition2_samples_ip    <- filter( metadata_condition2, method=='IP' )   $sample_id
#metadata_condition2_samples_input <- filter( metadata_condition2, method=='Input' )$sample_id
#
#condition1_ip_filenames    <- paste0( '03_aligned/', metadata_condition1_samples_ip   , '.star_aligned.bam' )
#condition1_input_filenames <- paste0( '03_aligned/', metadata_condition1_samples_input, '.star_aligned.bam' )
#condition2_ip_filenames    <- paste0( '03_aligned/', metadata_condition2_samples_ip   , '.star_aligned.bam' )
#condition2_input_filenames <- paste0( '03_aligned/', metadata_condition2_samples_input, '.star_aligned.bam' )
#peak_filename              <- '04_peakcalling/deq/MAVS_vs_d103-467_combined_peaks.bed'
#gtf_filename               <- '../reference/ucsc/hg38.gtf'
#
#read_length                <- 75

#output_filename <- '04_peakcalling/deq/MAVS_vs_d103-467.txt'

# SNAKEMAKE
condition1_ip_filenames    <- snakemake@input[['condition1_ips']]
condition1_input_filenames <- snakemake@input[['condition1_inputs']]
condition2_ip_filenames    <- snakemake@input[['condition2_ips']]
condition2_input_filenames <- snakemake@input[['condition2_inputs']]
peak_filename              <- snakemake@input[['peaks']]
gtf_filename               <- snakemake@input[['gtf']]
bam_infos                  <- snakemake@input[['bam_infos']]

read_length                <- snakemake@params[['read_length']]
avg_fragment_length        <- bam_infos %>%
	map_dfr( read_tsv ) %>%
	mutate( reads = `1st fragments`) %>%
	mutate( avg_insert_length   = `insert size average`) %>%
	mutate( avg_fragment_length = 2*read_length + avg_insert_length) %>%
	mutate( total_reads = sum(reads) ) %>%
	mutate( inte = (reads/total_reads) * avg_fragment_length ) %>%
	pull( inte ) %>%
	sum()

output_filename <- snakemake@output[['peaks']]

deq <- function(input.bams,ip.bams,treated.input.bams,treated.ip.bams,
				peak.files,gtf,paired.end=FALSE,outfi='deq_results.txt',
				tool='deq',compare.gene=TRUE,readlen=100,fraglen=100,nthreads=1){
	
	if (length(input.bams) != length(ip.bams) | length(treated.input.bams) != length(treated.ip.bams)){
		stop('number of IP bam files must equal number of input bam files')
	}
	n.c <- length(input.bams)
	n.t <- length(treated.input.bams)
	extension <- fraglen-readlen
	meta.data <- data.frame(Condition=c(rep('control',n.c*2),rep('treatment',n.t*2)),
							IP=c(rep('input',n.c),rep('IP',n.c),rep('input',n.t),rep('IP',n.t)))
	
	#load gtf annotations
	txdb <- GenomicFeatures::makeTxDbFromGFF(gtf,format='gtf')
	gtf.in <- rtracklayer::import(gtf)
	genenames <- unique(as.data.frame(gtf.in)[,c("gene_id","gene_name","strand")])
	anno <- list(txdb=txdb,genenames=genenames)
	annotation.order <- c("utr3","utr5","exon","intron","utr3*","utr5*")
	
	#load peaks and annotate
	peaks <- deq:::import.peaks(peak.files,anno,annotation.order) 
	
	#count reads
	all.bams <- c(input.bams,ip.bams,treated.input.bams,treated.ip.bams)
	peak.data <- peaks
	bamfiles <- all.bams
	peaks <- deq:::count.reads(peaks,all.bams,paired.end,extension,nthreads=32)
	peak.counts <- DESeq2::counts(peaks$peak.counts)
	
	#run DESeq2, edgeR, and QNB to predict changes in m6A methylation
	results <- peaks$peaks[,c("annot","main.gene")]
	results <- deq:::run.tools(results,peak.counts,meta.data,tool,input.bams,ip.bams,treated.input.bams,treated.ip.bams) 
	peaks$peak.de <- deq:::run.deseq2.4l2fc(peak.counts[,which(meta.data$IP == "IP")],
									  meta.data[which(meta.data$IP == "IP"),],'peak')
	results$peak.l2fc <- peaks$peak.de$peak.l2fc
	
	#calculate gene log2 fold change  
	if (compare.gene){
		peaks$gene.counts <- deq:::get.gene.counts(c(input.bams,treated.input.bams),
											 gtf,paired.end,extension,genenames)
		peaks$gene.de <- deq:::run.deseq2.4l2fc(peaks$gene.counts,meta.data[which(meta.data$IP == "input"),],'gene')
		results$gene.l2fc <- peaks$gene.de[results$main.gene,]$gene.l2fc
		results$gene.p <- peaks$gene.de[results$main.gene,]$gene.p
		results$gene.padj <- peaks$gene.de[results$main.gene,]$gene.padj
		results$diff.l2fc <- results$peak.l2fc - results$gene.l2fc 
	}
	
	results$start <- results$start-1
	colnames(peak.counts) <- make.names(paste0(meta.data$Condition,"_",meta.data$IP),unique=TRUE)
	write_tsv(as.data.frame(peak.counts),gsub('.txt','.counts.txt',outfi))
	write_tsv(results,outfi)
	#write.table(peak.counts,gsub('.txt','.counts.txt',outfi),quote = FALSE,sep='\t',row.names = TRUE,col.names=TRUE)
	#write.table(results,outfi,quote = FALSE,sep='\t',row.names = FALSE,col.names=TRUE)
	return(results)
}

deq( input.bams         = condition2_input_filenames
   , ip.bams            = condition2_ip_filenames
   , treated.input.bams = condition1_input_filenames
   , treated.ip.bams    = condition1_ip_filenames
   , peak.files         = peak_filename
   , gtf                = gtf_filename
   , readlen            = read_length
   , fraglen            = avg_fragment_length
   , paired.end         = TRUE
   , outfi              = output_filename            )


 input.bams         = condition2_input_filenames
 ip.bams            = condition2_ip_filenames
 treated.input.bams = condition1_input_filenames
 treated.ip.bams    = condition1_ip_filenames
 peak.files         = peak_filename
 gtf                = gtf_filename
 readlen            = read_length
 paired.end         = TRUE
 outfi              = output_filename            
