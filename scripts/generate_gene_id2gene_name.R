library(tidyverse)

output_file <- snakemake@output[[1]]

gene_id_type <- snakemake@params[['gene_id_type']]
assembly     <- snakemake@params[['assembly']]

if( assembly == 'hg38' ){
	mart <- biomaRt::useEnsembl('ensembl', 'hsapiens_gene_ensembl')
	if( gene_id_type == 'ensembl' ){
		missing_metadata <- biomaRt::getBM( attributes=c('ensembl_gene_id', 'external_gene_name', 'gene_biotype'), mart=mart, verbose=TRUE ) %>% 
			rename(gene_name='external_gene_name', gene_id='ensembl_gene_id')
		write_tsv( missing_metadata, output_file )
	}
}
