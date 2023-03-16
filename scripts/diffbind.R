library(DiffBind)
library(tidyverse)

# Baisically remove from bam files reads from contigs not in their default detected chromosome file
# https://support.bioconductor.org/p/9138617/#9148413
pv.countGreylistEdited <- function (bamfile, pv, ktype) {
  #gl <- new("GreyList", karyotype = ktype[pv$chrmap, ]) #previous line
  #edit to restrict to just chromosomes included in the ktype object
  gl <- new("GreyList", karyotype = ktype[intersect(pv$chrmap, names(ktype)), ])
  gl <- GreyListChIP::countReads(gl, bamfile)
  return(gl)
}
environment(pv.countGreylistEdited) <- asNamespace('DiffBind')
assignInNamespace("pv.countGreylist", pv.countGreylistEdited, ns = "DiffBind")

input_param_file  <- snakemake@input[['param_file']]
output_peak_file  <- snakemake@output[['peaks']]
output_cor_file   <- snakemake@output[['cor']]
output_norm_file  <- snakemake@output[['norm']]
control_condition <- snakemake@params[['condition2']]
threads           <- snakemake@threads

#input_param_file  <- '04_peakcalling/diffbind/MAVS_vs_d103-467_params.tsv'
#output_cor_file   <- '04_peakcalling/diffbind/MAVS_vs_d103-467_cor.pdf'
#output_peak_file  <- '04_peakcalling/diffbind/MAVS_vs_d103-467_diffpeaks.tsv'
#output_norm_file  <- '04_peakcalling/diffbind/MAVS_vs_d103-467_norms.tsv'
#control_condition <- 'd103-467'
#threads           <- 8

options(mc.cores = threads)

input_params <- read_tsv(input_param_file)
db <- dba(sampleSheet=input_params)

db <- dba.count(db)
pdf(output_cor_file)
plot(db)
dev.off()

db <- dba.normalize(db)
norm <- dba.normalize(db, bRetrieve=TRUE)
norm

df_norms <- data.frame( sample_id=input_params$SampleID, norm_factor=norm$norm.factors, lib_size=norm$lib.sizes )
write_tsv(df_norms, output_norm_file)

db <- dba.contrast(db, reorderMeta=list(Factor=control_condition))
db <- dba.analyze(db)

#plot(db, contrast=1)

df <- dba.report(db)
write_tsv(df, output_peak_file)