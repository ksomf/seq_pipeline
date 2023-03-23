install.packages(c( 'ggplot2', 'ggthemes', 'tidyverse', 'furrr', 'BiocManager', 'devtools' ))
BiocManager::install(c( 'GenomicRanges', 'Rsamtools', 'DiffBind', 'DESeq2', 'edgeR', 'GenomicFeatures' )) #NOTE(KIM): Might need to be run more than once. Yay software quality!




# DEQ
BiocManager::install(c( 'GenomicRanges', 'Rsamtools', 'DESeq2', 'edgeR', 'GenomicFeatures', 'Rsubread' ))

install.packages('https://www.bioconductor.org/packages//2.13/bioc/src/contrib/exomePeak_1.0.0.tar.gz', repos=NULL, type='source')
install.packages('https://www.bioconductor.org/packages//2.12/bioc/src/contrib/DESeq_1.12.1.tar.gz'   , repos=NULL, type='source')

devtools::install_github("lianliu09/QNB")
devtools::install_github("al-mcintyre/DEQ")
