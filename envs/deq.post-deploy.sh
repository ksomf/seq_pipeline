#!env bash
set -o pipefail

#R -e "install.packages('S4Vectors', repos='https://cloud.r-project.org')"
R -e "install.packages('https://www.bioconductor.org/packages//2.13/bioc/src/contrib/exomePeak_1.0.0.tar.gz', repos=NULL, type='source')"
R -e "install.packages('https://www.bioconductor.org/packages//2.12/bioc/src/contrib/DESeq_1.12.1.tar.gz'   , repos=NULL, type='source')"

#R -e "devtools::install_github('lianliu09/QNB')" # This would only install 0.99 which doesn't by default output p.adj values breaking DEQ
R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/QNB/QNB_1.1.11.tar.gz', repos = NULL, type='source')"
R -e "devtools::install_github('ksomf/DEQ')"
