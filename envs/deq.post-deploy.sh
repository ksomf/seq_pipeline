#!env bash
set -o pipefail

R -e "install.packages('https://www.bioconductor.org/packages//2.13/bioc/src/contrib/exomePeak_1.0.0.tar.gz', repos=NULL, type='source')"
R -e "install.packages('https://www.bioconductor.org/packages//2.12/bioc/src/contrib/DESeq_1.12.1.tar.gz'   , repos=NULL, type='source')"

R -e "devtools::install_github('lianliu09/QNB')"
R -e "devtools::install_github('al-mcintyre/DEQ')"
