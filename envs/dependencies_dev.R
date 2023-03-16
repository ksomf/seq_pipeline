install.packages(c( 'ggplot2', 'ggthemes', 'tidyverse', 'furrr', 'BiocManager' ))
BiocManager::install(c( 'GenomicRanges', 'Rsamtools', 'DiffBind' )) #NOTE(KIM): Might need to be run more than once. Yay software quality!
