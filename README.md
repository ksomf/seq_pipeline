# seq_pipeline

## Running the pipeline

The pipeline is run using conda environments on the local machine using:

> snakemake --cores {cores} --resources ncbi_connection=1 --use-conda --conda-frontend mamba all

## Config File `config.yml` and Metadata File `metadata.tsv`

For examples of different config files see [parameter_templates/](https://github.com/ksomf/seq_pipeline/tree/main/parameter_templates)

## Aligners

The pipeline supports either STAR or bowtie2 for aligning given in `config.yml` as `aligner: star` or `aligner: bowtie2`

### Input Files

Currently the pipeline only acceps pair end data, except for bam file inputs.

It expects `readlength: {the_readlength}` in `config.yml` and optionally for aligning adapters can be specified in `config.yml` see [parameter_templates/adapter_trimming.yml](https://github.com/ksomf/seq_pipeline/tree/main/parameter_templates/adapter_trimming.yml), with type of adapter present beind passed on to `cutadapt`.

#### SRR

With `metadata_files: srr` in `config.yml` and `SRRxxxxxxxx` for `sample_ids` in `metadata.tsv`

#### fastq

With `metadata_files: fastq` in `config.yml` and `R1` and `R2` columns with paths to the respective fastq files in `metadata.tsv`

#### bam

With `metadata_files: bam` in `config.yml` and `bam` column with paths to the respective bam files in `metadata.tsv`

## Pipelines

### Bulk

Runs qc, and generates a count matrix, super basic ready for downstream analysis.

### Ripseq

The pipeline runs peakcalling using MACS2, and then run through various peak detection methods:

- IDR
- PePr
- DEQ
- Thor
- Genrich

before saving pileup summaries of these analysis tools. The pipeline strictly requires unique matched input control for each condtion.

The `metadata.tsv` file requires one line per sample with the following columns:

- sample_id
- condition
- method: IP or Input
- matching_input_control: the sample_id of the matching input control

The config file requires the following terms in `config.yml`:

- `control_condition`: setting which condition is compared against, all other conditions compared against this

See [parameter_templates/ripseq_config.yml](https://github.com/ksomf/seq_pipeline/blob/main/parameter_templates/ripseq_config.yml) and [parameter_templates/ripseq_metadata.tsv](https://github.com/ksomf/seq_pipeline/blob/main/parameter_templates/ripseq_metadata.tsv) for an example template

### Stamp

The pipeline runs the Bullseye C to T editing pipeline an performs analysis on the data, and plotting pileups with optionally other bam files on edit genes.

The `metadata.tsv` file requires one line per sample with the following columns:

- sample_id
- condition
- method: IP or Input
- matching_{condition}: One column per condition providing the matching samples of the same condition

The config file requires the following terms in `config.yml`:

- `complex_comparisons`: allowing comparisons of one condition against multiple conditions
- `simple_comparisons`: allowing comparisons of pairs of conditions
- `display_order`: determining the order of conditions in the plots

See [parameter_templates/stamp_config.yml](https://github.com/ksomf/seq_pipeline/blob/main/parameter_templates/stampeq_config.yml) and [parameter_templates/stamp_metadata.tsv](https://github.com/ksomf/seq_pipeline/blob/main/parameter_templates/stamp_metadata.tsv) for an example template

## Dependencies

These will be downloaded automatically by the pipeline when running using conda.

#### Base

- [snakemake](https://doi.org/10.12688/f1000research.29032.1)
- [sra-tools](https://github.com/ncbi/sra-tools)
- [parallel-fastq-dump](https://github.com/rvalieris/parallel-fastq-dump)
- [samtools](https://doi.org/10.1093/gigascience/giab008)
- [bedtools](http://bioinformatics.oxfordjournals.org/content/26/6/841.short
- [pybedtools](http://bioinformatics.oxfordjournals.org/content/27/24/3423)
- [tidyverse](doi:10.21105/joss.01686)
- [pandas](https://doi.org/10.5281/zenodo.3509134)
- [STAR](https://www.ncbi.nlm.nih.gov/pubmed/23104886)
- bowtie2
- cutadapt
- [fastqc](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [multiqc](10.1093/bioinformatics/btw354)
- [biomaRt](10.18129/B9.bioc.biomaRt)
- [GenomicRanges](10.18129/B9.bioc.GenomicRanges)
- [Rsamtools](https://bioconductor.org/packages/Rsamtools)
- [rtracklayer](doi:10.1093/bioinformatics/btp328)
- [ggseqlogo](https://doi.org/10.1093/bioinformatics/btx469)

#### STAMP

- [BULLSEYE](https://doi.org/10.1016/j.molcel.2021.12.038)

#### Ripseq

- [MACS2](https://doi.org/10.1186/gb-2008-9-9-r137)
- [IDR](https://www.jstor.org/stable/23069353)
- [piranha](https://doi.org/10.1093%2Fbioinformatics%2Fbts569)
- [genrich](https://github.com/jsh58/Genrich)
- [THOR](https://doi.org/10.1186/s12859-023-05184-5)
- [PePr](https://doi.org/10.1093/bioinformatics/btu372)
- [DEQ](https://doi.org/10.1038/s41598-020-63355-3) modified to run on more modern versions of R
