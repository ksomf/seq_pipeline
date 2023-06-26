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

## Method section stub

The [Yet Another Sequencing Chewer (YASC) pipeline](https://github.com/ksomf/seq_pipeline) was run for (bulk/ripseq/stamp). 

Briefly it is a [snakemake](https://doi.org/10.12688/f1000research.29032.1) pipeline build on [samtools](), [bedtools](), [pybedtools]()

If using srr files: SRR files were downloaded using [sra-tools]() and extracted using [parallel-fastq-dump]().
If Trimming: The fastq files were trimmed using [cutadapt]() with `add_adapters_here`,

The sequences were aligned using [STAR]() or [bowtie2]() against `insert_genome` and filtered using samtools.

### Pipelines
#### Bulk

Count matrix was produced using (subread)[]

#### Ripseq

Peakcalling was performed using (MACS2)[]. Three in conditin peak finders were used: (IDR)[], (piranha)[], and (genrich)[]. Three differential peak callers were used: (DEQ)[], (THOR)[], and (PePr)[].

#### Stamp

Bam files were processed using the [bullseye]() pipeline: first bam files where parsed using `parseBAM.pl`, edit sites were found using `Find_edit_site.pl` before treatment conditions where compared to covariates using `summarise_sites.pl`. Bullseye was run both at standard settings and at relaxed settings, the latter was examined to find genes with multiple edits. Motif analysis was performed using [homer]()

- Seqlogos were generated from reference exons with edit sites using [ggseqlogo]()

### Continuing

QC was performed using [fastqc]() and [multiqc]()
Downstream analysis was performed using [tidyverse](), [pandas](), [biomaRt]()

- TODO for both ripseq and stamp: [GenomicRanges](), [Rsamtools](), [rtracklayer]()

### Citations

- Snakemake: [Mölder, F., Jablonski, K.P., Letcher, B., Hall, M.B., Tomkins-Tinch, C.H., Sochat, V., Forster, J., Lee, S., Twardziok, S.O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., Nahnsen, S., Köster, J., 2021. Sustainable data analysis with Snakemake. F1000Res 10, 33.](https://doi.org/10.12688/f1000research.29032.1)
