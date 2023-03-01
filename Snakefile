# snakemake --cores 32 --use-conda --conda-frontend mamba dev

import pandas as pd
from itertools import chain

configfile: 'config.yml'
defaults = { 'sra_dir'       : '01_sra_download'
           , 'fastq_dir'     : '02_fastq'
           , 'align_dir'     : '03_aligned'
           , 'metadata_file' : 'metadata.tsv'
           , 'reference_dir' : 'reference'       }
for k, v in defaults.items():
	if k not in config:
		config[k] = v

metadata = pd.read_csv(config['metadata'], sep='\t')

sample_ids = metadata['sample_id']
sample_readlist = list(chain(zip(metadata['R1'],metadata['R2'])))
sample_id2reads = dict(zip(metadata['sample_id'],zip(metadata['R1'],metadata['R2'])))

wildcard_constraints:
	sample_id='|'.join(sample_ids),
	aligner  ='|'.join(['star','bowtie2']),

include: 'rules/assemblies.smk'
include: 'rules/bam_utils.smk'
include: 'rules/qc.smk'
include: 'rules/align_star.smk'
include: 'rules/align_bowtie2.smk'
include: 'rules/sailor.smk'

rule all:
	input:
		aligned_bam=[f'{config["align_dir"]}/{sample_id}.star_aligned.bam' for sample_id in sample_ids]

#TODO: Make the blacklist use the standard chromosome names by using the chrom report to map
#TODO: Get QC working
