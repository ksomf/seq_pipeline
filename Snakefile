# snakemake --cores 32 --use-conda --conda-frontend mamba dev

import pandas as pd

configfile: 'config.yml'

metadata = pd.read_csv(config['metadata'], sep='\t')

sample_ids = metadata['sample_id']
sample_id2reads = dict(zip(metadata['sample_id'],zip(metadata['R1'],metadata['R2'])))

wildcard_constraints:
	sample_id='|'.join(sample_ids),

include: 'rules/assemblies.smk'
include: 'rules/align_star.smk'
include: 'rules/sailor.smk'

rule all:
	input:
		aligned_bam=[f'align/{sample_id}.star_aligned.bam' for sample_id in sample_ids]
