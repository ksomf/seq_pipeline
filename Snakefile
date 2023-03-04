# snakemake --cores 32 --use-conda --conda-frontend mamba dev

import pandas as pd
from itertools import chain

configfile: 'config.yml'
defaults = { 'sra_dir'         : '01_sra_download'
           , 'fastq_dir'       : '02_fastq'
           , 'align_dir'       : '03_aligned'
          , 'peakcalling_dir' : '04_peakcalling'
           , 'metadata_file'   : 'metadata.tsv'
           , 'reference_dir'   : 'reference'       }
for k, v in defaults.items():
	if k not in config:
		config[k] = v

metadata = pd.read_csv(config['metadata'], sep='\t')
metadata_ip_only    = metadata[ metadata['method']=='IP' ]
metadata_input_only = metadata[ metadata['method']=='Input' ]
condition2sample_ids = { g:df['sample_id'].to_list() for g, df in metadata_ip_only   .groupby(['condition']) }
condition2input_ids  = { g:df['sample_id'].to_list() for g, df in metadata_input_only.groupby(['condition']) }
conditions = list(set(metadata['condition']))


sample_ids = metadata['sample_id']
sample_readlist = list(chain(metadata['R1'],metadata['R2']))
sample_id2reads = dict(zip(metadata['sample_id'],zip(metadata['R1'],metadata['R2'])))
ip_sample_id2input_sample_id = dict(zip(metadata['sample_id'],metadata['matching_input_control']))
sample_ids_ip = metadata[ metadata['method']=='IP' ]['sample_id']


wildcard_constraints:
	sample_id  = '|'.join(sample_ids),
	aligner    = '|'.join(['star','bowtie2']),
	condition  = '|'.join(conditions),
	condition1 = '|'.join(conditions),
	condition2 = '|'.join(conditions),

include: 'rules/assemblies.smk'
include: 'rules/bam_utils.smk'
include: 'rules/qc.smk'
include: 'rules/align_star.smk'
include: 'rules/align_bowtie2.smk'
include: 'rules/macs2.smk'
include: 'rules/pepr.smk'
include: 'rules/genrich.smk'
include: 'rules/sailor.smk'

rule all:
	input:
		peak_calls=[os.path.join(config['peakcalling_dir'], f'{sample_id}_peaks.narrowPeak') for sample_id in sample_ids_ip],
		diffbind_pepr=os.path.join(config["peakcalling_dir"],'pepr','MAVSvsd103-467_PePr_chip1_peaks.bed'),
		genrich_mavs=os.path.join(config["peakcalling_dir"],'genrich','MAVS_peaks.narrowPeak'),
		genrich_ds=os.path.join(config["peakcalling_dir"],'genrich','d103-467_peaks.narrowPeak'),

		#aligned_bam=[os.path.join(config['align_dir'], f'{sample_id}.star_aligned.bam') for sample_id in sample_ids]

#TODO: Make the blacklist use the standard chromosome names by using the chrom report to map
#TODO: Get QC working
#TODO: MACS Peakcalling
