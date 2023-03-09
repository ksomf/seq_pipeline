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

print('Conditin2samples')
print(condition2sample_ids)
print('Conditin2inputs')
print(condition2input_ids)
print('Sample_id2input_id')
print(ip_sample_id2input_sample_id)

need_to_download = config['metadata_files'] == 'srr'
need_to_align    = need_to_download | (config['metadata_files'] == 'fastq')
need_to_peakcall = config['pipeline'] == 'ripseq'

multiqc_inputs = []
#generate multiqc files
if need_to_align:
	multiqc_inputs += [ os.path.join(config['fastq_dir'], os.path.basename(f).replace('.fastq.gz','_fastqc.zip' ) ) for f in sample_readlist ]
if need_to_align:
	multiqc_inputs += [ os.path.join(config['align_dir'], f'{sample_id}.star_aligned.unsorted.Log.final.out') for sample_id in sample_ids ]
if need_to_peakcall:
	multiqc_inputs += [ os.path.join(config["peakcalling_dir"], f'{sample_id}_peaks.xls') for sample_id in sample_ids_ip ]

wildcard_constraints:
	sample_id   = '|'.join(sample_ids),
	sample1_id  = '|'.join(sample_ids),
	sample2_id  = '|'.join(sample_ids),
	aligner     = '|'.join(['star','bowtie2']),
	condition   = '|'.join(conditions),
	condition1  = '|'.join(conditions),
	condition2  = '|'.join(conditions),

include: 'rules/assemblies.smk'
include: 'rules/bam_utils.smk'
include: 'rules/qc.smk'
include: 'rules/align_star.smk'
include: 'rules/align_bowtie2.smk'
include: 'rules/macs2.smk'
include: 'rules/pepr.smk'
include: 'rules/genrich.smk'
include: 'rules/sailor.smk'
include: 'rules/pileups.smk'

rule all:
	input:
		peak_calls=[os.path.join(config['peakcalling_dir'], f'{sample_id}_peaks.narrowPeak') for sample_id in sample_ids_ip],
		diffbind_pepr=os.path.join(config["peakcalling_dir"],'pepr','MAVSvsd103-467_PePr_chip1_peaks.bed'),
		genrich_mavs=os.path.join(config["peakcalling_dir"],'genrich','MAVS_peaks.narrowPeak'),
		genrich_ds=os.path.join(config["peakcalling_dir"],'genrich','d103-467_peaks.narrowPeak'),
		qc='multiqc_report.html',

rule dev:
	input:
		plot_dir=directory(os.path.join(config['peakcalling_dir'],'analysis','plots')),
		#aligned_bam=[os.path.join(config['align_dir'], f'{sample_id}.star_aligned.bam') for sample_id in sample_ids]

#TODO: Make the blacklist use the standard chromosome names by using the chrom report to map
#TODO: Get QC working
