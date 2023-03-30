# snakemake --cores 32 --use-conda --conda-frontend mamba dev

#TODO(KIM): Should I reduce all genomics to the standard chromosomes?
#TODO(KIM): Decide if I should delete diffreps and multigps
#TODO(KIM): Make the blacklist use the standard chromosome names by using the chrom report to map


import numpy      as np
import pandas     as pd
import pybedtools as bedtools
from itertools import chain

configfile: 'config.yml'
defaults = { 'sra_dir'              : '01_sra_download'
           , 'fastq_dir'            : '02_fastq'
           , 'align_dir'            : '03_aligned'
           , 'peakcalling_dir'      : '04_peakcalling'
           , 'metadata_file'        : 'metadata.tsv'
           , 'reference_dir'        : 'reference'
           , 'database'             : 'ensembl'
           , 'treatment_conditions' : []
           , 'control_condition'    : None }
for k, v in defaults.items():
	if k not in config:
		config[k] = v
need_to_download = config['metadata_files'] == 'srr'
need_to_align    = need_to_download | (config['metadata_files'] == 'fastq')
need_to_peakcall = config['pipeline'] == 'ripseq'
need_to_stamp    = config['pipeline'] == 'stamp'

metadata = pd.read_csv(config['metadata'], sep='\t')
conditions = config['treatment_conditions'] + [config['control_condition']]
sample_ids = metadata['sample_id']

ip_sample_id2input_sample_id = dict()
metadata_ip_only    = []
metadata_input_only = []
sample_ids_ip       = []
condition2sample_ids = { g:df['sample_id'].to_list() for g, df in metadata.groupby(['condition']) }
condition2input_ids  = {}
if config['pipeline'] == 'ripseq':
	ip_sample_id2input_sample_id = dict(filter(lambda xs: xs[0] != xs[1], zip(metadata['sample_id'],metadata['matching_input_control'])))
	metadata_ip_only    = metadata[ metadata['method']=='IP' ]
	metadata_input_only = metadata[ metadata['method']=='Input' ]
	sample_ids_ip       = metadata_ip_only['sample_id']
	condition2sample_ids = { g:df['sample_id'].to_list() for g, df in metadata_ip_only.groupby(['condition']) }
	condition2input_ids  = { g:[ip_sample_id2input_sample_id[s] for s in condition2sample_ids[g]] for g in condition2sample_ids.keys() }

sample_readlist = []
sample_id2reads = {}
sample_id2bam   = { s:os.path.join(config['align_dir'], f'{s}.{config["aligner"]}_aligned.bam') for s in sample_ids }
if need_to_align:
	sample_readlist = list(chain(metadata['R1'],metadata['R2']))
	sample_id2reads = dict(zip(metadata['sample_id'],zip(metadata['R1'],metadata['R2'])))
else:
	sample_id2bam   = dict(zip(metadata['sample_id'],metadata['bam']))
print('Conditions')
print(conditions)
print('Conditin2samples')
print(condition2sample_ids)
print('Conditin2inputs')
print(condition2input_ids)
print('Sample_id2input_id')
print(ip_sample_id2input_sample_id)

multiqc_inputs = []
#generate multiqc files
if need_to_align:
	multiqc_inputs += [ os.path.join(config['fastq_dir'], os.path.basename(f).replace('.fastq.gz','_fastqc.zip' ).replace('.fq.gz','_fastqc.zip') ) for f in sample_readlist ]
	if config['aligner'] == 'star':
		multiqc_inputs += [ os.path.join(config['align_dir'], f'{sample_id}.star_aligned.unsorted.Log.final.out') for sample_id in sample_ids ]
if need_to_peakcall:
	multiqc_inputs += [ os.path.join(config["peakcalling_dir"], f'{sample_id}_full_peaks.xls') for sample_id in sample_ids_ip ]

wildcard_constraints:
	sample_id  = '|'.join(sample_ids),
	sample1_id = '|'.join(sample_ids),
	sample2_id = '|'.join(sample_ids),
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
include: 'rules/piranha.smk'
include: 'rules/pepr.smk'
include: 'rules/thor.smk'
include: 'rules/genrich.smk'
include: 'rules/diffbind.smk'
#include: 'rules/multigps.smk'
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
		plot_dir=os.path.join(config['peakcalling_dir'],'analysis','summary.txt'),
		#aligned_bam=[os.path.join(config['align_dir'], f'{sample_id}.star_aligned.bam') for sample_id in sample_ids]
		qc='multiqc_report.html',

rule dev2:
	input:
		qc='multiqc_report.html',
