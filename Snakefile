# snakemake --cores 32 --use-conda --conda-frontend mamba all

#TODO(KIM): Should I reduce all genomics to the standard chromosomes?

import numpy      as np
import pandas     as pd
import pybedtools as bedtools

import yaml

from itertools   import chain
from collections import Counter

configfile: 'config.yml'
defaults = { 'sra_dir'              : '01_sra_download'
           , 'fastq_dir'            : '02_fastq'
           , 'align_dir'            : '03_aligned'
           , 'peakcalling_dir'      : '04_peakcalling'
           , 'stamp_dir'            : '04_stamp'
           , 'metadata_file'        : 'metadata.tsv'
           , 'reference_dir'        : 'reference'
           , 'database'             : 'ucsc'
           , 'treatment_conditions' : []
           , 'simple_comparisons'   : []
           , 'complex_comparisons'  : {}
           , 'manual_peak_files'    : {}
           , 'control_condition'    : None }
for k, v in defaults.items():
	if k not in config:
		config[k] = v

metadata = pd.read_csv(config['metadata'], sep='\t')

# Setup helper list/maps
multiqc_inputs = [] # The files to run multiqc on

sample_ids      = metadata['sample_id']
sample_id2fastq = {} # The location of fastq files , whether given or generated
sample_id2bam   = {} # The location of the bam file, whether given or generated

conditions = list(set(metadata['condition']))
condition2sample_ids = { g[0]:df['sample_id'].to_list() for g, df in metadata.groupby(['condition']) }

# Setups based on input filetype
if   config['metadata_files'] == 'srr':
	sample_id2fastq = { sample_id:[os.path.join(config["fastq_dir"], f'{sample_id}_{read}.fastq.gz') for read in ['r1', 'r2']] for sample_id in sample_ids }
	sample_id2bam   = { sample_id:os.path.join(config["align_dir"], f'{sample_id}.{config["aligner"]}_aligned.bam') for sample_id in sample_ids }

	include: 'rules/srr.smk'
elif config['metadata_files'] == 'fastq':
	sample_id2fastq = { sample_id:[r1,r2] for sample_id,r1,r2 in zip( metadata['sample_id'], metadata['R1'], metadata['R2'] ) }
	sample_id2bam   = { sample_id:os.path.join(config["align_dir"], f'{sample_id}.{config["aligner"]}_aligned.bam') for sample_id in sample_ids }
elif config['metadata_files'] == 'bam':
	sample_id2bam   = dict(zip(metadata['sample_id'], metadata['bam']))
else:
	print(f'Unsuported filetype, config["metadata_files"]: {config["metadata_files"]}')
	exit(1)
if 'adapters' in config:
	sample_id2fastq = { k:[v.replace('fastq.gz', 'trim.fastq.gz').replace('fq.gz','.trim.fq.gz') for v in vs] for k, vs in sample_id2fastq.items() }
fastq_files = [ v for k,vs in sample_id2fastq.items() for v in vs ]

# Setup for multiqc files
multiqc_inputs += [ os.path.join(config['fastq_dir'], os.path.basename(f).replace('.fastq.gz','_fastqc.zip' ).replace('.fq.gz','_fastqc.zip') ) for k,vs in sample_id2fastq.items() for f in vs ]
if config['metadata_files'] != 'bam':
	if config['aligner'] == 'star':
		multiqc_inputs += [ os.path.join(config['align_dir'], f'{sample_id}.star_aligned.unsorted.Log.final.out') for sample_id in sample_ids ]
	elif config['aligner'] == 'bowtie2':
		print('Need to determine bowtie2 multiqc files')
	else:
		print(f'Unsuported aligner, config["aligner"]: {config["aligner"]}')
		exit(1)

# Common Rules
wildcard_constraints:
	sample_id        = '|'.join(sample_ids),
	sample1_id       = '|'.join(sample_ids),
	sample2_id       = '|'.join(sample_ids),
	aligner          = '|'.join(['star','bowtie2']),
	condition        = '|'.join(conditions),
	condition1       = '|'.join(conditions),
	condition2       = '|'.join(conditions),
	named_comparison = '|'.join(config['complex_comparisons'].keys())
include: 'rules/assemblies.smk'
include: 'rules/bam_utils.smk'
include: 'rules/qc.smk'
include: 'rules/align.smk'
include: 'rules/align_star.smk'
include: 'rules/align_bowtie2.smk'

# The pipelines
if config['pipeline'] == 'bulk':
	config['use_whitelist'] = False
	multiqc_inputs += [ os.path.join( config['align_dir'], 'counts.raw_feature_counts.tsv.summary' ) ]

	rule all:
		input:
			qc=rules.multiqc_report.output.report,
			counts=rules.count_matrix.output,

elif config['pipeline'] == 'ripseq':
	config['treatment_conditions'] = set(conditions) - set(config['control_condition'])
	ip_sample_id2input_sample_id = dict(filter(lambda xs: xs[0] != xs[1], zip(metadata['sample_id'],metadata['matching_input_control'])))
	metadata_ip_only    = metadata[ metadata['method']=='IP' ]
	metadata_input_only = metadata[ metadata['method']=='Input' ]
	sample_ids_ip       = metadata_ip_only['sample_id']
	condition2sample_ids = { g[0]:df['sample_id'].to_list() for g, df in metadata_ip_only.groupby(['condition']) }
	condition2input_ids  = { g:[ip_sample_id2input_sample_id[s] for s in condition2sample_ids[g]] for g in condition2sample_ids.keys() }
	shared_input_controls = [ sample_id for sample_id in metadata_ip_only['matching_input_control'] if len(metadata_ip_only[metadata_ip_only['matching_input_control'] == sample_id]) > 1 ]
	unique_file_bam = {}
	for s in metadata_input_only['sample_id']:
		res = sample_id2bam[s]
		if np.isin( s, shared_input_controls ):
			res = res.replace('.bam','copy.bam')
		unique_file_bam[s] = res
	print('Conditions')
	print(conditions)
	print('Conditin2samples')
	print(condition2sample_ids)
	print('Conditin2inputs')
	print(condition2input_ids)
	print('Sample_id2input_id')
	print(ip_sample_id2input_sample_id)
	multiqc_inputs += [ os.path.join(config["peakcalling_dir"], f'{sample_id}_full_peaks.xls') for sample_id in sample_ids_ip ]

	include: 'rules/macs2.smk'
	include: 'rules/piranha.smk'
	include: 'rules/pepr.smk'
	include: 'rules/deq.smk'
	include: 'rules/thor.smk'
	include: 'rules/genrich.smk'
	include: 'rules/pileups.smk'

	rule all:
		input:
			pileups=rules.plot_pileups.output.summary_file,
			qc     =rules.multiqc_report.output.report,

elif config['pipeline'] == 'stamp':
	config['use_whitelist'] = False
	conditions = list(set(metadata['condition']))
	sample_id2matching_ids = { d['sample_id']:{ condition:d[f'matching_{condition}'] for condition in conditions if condition != d['condition']} for d in metadata.to_dict(orient='records') }

	print('Conditions')
	print(conditions)
	print('Conditin2samples')
	print(condition2sample_ids)
	print('Sample_id2control_id')
	print(sample_id2matching_ids)
	print('Simple Comparisons')
	print(config['simple_comparisons'])
	print('Complex Comparisons')
	print(config['complex_comparisons'])

	include: 'rules/bullseye.smk'

	rule all:
		input:
			simple_bullseye=[os.path.join(config["stamp_dir"], f'simple_normal_condition_{condition1}_vs_{condition2}.tsv') for condition1, condition2 in config["simple_comparisons"]],
			complex_bullseye=[os.path.join(config["stamp_dir"], f'complex_normal_condition_{name}.tsv') for name in config["complex_comparisons"]],
			qc=rules.multiqc_report.output.report,
			simple_gbullseye=[os.path.join(config["stamp_dir"], f'simple_relaxed_condition_{condition1}_vs_{condition2}_edited_genes.tsv') for condition1, condition2 in config["simple_comparisons"]],
			complex_gbullseye=[os.path.join(config["stamp_dir"], f'complex_relaxed_condition_{name}_edited_genes.tsv') for name in config["complex_comparisons"]],
			plots=os.path.join(config["stamp_dir"], 'plots.flag'),
			motifs=[os.path.join(config["stamp_dir"], f'complex_normal_condition_{name}_seqlogo.svg') for name in config["complex_comparisons"]],
			homer=[os.path.join(config["stamp_dir"], f'complex_condition_{name}_homer.txt') for name in config["complex_comparisons"]],
