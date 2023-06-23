import pandas as pd

rule sorted_bam:
	input:
		aligned_file=os.path.join(config['align_dir'], '{sample_id}.{aligner}_aligned.unfiltered.unsorted.sam'),
	output:
		sorted_file=os.path.join(config['align_dir'], '{sample_id}.{aligner}_aligned.unfiltered.bam'),
	params:
		work_dir=lambda wildcards, output: os.path.dirname(output['sorted_file']),
	conda: '../envs/align.yml'
	threads: 16
	shell: 'samtools sort --output-fmt bam --threads {threads} -T {params.work_dir}/{wildcards.sample_id} -o {output.sorted_file} {input.aligned_file}'

def get_cutadapt_params(d):
	adapters = d.get('adapters', {})
	res = []
	if 'r1_adapter' in adapters:
		res += ['--adapter', adapters["r1_adapter"]]
	if 'r1_front' in adapters:
		res += ['--front'  , adapters["r1_front"]]
	if 'r1_anywhere' in adapters:
		res += ['--anywhere'  , adapters["r1_anywhere"]]
	if 'r2_adapter' in adapters:
		res += ['-A', adapters["r2_adapter"]]
	if 'r2_front' in adapters:
		res += ['-G'  , adapters["r2_front"]]
	if 'r1_anywhere' in adapters:
		res += ['-B'  , adapters["r2_anywhere"]]
	return ' '.join(res)

rule cutadapter:
	input:  expand('{file}_{read}.fastq.gz'     , read=['r1', 'r2'], allow_missing=True),
	output: expand('{file}_{read}.trim.fastq.gz', read=['r1', 'r2'], allow_missing=True),
	params:
		adapter_str=get_cutadapt_params(config)
	conda: '../envs/align.yml'
	threads: 16
	shell:
		'''
			cutadapt --cores {threads}           \
			         {params.adapter_str}        \
			         --output {output[0]}        \
			         --paired-output {output[1]} \
			         --error-rate 0.1            \
			         --minimum-length 5          \
			         {input}
		'''

use rule cutadapter as cutadapterfq with:
	input:  expand('{file}_{read}.fq.gz'     , read=['r1', 'r2'], allow_missing=True)
	output: expand('{file}_{read}.trim.fq.gz', read=['r1', 'r2'], allow_missing=True)
	

rule count_matrix_feature_counts:
	input: 
		bams = [ sample_id2bam[s] for s in sample_ids ],
		gtf  = os.path.join(config['reference_dir'], f'{config["database"]}', f'{config["assembly"]}.gtf'),
	threads: workflow.cores
	output: 
		counts         = os.path.join( config['align_dir'], 'counts.raw_feature_counts.tsv' ),
		counts_summary = os.path.join( config['align_dir'], 'counts.raw_feature_counts.tsv.summary' ),
	conda: '../envs/align.yml'
	shell:
		'''
			featureCounts -p                 \
			              --countReadPairs   \
			              -T {threads}       \
			              -t exon            \
			              -g gene_id         \
			              -a {input.gtf}     \
			              -o {output.counts} \
			              {input.bams}
		'''

rule count_matrix:
	input:  rules.count_matrix_feature_counts.output.counts,
	output: rules.count_matrix_feature_counts.output.counts.replace('raw_feature_counts.',''),
	params:
		sample_id2bams=sample_id2bam
	run:
		bam2sample_id = dict( [('Geneid','gene_id')] + [ (b,a) for a,b in params.sample_id2bams.items() ] )
		result_columns = list(bam2sample_id.values())

		df = pd.read_csv( input[0], sep='\t', header=1 )
		df = df.rename(columns=bam2sample_id)
		df = df[result_columns]
		df.to_csv( output[0], sep='\t', index=False )
