import pandas as pd

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
