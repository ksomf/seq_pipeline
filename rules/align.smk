rule count_matrix:
	input: 
		bams = [ sample_id2bam[s] for s in sample_ids ],
		gtf  = os.path.join(config['reference_dir'], f'{config["database"]}', f'{config["assembly"]}.gtf'),
	threads: lambda wildcards, input: len(input['bams']),
	output: 
		counts         = os.path.join( config['align_dir'], 'counts.tsv' ),
		counts_summary = os.path.join( config['align_dir'], 'counts.tsv.summary' ),
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
	


