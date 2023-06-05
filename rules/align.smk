rule count_matrix:
	input: 
		bams = [ sample_id2bam[s] for s in sample_ids ],
		gtf  = os.path.join(config['reference_dir'], f'{config["database"]}', f'{config["assembly"]}.gtf'),
	output: os.path.join( config['align_dir'], 'counts.tsv' ),
	conda: '../envs/align.yml'
	shell:
		'''
			featureCounts -p                \
			              --countReadPairs  \
			              -t exon           \
			              -g gene_id        \
			              -a {input.gtf}    \
			              -o {output}       \
			              {input.bams}
		'''
	


