rule piranha_diffbind:
	input:
		ip_seq    = lambda wildcards: sample_id2bam[wildcards.sample_id]                              .replace('.bam','.bam2bed.bed'),
		input_seq = lambda wildcards: sample_id2bam[ip_sample_id2input_sample_id[wildcards.sample_id]].replace('.bam','.bam2bed.bed'),
	output: os.path.join(config["peakcalling_dir"], 'piranha', '{sample_id}.bed')
	conda: '../envs/piranha.yml'
	shell:
		'''
			Piranha -output        {output}         \
			        -bin_size_both 100              \
			        {input.ip_sep} {input.input_seq}
		'''
