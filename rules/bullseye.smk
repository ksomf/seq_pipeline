rule bullseye_parse_bam:
	input: lambda wildcards: sample_id2bam[wildcards.sample_id],
	output: os.path.join(config["stamp_dir"], '{sample_id}.matrix'),
	params:
		min_coverage = 10
	conda: '../envs/bullseye.yml'
	threads: 8
	shell: 'perl parseBAM.pl --input {input} --output {output} --cpu {threads} --minCoverage {params.min_coverage} --removeDuplicates'
