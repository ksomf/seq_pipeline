rule download_sra:
	output: os.path.join(config["sra_dir"], '{sample_id}.sra')
	conda: "../envs/srr.yml"
	resources:
		ncbi_connection=1
	shell: "prefetch --output-file {output} {wildcards.sample_id}"

rule sra2fastq:
	input: rules.download_sra.output
	output: expand(os.path.join(config["fastq_dir"], '{sample_id}_{read}.fastq.gz'), read=['r1', 'r2'], allow_missing=True)
	params: 
		output_dir = config["fastq_dir"]
	conda: "../envs/srr.yml"
	threads: 8
	shell: 
		'''
			wd=$(mktemp -d -t parallel-fastq-dump-$(date +%Y-%m-%d-%H-%M-%S)-XXXXXXXXXXXXX --tmpdir={resources.tmpdir})

			parallel-fastq-dump --threads {threads} --sra-id {input} --outdir {params.output_dir} --tmpdir $wd --gzip --split-3 --dumpbase --readids --clip --read-filter pass --defline-qual '+' --defline-seq '@$ac.$si.$sg/$ri' --skip-technical

			mv $(dirname {output[0]})/{wildcards.sample_id}_pass_1.fastq.gz {output[0]}
			mv $(dirname {output[1]})/{wildcards.sample_id}_pass_2.fastq.gz {output[1]}

			rm -rf $wd
		'''
