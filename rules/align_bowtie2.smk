bowtie2_index_extensions = [ f'{s}.bt2' for s in ['1','2','3','4','rev.1','rev.2'] ]

rule index_bowtie2:
	input:         os.path.join(config['reference_dir'],'{database}','{assembly}.fasta'),
	output: expand(os.path.join(config['reference_dir'],'{database}','{assembly}.{ext}'), ext=bowtie2_index_extensions, allow_missing=True),
	params:
		prefix = lambda wildcards: os.path.join(config["reference_dir"], wildcards.database, wildcards.assembly),
	conda: '../envs/align_bowtie2.yml'
	threads: workflow.cores
	shell: 'bowtie2-build --threads {threads} {input} {params.prefix}'

rule align_bowtie2:
	input:
		fastq=lambda wildcards: sample_id2fastq[wildcards.sample_id],
		index=ancient(expand(os.path.join(config['reference_dir'], config['database'], '{assembly}.{ext}'), assembly=config["assembly"], ext=bowtie2_index_extensions, allow_missing=True)),
	output: 
		bam=temp(local(os.path.join(config['align_dir'],'{sample_id}.bowtie2_aligned.unfiltered.unsorted.sam'))),
	params:
		prefix = os.path.join(config["reference_dir"], config["database"], config["assembly"]),
	conda: '../envs/align_bowtie2.yml'
	threads: 16
	shell: 'bowtie2 --threads {threads} --very-sensitive --maxins 2000 -x {params.prefix} -1 {input.fastq[0]} -2 {input.fastq[1]} -S {output}'
#			samblaster --addMateTags 2> {params.samblaster_log} | \
#			samtools sort --output-fmt bam -T {params.work_dir}/{wildcards.sample_id} -o {output.bam} - 2> {params.sort_log}
