import os

bowtie2_index_extensions = [ f'{s}.bt2' for s in ['1','2','3','4','rev.1','rev.2'] ]

rule index_bowtie2:
	input: os.path.join(config['reference_dir'],'{database}','{assembly}.fasta'),
	output: expand(os.path.join(config['reference_dir'],'{database}','{assembly}.{ext}'), ext=bowtie2_index_extensions, allow_missing=True)
	conda: '../envs/align_bowtie2.yml'
	threads: workflow.cores
	shell: 'bowtie2-build --threads {threads} {input} reference/{wildcards.assembly}'

rule align_bowtie2:
	input:
		fastq=lambda w: sample_id2reads[w.sample_id],
		index=expand(os.path.join(config['reference_dir'],config['database'],'{assembly}.{ext}'), ext=bowtie2_index_extensions, allow_missing=True)
	output: 
		bam=temp(local(os.path.join(config['align_dir'],'{sample_id}.bw2_aligned.unfiltered.unsorted.sam')))
	conda: '../envs/align_bowtie2.yml'
	threads: workflow.cores
	shell: 'bowtie2 --threads {threads} --very-sensitive --maxins 2000 -x reference/{assembly} -1 {input.fastq[0]} -2 {input.fastq[1]} -S {output}'
#			samblaster --addMateTags 2> {params.samblaster_log} | \
#			samtools sort --output-fmt bam -T {params.work_dir}/{wildcards.sample_id} -o {output.bam} - 2> {params.sort_log}
