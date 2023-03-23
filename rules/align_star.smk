import os

star_index_files = [ 'chrLength.txt'
                   , 'chrNameLength.txt'
                   , 'chrName.txt'
                   , 'chrStart.txt'
                   , 'exonGeTrInfo.tab'
                   , 'exonInfo.tab'
                   , 'geneInfo.tab'
                   , 'Genome'
                   , 'genomeParameters.txt'
                   , 'Log.out'
                   , 'SA'
                   , 'SAindex'
                   , 'sjdbInfo.txt'
                   , 'sjdbList.fromGTF.out.tab'
                   , 'sjdbList.out.tab'
                   , 'transcriptInfo.tab' ]

rule star_generate_index:
	input:
		fasta=os.path.join(config['reference_dir'],'{assembly}.fasta'),
		gtf  =os.path.join(config['reference_dir'],'{assembly}.gtf'),
	output:
		index=expand(os.path.join(config['reference_dir'],'{assembly}_star_index/{f}'), f=star_index_files, allow_missing=True),
	params:
		index_folder=lambda wildcards, output: os.path.dirname(output['index'][0]),
		overhang=100, #maxreadlength - 1 ideally
	threads: workflow.cores
	conda: '../envs/align_star.yml'
	shell: 
		'''
			test -d {params.index_folder} && rm -r {params.index_folder}

			STAR --runThreadN       {threads}             \
			     --runMode          genomeGenerate        \
			     --genomeDir        {params.index_folder} \
			     --genomeFastaFiles {input.fasta}         \
			     --sjdbGTFfile      {input.gtf}           \
			     --sjdbOverhang     {params.overhang}
		'''

rule star_align_pair_end:
	input:
		index=expand(os.path.join(config['reference_dir'],'{assembly}_star_index/{f}'), assembly=config['assembly'], f=star_index_files, allow_missing=True),
		fastq=lambda w: sample_id2reads[w.sample_id]
	output:
		aligned_file=temp(local(os.path.join(config['align_dir'], '{sample_id}.star_aligned.unsorted.Aligned.out.bam'))),
		log_file    =           os.path.join(config['align_dir'], '{sample_id}.star_aligned.unsorted.Log.final.out'),
	params:
		index_folder=lambda wildcards, input: os.path.dirname(input['index'][0]),
		output_prefix=lambda wildcards, output: output['aligned_file'].replace('Aligned.out.bam',''),
		output_path=lambda wildcards, output: os.path.dirname(output['aligned_file']),
		overhang=100, #maxreadlength - 1 ideally
	threads: min(workflow.cores,32)
	conda: '../envs/align_star.yml'
	shell:
		'''
			mkdir -p {params.output_path}

			STAR --runThreadN        {threads}              \
			     --genomeDir         {params.index_folder}  \
			     --readFilesIn       {input.fastq}          \
			     --readFilesCommand  zcat                   \
			     --outFileNamePrefix {params.output_prefix} \
			     --outSAMtype        BAM Unsorted           \
			     --outSAMunmapped    Within
		'''

rule star_sorted_bam:
	input:
		aligned_file=os.path.join(config['align_dir'], '{sample_id}.star_aligned.unsorted.Aligned.out.bam'),
	output:
		sorted_file=os.path.join(config['align_dir'], '{sample_id}.star_aligned.unfiltered.bam'),
	params:
		work_dir=lambda wildcards, output: os.path.dirname(output['sorted_file']),
	conda: '../envs/align_star.yml'
	threads: 16
	shell: 'samtools sort --output-fmt bam --threads {threads} -T {params.work_dir}/{wildcards.sample_id} -o {output.sorted_file} {input.aligned_file}'

