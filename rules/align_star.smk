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
		fasta=os.path.join(config['reference_dir'],'{database}','{assembly}.fasta'),
		gtf  =os.path.join(config['reference_dir'],'{database}','{assembly}.gtf'),
	output:
		index=expand(os.path.join(config['reference_dir'],'{database}','{assembly}_star_index/{f}'), f=star_index_files, allow_missing=True),
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

rule star_align_genome_shared_memory:
	input:
		index=[ os.path.join(config['reference_dir'],config['database'],f'{config["assembly"]}_star_index/{f}') for f in star_index_files ],
	output: service(os.path.join(config['align_dir'],f'config["database"]_star_index.service')),
	conda: '../envs/align_star.yml'
	params:
		index_folder=lambda wildcards, input: os.path.dirname(input['index'][0]),
		output_path=lambda wildcards, output: os.path.dirname(output[0]),
	shell:
		'''
			mkdir -p {params.output_path}

			STAR --genomeDir  {params.index_folder} \
			     --genomeLoad LoadAndExit

			touch (output}

			function unload_star() {
				STAR --genomeDir  {params.index_folder} \
				     --genomeLoad Remove
			}
			trap "unload_star" exit
		'''

rule star_align_pair_end:
	input:
		index=[ os.path.join(config['reference_dir'],config['database'],f'{config["assembly"]}_star_index/{f}') for f in star_index_files ],
		fastq=lambda w: sample_id2reads[w.sample_id],
		#shared=os.path.join(config['align_dir'],f'config["database"]_star_index.service'),
	output:
		align_chrom         = temp(local(os.path.join(config['align_dir'], '{sample_id}.star_aligned.unsorted.Aligned.out.bam'))),
		align_transcriptome = temp(local(os.path.join(config['align_dir'], '{sample_id}.star_aligned.unsorted.Aligned.toTranscriptome.out.bam'))),
		log_file            =            os.path.join(config['align_dir'], '{sample_id}.star_aligned.unsorted.Log.final.out'),
	params:
		index_folder=lambda wildcards, input: os.path.dirname(input['index'][0]),
		output_prefix=lambda wildcards, output: output['align_chrom'].replace('Aligned.out.bam',''),
		output_path=lambda wildcards, output: os.path.dirname(output['align_chrom']),
		fastq_safe=lambda wildcards, input: [ '"'+s+'"' for s in input['fastq'] ],
		overhang=config['readlength'] - 1,
	threads: 8
	conda: '../envs/align_star.yml'
	shell:
		'''
			mkdir -p {params.output_path}

			STAR --runThreadN        {threads}              \
			     --genomeDir         {params.index_folder}  \
			     --genomeLoad        LoadAndRemove          \
			     --readFilesIn       {params.fastq_safe}    \
			     --readFilesCommand  zcat                   \
			     --quantMode         TranscriptomeSAM       \
			     --outFileNamePrefix {params.output_prefix} \
			     --outSAMtype        BAM Unsorted
		'''

use rule sorted_bam as star_sorted_bam with:
	input:
		aligned_file=os.path.join(config['align_dir'], '{sample_id}.star_aligned.unsorted.Aligned.out.bam'),
	output:
		sorted_file=os.path.join(config['align_dir'], '{sample_id}.star_aligned.unfiltered.bam'),

use rule sorted_bam as star_sorted_transcriptome_bam with:
	input:
		aligned_file=os.path.join(config['align_dir'], '{sample_id}.star_aligned.unsorted.Aligned.toTranscriptome.out.bam'),
	output:
		sorted_file=os.path.join(config['align_dir'], '{sample_id}.star_aligned.transcriptome.unfiltered.bam'),

#rule rsem_bam_transcipt2genome_coords: #TODO switch star index generation to the one from rsem-prepare-reference
#	input:
#		tbam  = os.path.join(config['align_dir'], '{sample_id}.star_aligned.transcriptome.bam'),
#		fasta = os.path.join(config['reference_dir'],config['database'],'{config["assembly"]}.fasta'),
#		gtf   = os.path.join(config['reference_dir'],config['database'],'{config["assembly"]}.gtf'),
#	output: 
#		gbam = os.path.join(config['align_dir'], '{sample_id}.star_aligned.transcriptome_genome.bam'),
#	conda: '../envs/align_rsem.yml'
#	threads: 16
#	shell: 'rsem-tbam2gbam '

