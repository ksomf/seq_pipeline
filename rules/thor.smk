rule thor_diffbind:
	input:
		cond1_ips          = lambda wildcards: [ os.path.join(config['align_dir'], f'{sample_id}.' + config["aligner"] + '_aligned.bam'    ) for sample_id in condition2sample_ids[wildcards.condition1] ],
		cond1_ips_index    = lambda wildcards: [ os.path.join(config['align_dir'], f'{sample_id}.' + config["aligner"] + '_aligned.bam.bai') for sample_id in condition2sample_ids[wildcards.condition1] ],
		cond1_inputs       = lambda wildcards: [ os.path.join(config['align_dir'], f'{sample_id}.' + config["aligner"] + '_aligned.bam'    ) for sample_id in condition2input_ids [wildcards.condition1] ],
		cond1_inputs_index = lambda wildcards: [ os.path.join(config['align_dir'], f'{sample_id}.' + config["aligner"] + '_aligned.bam.bai') for sample_id in condition2input_ids [wildcards.condition1] ],
		cond2_ips          = lambda wildcards: [ os.path.join(config['align_dir'], f'{sample_id}.' + config["aligner"] + '_aligned.bam'    ) for sample_id in condition2sample_ids[wildcards.condition2] ],
		cond2_ips_index    = lambda wildcards: [ os.path.join(config['align_dir'], f'{sample_id}.' + config["aligner"] + '_aligned.bam.bai') for sample_id in condition2sample_ids[wildcards.condition2] ],
		cond2_inputs       = lambda wildcards: [ os.path.join(config['align_dir'], f'{sample_id}.' + config["aligner"] + '_aligned.bam'    ) for sample_id in condition2input_ids [wildcards.condition2] ],
		cond2_inputs_index = lambda wildcards: [ os.path.join(config['align_dir'], f'{sample_id}.' + config["aligner"] + '_aligned.bam.bai') for sample_id in condition2input_ids [wildcards.condition2] ],
		chromosome_sizes   = os.path.join(config['reference_dir'],f'{config["assembly"]}_chromosome_sizes.tsv'),
		genome             = os.path.join(config['reference_dir'],f'{config["assembly"]}.fasta'),
	output:
		peaks  = os.path.join(config["peakcalling_dir"],'thor','run','{condition1}vs{condition2}-diffpeaks.bed'),
		npeaks = os.path.join(config["peakcalling_dir"],'thor','run','{condition1}vs{condition2}-diffpeaks.narrowPeaks'),
		info   = os.path.join(config["peakcalling_dir"],'thor','run','{condition1}vs{condition2}-setup.info'),
		config = os.path.join(config["peakcalling_dir"],'thor'      ,'{condition1}vs{condition2}-setup.config'), #thor fails if there are any files with the prefix in the directory
	params:
		prefix       = lambda wildcards: 'vs'.join([wildcards.condition1, wildcards.condition2]),
		output_dir   = lambda wildcards, output: os.path.dirname(output.peaks),
		cond1_nl     = lambda wildcards, input: '\n'.join(input.cond1_ips),
		cond2_nl     = lambda wildcards, input: '\n'.join(input.cond2_ips),
		cond1_inl    = lambda wildcards, input: '\n'.join(input.cond1_inputs),
		cond2_inl    = lambda wildcards, input: '\n'.join(input.cond2_inputs),
	threads: 8
	conda: '../envs/thor.yml'
	shell:
		'''
			mkdir -p {params.output_dir}
			echo '#rep1'                    >  {output.config}
			echo '{params.cond1_nl}'        >> {output.config}
			echo '#rep2'                    >> {output.config}
			echo '{params.cond2_nl}'        >> {output.config}
			echo '#inputs1'                 >> {output.config}
			echo '{params.cond1_inl}'       >> {output.config}
			echo '#inputs2'                 >> {output.config}
			echo '{params.cond2_inl}'       >> {output.config}
			echo '#chrom_sizes'             >> {output.config}
			echo '{input.chromosome_sizes}' >> {output.config}
			echo '#genome'                  >> {output.config}
			echo '{input.genome}'           >> {output.config}

			cat {output.config}

			rgt-THOR --name       {params.prefix}     \
			         --output-dir {params.output_dir} \
			         {output.config}

			'''
