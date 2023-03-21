rule thor_diffbind:
	input:
		cond1_ips          = lambda wildcards: [ sample_id2bam[sample_id]                            for sample_id in condition2sample_ids[wildcards.condition1] ],
		cond1_ips_index    = lambda wildcards: [ sample_id2bam[sample_id].replace('.bam','.bam.bai') for sample_id in condition2sample_ids[wildcards.condition1] ],
		cond1_inputs       = lambda wildcards: [ sample_id2bam[sample_id]                            for sample_id in condition2input_ids [wildcards.condition1] ],
		cond1_inputs_index = lambda wildcards: [ sample_id2bam[sample_id].replace('.bam','.bam.bai') for sample_id in condition2input_ids [wildcards.condition1] ],
		cond2_ips          = lambda wildcards: [ sample_id2bam[sample_id]                            for sample_id in condition2sample_ids[wildcards.condition2] ],
		cond2_ips_index    = lambda wildcards: [ sample_id2bam[sample_id].replace('.bam','.bam.bai') for sample_id in condition2sample_ids[wildcards.condition2] ],
		cond2_inputs       = lambda wildcards: [ sample_id2bam[sample_id]                            for sample_id in condition2input_ids [wildcards.condition2] ],
		cond2_inputs_index = lambda wildcards: [ sample_id2bam[sample_id].replace('.bam','.bam.bai') for sample_id in condition2input_ids [wildcards.condition2] ],
		chromosome_sizes   = os.path.join(config['reference_dir'],f'{config["assembly"]}_chromosome_sizes.tsv'),
		genome             = os.path.join(config['reference_dir'],f'{config["assembly"]}.fasta'),
	output:
		peaks  = os.path.join(config["peakcalling_dir"],'thor','run','{condition1}_vs_{condition2}-diffpeaks.bed'),
		#npeaks = os.path.join(config["peakcalling_dir"],'thor','run','{condition1}_vs_{condition2}-diffpeaks.narrowPeak'), # thor currently fails to produce this file due to internal file handle problems
		info   = os.path.join(config["peakcalling_dir"],'thor','run','{condition1}_vs_{condition2}-setup.info'),
		config = os.path.join(config["peakcalling_dir"],'thor'      ,'{condition1}_vs_{condition2}-setup.config'), #thor fails if there are any files with the prefix in the output directory
	params:
		prefix       = lambda wildcards: '_vs_'.join([wildcards.condition1, wildcards.condition2]),
		output_dir   = lambda wildcards, output: os.path.dirname(output.peaks),
		cond1_nl     = lambda wildcards, input: '\n'.join(input.cond1_ips),
		cond2_nl     = lambda wildcards, input: '\n'.join(input.cond2_ips),
		cond1_inl    = lambda wildcards, input: '\n'.join(input.cond1_inputs),
		cond2_inl    = lambda wildcards, input: '\n'.join(input.cond2_inputs),
	threads: 8
	conda: '../envs/thor.yml'
	shell:
		'''
			rm -r {params.output_dir}
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
			         {output.config} || exit 0 # This command will fail because of it's internal file handles
			'''
