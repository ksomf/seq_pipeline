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
		chromosome_sizes   = os.path.join(config['reference_dir'],config['database'],f'{config["assembly"]}_chromosome_sizes.tsv'),
		genome             = os.path.join(config['reference_dir'],config['database'],f'{config["assembly"]}.fasta'),
	output:
		peaks  = os.path.join(config["peakcalling_dir"],'thor','run','{condition1}_vs_{condition2}-diffpeaks.bed'),
		npeaks = os.path.join(config["peakcalling_dir"],'thor','run','{condition1}_vs_{condition2}-diffpeaks.narrowPeak'),
		info   = os.path.join(config["peakcalling_dir"],'thor','run','{condition1}_vs_{condition2}-setup.info'),
		config = os.path.join(config["peakcalling_dir"],'thor'      ,'{condition1}_vs_{condition2}-setup.config'), #thor fails if there are any files with the prefix in the output directory
	params:
		prefix         = lambda wildcards: '_vs_'.join([wildcards.condition1, wildcards.condition2]),
		output_dir     = lambda wildcards, output: os.path.dirname(output.peaks),
		cond1_nl       = lambda wildcards, input: '\n'.join(input.cond1_ips),
		cond2_nl       = lambda wildcards, input: '\n'.join(input.cond2_ips),
		cond1_inl      = lambda wildcards, input: '\n'.join(input.cond1_inputs),
		cond2_inl      = lambda wildcards, input: '\n'.join(input.cond2_inputs),
		pvalue_cuttoff = 0.01
	threads: len(sample_ids)
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

			export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${{LD_LIBRARY_PATH:-}}" # For libpng12.so.0

			rgt-THOR --name       {params.prefix}         \
			         --output-dir {params.output_dir}     \
			         --pvalue     {params.pvalue_cuttoff} \
			         {output.config}
			'''

rule thor2tsv:
	input:
		thor_peaks = [ os.path.join( config['peakcalling_dir'], 'thor', 'run', f'{condition}_vs_{config["control_condition"]}-diffpeaks.narrowPeak' ) for condition in config['treatment_conditions'] ],
	output:
		diffbind_peaks = os.path.join(config['peakcalling_dir'],'thor','merged_peaks.tsv'),
	params:
		thor_conditions = config['treatment_conditions'],
		thor_cuttoff    = 0.01,
	run:
		broadpeak_colnames = [ 'chrom', 'start', 'end', 'name', 'always_zero', 'up_or_down', 'always_zero2', 'stat', 'always_zero3', 'always_minus_one' ]
		res = []
		for thor_filename, thor_condition in zip(input.thor_peaks, params.thor_conditions):
			df                = pd.read_csv(thor_filename, sep='\t', names=broadpeak_colnames)
			df['name']        = [ '_'.join(['thor', thor_condition, name]) for name in df['name'] ]
			df['method']      = 'thor'
			df['condition']   = thor_condition
			df['stat']        = 10**(-df['stat'])
			df['significant'] = df['stat'] < params.thor_cuttoff
			df['stat_type']   = 'pvalues'
			df['strand']      = '*'
			res.append(df)
		res = pd.concat(res)
		res = res[['chrom', 'start', 'end', 'strand', 'name', 'method', 'condition', 'stat', 'stat_type', 'significant' ]]
		res = res[res['significant']]
		res = res.sort_values('stat')
		res.to_csv(output.diffbind_peaks, sep='\t', index=False)

