rule pepr_diffbind:
	input:
		cond1_ips        = lambda wildcards: [ sample_id2bam[sample_id].replace('.bam','.sorted_by_name.bam') for sample_id in condition2sample_ids[wildcards.condition1] ],
		cond1_inputs     = lambda wildcards: [ sample_id2bam[sample_id].replace('.bam','.sorted_by_name.bam') for sample_id in condition2input_ids [wildcards.condition1] ],
		cond2_ips        = lambda wildcards: [ sample_id2bam[sample_id].replace('.bam','.sorted_by_name.bam') for sample_id in condition2sample_ids[wildcards.condition2] ],
		cond2_inputs     = lambda wildcards: [ sample_id2bam[sample_id].replace('.bam','.sorted_by_name.bam') for sample_id in condition2input_ids [wildcards.condition2] ],
	output:
		pepr_cond1  = os.path.join(config["peakcalling_dir"],'pepr','{condition1}_vs_{condition2}__PePr_chip1_peaks.bed'),
		pepr_cond2  = os.path.join(config["peakcalling_dir"],'pepr','{condition1}_vs_{condition2}__PePr_chip2_peaks.bed'),
		pepr_params = os.path.join(config["peakcalling_dir"],'pepr','{condition1}_vs_{condition2}__PePr_parameters.txt'),
	params:
		prefix       = lambda wildcards: '_vs_'.join([wildcards.condition1, wildcards.condition2]),
		output_dir   = lambda wildcards, output: os.path.dirname(output.pepr_cond1),
		cond1_ips    = lambda wildcards, input: ','.join(input.cond1_ips),
		cond1_inputs = lambda wildcards, input: ','.join(input.cond1_inputs),
		cond2_ips    = lambda wildcards, input: ','.join(input.cond2_ips),
		cond2_inputs = lambda wildcards, input: ','.join(input.cond2_inputs),
	threads: 8
	conda: '../envs/pepr.yml'
	shell:
		'''
			PePr --input1           {params.cond1_inputs} \
			     --chip1            {params.cond1_ips}    \
			     --input2           {params.cond2_inputs} \
			     --chip2            {params.cond2_ips}    \
			     --name             {params.prefix}       \
			     --file-format      bampe                 \
			     --peaktype         sharp                 \
			     --num-processors   {threads}             \
			     --output-directory {params.output_dir}   \
			     --diff
		'''

rule pepr2tsv:
	input:
		pepr_peaks = [ os.path.join(config['peakcalling_dir'],'pepr',f'{condition}_vs_{config["control_condition"]}__PePr_chip1_peaks.bed') for condition in config['treatment_conditions'] ],
	output:
		diffbind_peaks = os.path.join(config['peakcalling_dir'],'pepr','merged_peaks.tsv'),
	params:
		pepr_conditions = config['treatment_conditions']
		pepr_cuttoff = 1e-15,
	run:
		broadpeak_colnames <- [ 'chrom', 'start', 'end', 'name', 'score', 'strand', 'enrichment', 'p_value', 'q_value' ]
		res = []
		for pepr_filename, pepr_condition in zip(input.pepr_peaks, params.pepr_conditions):
			df                = pd.read_csv(pepr_filename, sep='\t', names=pepr_colnames)
			df['strand']      = np.where( df['strand'] == '.', '*', df['strand'] )
			df['name']        = [ '_'.join(['pepr', pepr_condition, i]) for i,name in enumerate(df['name']) ]
			df['method']      = 'pepr'
			df['condition']   = pepr_condition
			df['significant'] = df['stat'] < params.pepr_cuttoff
			df['stat']        = 'pvalue'
			res.append(df)
		res = pd.concat(res)
		res = res[['chrom', 'start', 'end', 'strand', 'name', 'method', 'condtion', 'stat', 'stat_type', 'significant' ]]
		res.to_csv( output.condition_peaks, sep='\t', index=False )
