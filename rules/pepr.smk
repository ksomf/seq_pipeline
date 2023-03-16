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
