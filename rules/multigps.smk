rule multigps_diffbind_input:
	input:
		condition1_ips        = lambda wildcards: [ sample_id2bam[sample_id]                                                      for sample_id in condition2sample_ids[wildcards.condition1] ],
		condition1_inputs     = lambda wildcards: [ sample_id2bam[sample_id]                                                      for sample_id in condition2input_ids [wildcards.condition1] ],
		condition2_ips        = lambda wildcards: [ sample_id2bam[sample_id]                                                      for sample_id in condition2sample_ids[wildcards.condition2] ],
		condition2_inputs     = lambda wildcards: [ sample_id2bam[sample_id]                                                      for sample_id in condition2input_ids [wildcards.condition2] ],
	output:
		param_file = os.path.join( config['peakcalling_dir'], 'multigps', '{condition1}_vs_{condition2}_params.tsv' ),
	params:
		condition1_ips   = lambda wildcards: condition2sample_ids[wildcards.condition1],
		condition1_input = lambda wildcards: condition2input_ids [wildcards.condition1],
		condition2_ips   = lambda wildcards: condition2sample_ids[wildcards.condition2],
		condition2_input = lambda wildcards: condition2input_ids [wildcards.condition2],
	run:
		# Filename	[signal,control]	file_fmt	condition	rep

		df_ip              = pd.DataFrame()
		df_ip['filename']  = list(chain(input.condition1_ips, input.condition2_ips))
		df_ip['input']     = 'signal'
		df_ip['file_fmt']  = 'bam'
		df_ip['condition'] = [wildcards.condition1]*len(params.condition1_ips) + [wildcards.condition2]*len(params.condition2_ips)
		df_ip['rep']       = list(range(1,1+len(df_ip)))

		df_input              = pd.DataFrame()
		df_input['filename']  = list(chain(input.condition1_inputs, input.condition2_inputs))
		df_input['input']     = 'control'
		df_input['file_fmt']  = 'bam'
		df_input['condition'] = [wildcards.condition1]*len(params.condition1_input) + [wildcards.condition2]*len(params.condition2_input)
		df_input['rep']       = list(range(1,1+len(df_input)))

		df = pd.concat([ df_ip, df_input ])
		df.to_csv( output.param_file, header=False, sep='\t', index=False )

rule multigps_diffbind:
	input:
		param_file       = os.path.join( config['peakcalling_dir'], 'multigps', '{condition1}_vs_{condition2}_params.tsv' ),
		chromosome_sizes = os.path.join(config['reference_dir'],f'{config["assembly"]}_chromosome_sizes.tsv'),
	output:
		peaks = os.path.join( config['peakcalling_dir'], 'multigps', '{condition1}_vs_{condition2}_diffpeaks.tsv' ),
		cor   = os.path.join( config['peakcalling_dir'], 'multigps', '{condition1}_vs_{condition2}_cor.pdf' ),
		norm  = os.path.join( config['peakcalling_dir'], 'multigps', '{condition1}_vs_{condition2}_norm.tsv' ),
	params:
		control_condition = lambda wildcards: wildcards.condition2
	threads: len(sample_ids)
	conda: '../envs/multigps.yml'
	shell:
		'''
			java -Xmx20G -jar multigps.jar          \
			     --geninfo {input.chromosome_sizes} \
			     --threads {threads}                \
			     --design  {input.param_file}
		'''

