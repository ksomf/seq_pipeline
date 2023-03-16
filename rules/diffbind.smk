rule diffbind_diffbind_input:
	input:
		ips    = [ sample_id2bam[sample_id]                                                      for condition in conditions for sample_id in condition2sample_ids[condition] ],
		inputs = [ sample_id2bam[sample_id]                                                      for condition in conditions for sample_id in condition2input_ids [condition] ],
		peaks  = [ os.path.join(config["peakcalling_dir"], f'{sample_id}_full_peaks.narrowPeak') for condition in conditions for sample_id in condition2sample_ids[condition] ],
	output:
		param_file = os.path.join( config['peakcalling_dir'], 'diffbind', 'params.tsv' ),
	params:
		conditions = [ condition  for condition in conditions for sample_id in condition2sample_ids[condition] ],
		ip_ids     = [ sample_ids for condition in conditions for sample_id in condition2sample_ids[condition] ],
		input_ids  = [ sample_ids for condition in conditions for sample_id in condition2input_ids [condition] ],
	run:
		df               = pd.DataFrame()
		df['SampleID'  ] = params.ip_ids
		df['Factor'    ] = params.conditions
		df['bamReads'  ] = list(input.ips)
		df['ControlID' ] = params.input_ids
		df['bamControl'] = list(input.inputs)
		df['Peaks'     ] = list(input.peaks)
		df['PeakCaller'] = 'narrow'
		df['Replicate' ] = df.groupby('Factor').rank().astype(int)

		print(df)

		df.to_csv( output.param_file, sep='\t', index=False )

rule diffbind_diffbind:
	input:
		param_file = os.path.join( config['peakcalling_dir'], 'diffbind', 'params.tsv' ),
	output:
		peaks = os.path.join( config['peakcalling_dir'], 'diffbind', 'diffpeaks.tsv' ),
		cor   = os.path.join( config['peakcalling_dir'], 'diffbind', 'cor.pdf' ),
		norm  = os.path.join( config['peakcalling_dir'], 'diffbind', 'norm.tsv' ),
	params:
		control_condition = config['control_condition']
	conda: '../envs/diffbind.yml'
	script: '../scripts/diffbind.R'
