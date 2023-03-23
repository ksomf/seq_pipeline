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
		df['i']          = 1
		df['Replicate' ] = df.groupby('Factor')['i'].cumsum().sub(1)
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
	threads: int(1.5*len(sample_ids))
	conda: '../envs/diffbind.yml'
	script: '../scripts/diffbind.R'

rule diffbind2tsv:
	input:
		diffbind_peaks = os.path.join(config['peakcalling_dir'],'diffbind','diffpeaks.tsv'),
	output:
		diffbind_peaks = os.path.join(config['peakcalling_dir'],'diffbind','merged_peaks.tsv'),
	params:
		diffbind_cuttoff     = 0.05,
	run:
		broadpeak_colnames = [ 'chrom', 'start', 'end', 'name', 'score', 'strand', 'enrichment', 'p_value', 'q_value' ]
		res = []
		for diffbind_filename, diffbind_condition in zip(input.diffbind_peaks, params.diffbind_conditions):
			df                = pd.read_csv(diffbind_filename, sep='\t', names=diffbind_colnames)
			df['strand']      = np.where( df['strand'] == '.', '*', df['strand'] )
			df['name']        = [ '_'.join(['diffbind', diffbind_condition, str(i)]) for i,name in enumerate(df['name']) ]
			df['method']      = 'diffbind'
			df['condition']   = diffbind_condition
			df['significant'] = df['stat'] < params.diffbind_cuttoff
			df['stat']        = 'pvalue'
			res.append(df)
		res = pd.concat(res)
		res = res[['chrom', 'start', 'end', 'strand', 'name', 'method', 'condtion', 'stat', 'stat_type', 'significant' ]]
		res.to_csv( output.diffbind_peaks, sep='\t', index=False )
