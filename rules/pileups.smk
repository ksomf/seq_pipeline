rule join_peakcall_data:
	input:
		idr_peaks     = rules.idr2tsv    .output.condition_peaks,
		genrich_peaks = rules.genrich2tsv.output.condition_peaks,
		#piranha_peaks = rules.piranha2tsv.output.condition_peaks,
	output:
		condition_peaks = os.path.join(config['peakcalling_dir'],'analysis','condition_peaks.tsv'),
	run:
		df_idr     = pd.read_csv( input.idr_peaks    , sep='\t' )
		df_genrich = pd.read_csv( input.genrich_peaks, sep='\t' )
		#df_piranha = pd.read_csv( input.piranha_peaks, sep='\t' )
		res = pd.concat([df_idr, df_genrich]) #, df_piranha])
		res.to_csv( output.condition_peaks, sep='\t', index=False )

rule join_diffbind_data:
	input:
		pepr_peaks     = rules.pepr2tsv.output.diffbind_peaks,
		thor_peaks     = rules.thor2tsv.output.diffbind_peaks,
		deq_peaks      = rules.deq2tsv .output.diffbind_peaks,
		#diffbind_peaks = rules.diffbind2tsv.output.diffbind_peaks,
	output:
		diffbind_peaks = os.path.join(config['peakcalling_dir'],'analysis','diffbind_peaks.tsv'),
	run:
		df_pepr = pd.read_csv( input.pepr_peaks, sep='\t' )
		df_thor = pd.read_csv( input.thor_peaks, sep='\t' )
		df_deq  = pd.read_csv( input.deq_peaks , sep='\t' )
		res = pd.concat([ df_pepr, df_thor, df_deq ])
		res.to_csv( output.diffbind_peaks, sep='\t', index=False )

rule get_library_sizes:
	input:  [ sample_id2bam[sample_id].replace('.bam','.bam2counts.txt') for sample_id in sample_ids ],
	output: os.path.join(config['peakcalling_dir'],'analysis','library_sizes.tsv'),
	params:
		sample_ids = sample_ids
	run:
		res = []
		for sample_id, file in zip(params.sample_ids, input):
			res.append({'sample_id':sample_id, 'library_size':int(open(file, 'r').read())})
		pd.DataFrame(res).to_csv( output[0], sep='\t', index=False )

rule plot_pileups:
	input:
		condition_peaks = rules.join_peakcall_data.output.condition_peaks,
		diffbind_peaks  = rules.join_diffbind_data.output.diffbind_peaks ,
		bam_files       = [ sample_id2bam[sample_id]                            for sample_id in sample_ids ],
		bam_index_files = [ sample_id2bam[sample_id].replace('.bam','.bam.bai') for sample_id in sample_ids ],
		library_sizes   = rules.get_library_sizes.output,
		gtf             = os.path.join(config['reference_dir'], config['database'], config['assembly']+'.gtf'),
	output:
		plot_dir     = directory(os.path.join(config['peakcalling_dir'],'analysis','plots')),
		summary_file =           os.path.join(config['peakcalling_dir'],'analysis','summary.txt'),
	params:
		treatment_conditions = config['treatment_conditions'],
		control_condition    = config['control_condition'],
		metadata             = config['metadata'],
		sample_ids           = sample_ids,
	retries: 3
	conda: '../envs/pileups.yml'
	script: '../scripts/pileups.R'
