rule annotate_peaks:
	input:
		diffbind_peaks  = os.path.join('{path}','{name}.unnamed.tsv'),
		gtf             = os.path.join(config['reference_dir'], config['database'], config['assembly']+'.gtf'),
		gene_conversion = os.path.join(config['reference_dir'], config['database'], config['assembly'] + '_gene_id2gene_name.tsv')
	output:
		diffbind_peaks = os.path.join('{path}','{name}.tsv')
	conda: '../envs/deq.yml'
	script: '../scripts/annotate_peaks.R'

rule join_peakcall_data:
	input:
		idr_peaks     = rules.idr2tsv    .output.condition_peaks.replace('.unnamed.tsv','.tsv'),
		genrich_peaks = rules.genrich2tsv.output.condition_peaks.replace('.unnamed.tsv','.tsv'),
		#piranha_peaks = rules.piranha2tsv.output.condition_peaks.replace('.unnamed.tsv','.tsv'),
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
		pepr_peaks     = rules.pepr2tsv.output.diffbind_peaks.replace('.unnamed.tsv','.tsv'),
		thor_peaks     = rules.thor2tsv.output.diffbind_peaks.replace('.unnamed.tsv','.tsv'),
		deq_peaks      = rules.deq2tsv .output.diffbind_peaks,
		#diffbind_peaks = rules.diffbind2tsv.output.diffbind_peaks.replace('.unnamed.tsv','.tsv'),
	output:
		diffbind_peaks = os.path.join(config['peakcalling_dir'],'analysis','diffbind_peaks.tsv'),
	run:
		df_pepr = pd.read_csv( input.pepr_peaks, sep='\t' )
		df_thor = pd.read_csv( input.thor_peaks, sep='\t' )
		df_deq  = pd.read_csv( input.deq_peaks , sep='\t' )
		res = pd.concat([ df_pepr, df_thor, df_deq ])
		res.to_csv( output.diffbind_peaks, sep='\t', index=False )

rule generate_manual_data:
	output:
		manual_peaks = os.path.join(config['peakcalling_dir'],'analysis','manual_peaks.tsv')
	run:
		#convert to peaks
		columns = [ 'chrom', 'start', 'end', 'strand', 'name', 'method', 'condition', 'stat', 'stat_type', 'significant' ]
		res = []
		for peak_condition, peak_file in config['manual_peak_files'].items():
			peaks = pd.read_csv( peak_file, sep='\t' )
			peaks['condition'] = peak_condition
			res.append( peaks )
		res                 = pd.concat(res)
		res                 = res[['chrom', 'start', 'end', 'condition', 'gene_id', 'gene_name']]
		res['strand']       = '*'
		res['name']         = [ f'manual_{i}' for i in range(len(res)) ]
		res['method']       = 'manual'
		res['stat']         = 0
		res['stat_type']    = 'manual'
		res['significant']  = True
		res['annot']        = 'none (bullseye)'
		res['gene_biotype'] = 'none (bullseye)'
		res.to_csv( output.manual_peaks, sep='\t', index=False )

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
		manual_peaks    = rules.generate_manual_data.output.manual_peaks ,
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
	conda: '../envs/pileups.yml'
	script: '../scripts/pileups.R'
