rule join_peakcall_data:
	input:
		idr_peaks     = os.path.join(config['peakcalling_dir'],'idr'    ,'merged_true_peaks.tsv'),
		genrich_peaks = os.path.join(config['peakcalling_dir'],'genrich','merged_peaks.tsv'),
		piranha_peaks = os.path.join(config['peakcalling_dir'],'piranha','merged_peaks.tsv'),
	output:
		condition_peaks = os.path.join(config['peakcalling_dir'],'analysis','condition_peaks.tsv'),
	params:
		piranha_conditions = [ condition for condition in conditions for sample_id in condition2sample_ids[condition] ],
	run:
		df_idr     = pd.read_csv( input.idr_peaks    , sep='\t' )
		df_genrich = pd.read_csv( input.genrich_peaks, sep='\t' )
		df_piranha = pd.read_csv( input.piranha_peaks, sep='\t' )
		res = pd.concat([df_idr, df_genrich, df_piranha])
		res.to_csv( output.condition_peaks, sep='\t', index=False )

rule join_diffbind_data:
	input:
		pepr_peaks     = os.path.join(config['peakcalling_dir'],'pepr'    ,'merged_peaks.tsv'),
		thor_peaks     = os.path.join(config['peakcalling_dir'],'thor'    ,'merged_peaks.tsv'),
		diffbind_peaks = os.path.join(config['peakcalling_dir'],'diffbind','merged_peaks.tsv'),
	output:
		diffbind_peaks = os.path.join(config['peakcalling_dir'],'analysis','diffbind_peaks.tsv'),
	run:
		df_pepr = pd.read_csv( input.pepr_peaks, sep='\t' )
		df_thor = pd.read_csv( input.thor_peaks, sep='\t' )
		res = pd.concat([ df_pepr, df_thor ])
		res.to_csv( output.diffbind_peaks, sep='\t', index=False )

rule plot_pileups:
	input:
		condition_peaks = os.path.join(config['peakcalling_dir'],'analysis','condition_peaks.tsv'),
		diffbind_peaks  = os.path.join(config['peakcalling_dir'],'analysis','diffbind_peaks.tsv'),
		bam_files         = [ sample_id2bam[sample_id]                            for sample_id in sample_ids ],
		bam_index_files   = [ sample_id2bam[sample_id].replace('.bam','.bam.bai') for sample_id in sample_ids ],
		library_sizes    =   os.path.join(config['peakcalling_dir'],'diffbind',f'norm.tsv'),
		gff               =   os.path.join(config['reference_dir'], config['assembly']+'.gff'),
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
