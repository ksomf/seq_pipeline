rule plot_pileups:
	input:
		pepr_peaks        = [ os.path.join(config['peakcalling_dir'],'pepr'           ,f'{condition}_vs_{config["control_condition"]}__PePr_chip1_peaks.bed') for condition in config['treatment_conditions'] ],
		thor_peaks        = [ os.path.join(config['peakcalling_dir'],'thor'    ,'run' ,f'{condition}_vs_{config["control_condition"]}-diffpeaks.narrowPeak')  for condition in config['treatment_conditions'] ],
		multigps_peaks    = [ os.path.join(config['peakcalling_dir'],'multigps'       ,f'{condition}_vs_{config["control_condition"]}_diffpeaks.tsv')         for condition in config['treatment_conditions'] ],
		diffbind_peaks    =   os.path.join(config['peakcalling_dir'],'diffbind'       ,f'diffpeaks.tsv'),
		diffbind_norms    =   os.path.join(config['peakcalling_dir'],'diffbind'       ,f'norm.tsv'),
		genrich_filenames = [ os.path.join(config['peakcalling_dir'],'genrich'        ,f'{condition}_peaks.narrowPeak') for condition in conditions ],
		idr_filenames     = [ os.path.join(config['peakcalling_dir'],'idr'            ,f'{condition}_{sample1_id}_{sample2_id}_true.tsv') for condition in conditions for sample1_id, sample2_id in combinations(condition2sample_ids[condition], 2) ],
		bam_files         = [ sample_id2bam[sample_id]                            for sample_id in sample_ids ],
		bam_index_files   = [ sample_id2bam[sample_id].replace('.bam','.bam.bai') for sample_id in sample_ids ],
		gff               =   os.path.join(config['reference_dir'], config['assembly']+'.gff'),
	output:
		plot_dir     = directory(os.path.join(config['peakcalling_dir'],'analysis','plots')),
		summary_file = directory(os.path.join(config['peakcalling_dir'],'analysis','summary.txt')),
	params:
		genrich_conditions   = conditions,
		idr_conditions       = [ condition                          for condition in conditions for sample1_id, sample2_id in combinations(condition2sample_ids[condition], 2) ],
		idr_elements         = [ '_'.join([sample1_id, sample2_id]) for condition in conditions for sample1_id, sample2_id in combinations(condition2sample_ids[condition], 2) ],
		sample_ids           = sample_ids,
		treatment_conditions = config['treatment_conditions'],
		control_condition    = config['control_condition'],
		metadata             = config['metadata'],
		pepr_cuttoff         = 1e-15,
		genrich_cuttoff      = 1e-5,
		idr_cuttoff          = 0.01,
		diffbind_cuttoff     = 0.05,
		multigps_cuttoff     = 0.05,
	conda: '../envs/pileups.yml'
	script: '../scripts/pileups.R'
