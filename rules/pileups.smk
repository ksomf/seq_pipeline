rule plot_pileups:
	input:
		pepr_peaks      =   os.path.join(config['peakcalling_dir'],'pepr','MAVSvsd103-467_PePr_chip1_peaks.bed'),
		bam_files       = [ os.path.join(config['align_dir'], f'{sample_id}.{config["aligner"]}_aligned.bam')     for sample_id in sample_ids ],
		bam_index_files = [ os.path.join(config['align_dir'], f'{sample_id}.{config["aligner"]}_aligned.bam.bai') for sample_id in sample_ids ],
		gff             =   os.path.join(config['reference_dir'], config['assembly']+'.gff'),
	output:
		plot_dir=directory(os.path.join(config['peakcalling_dir'],'analysis','plots')),
	conda: '../envs/pileups.yml'
	script: '../scripts/pileups.R'

