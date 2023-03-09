rule plot_pileups:
	input:
		pepr_peaks=os.path.join(config['peakcalling_dir'],'pepr','MAVSvsd103-467_PePr_chip1_peaks.bed'),
		gff=os.path.join(config['reference_dir'], config['assembly']+'.gff'),
	output:
		plot_dir=directory(config['peakcalling_dir'],'analysis','plots'),
	conda: '../envs/pileups.yml'
	script: '../scripts/pileups.R'

