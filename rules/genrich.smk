rule genrich_consitent_peaks:
	input:
		bam_ips    = lambda wildcards: [ os.path.join(config['align_dir'], f'{sample_id}.' + config["aligner"] + '_aligned.sorted_by_name.bam') for sample_id in condition2sample_ids[wildcards.condition] ],
		bam_inputs = lambda wildcards: [ os.path.join(config['align_dir'], f'{sample_id}.' + config["aligner"] + '_aligned.sorted_by_name.bam') for sample_id in condition2input_ids [wildcards.condition] ],
	output:
		narrowpeaks = os.path.join(config["peakcalling_dir"],'genrich','{condition}_peaks.narrowPeak'),
		peaks       = os.path.join(config["peakcalling_dir"],'genrich','{condition}_peaks.bed'),
		pileups     = os.path.join(config["peakcalling_dir"],'genrich','{condition}_pileups.bed'),
	conda: '../envs/genrich.yml'
	shell:
		'''
			Genrich -t "{input.bam_ips}"    \
			        -c "{input.bam_inputs}" \
			        -o {output.narrowpeaks} \
			        -k {output.pileups}     \
			        -b {output.peaks}
		'''
