rule macs2_peakcalling:
	input:
		ip_seq    =                   os.path.join(config('align_dir'],                                      '{sample_id}' + f'{config["aligner"]}_aligned.bam'))
		input_seq = lambda wildcards: os.path.join(config['align_dir'], id_sample_id2input_sample_id(wildcards.sample_id)  + f'{config["aligner"]}_aligned.bam')
	output: multiext(f'{config["peakcalling_dir"]/{sample_id}', '_peaks.narrowPeak', '_summits.bed', '_treat_pileup.bdg', '_control_lambda.bdg')
	params:
		outdir=lambda wildcards, output: os.path.dirname(output[0])
		extsize=150
	conda: '../envs/macs2.smk'
	shell:
	 '''
		macs2 callpeak --treatment {input.ip_seq}        \
		               --control   {input.input_seq}     \
		               --format    BAMPE                 \
		               --gsize     hs                    \
		               --keep-dup  auto                  \
		               --outdir    {params.outdir}       \
		               --name      {wildcards.sample_id} \
		               --bdg                             \
		               --SPMR                            \
		               --call-summits
	 '''
#TODO appropriate here?
#		               --nomodel                         \
#		               --extsize   {params.extsize}      \
