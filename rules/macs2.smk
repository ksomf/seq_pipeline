#NOTE(kim) for atac-seq
# params.extsize=150
#		               --nomodel                         \
#		               --extsize   {params.extsize}      \
rule macs2_peakcalling:
	input:
		ip_seq    =                   os.path.join(config['align_dir'],                                      '{sample_id}' + f'.{config["aligner"]}_aligned.bam'),
		input_seq = lambda wildcards: os.path.join(config['align_dir'], ip_sample_id2input_sample_id[wildcards.sample_id]  + f'.{config["aligner"]}_aligned.bam'),
	output: multiext(config["peakcalling_dir"]+'/{sample_id}', '_peaks.narrowPeak', '_summits.bed', '_treat_pileup.bdg', '_control_lambda.bdg')
	params:
		pvalue=1e-3, #Relaxed to get 'noisy' peaks for IDR
		outdir=lambda wildcards, output: os.path.dirname(output[0]),
		slocal=2000 #Average fragment length ~1400
	conda: '../envs/macs2.yml'
	shell:
	 '''
		macs2 callpeak --treatment {input.ip_seq}        \
		               --control   {input.input_seq}     \
		               --format    BAMPE                 \
		               --gsize     hs                    \
		               --keep-dup  auto                  \
		               --outdir    {params.outdir}       \
		               --name      {wildcards.sample_id} \
		               --pvalue    {params.pvalue}       \
		               --slocal    {params.slocal}       \
		               --bdg                             \
		               --SPMR                            \
		               --call-summits
	 '''

rule sort_narrow_peaks:
	input:             os.path.join(config["peakcalling_dir"], '{file}_peaks.narrowPeak'),
	output: temp(local(os.path.join(config["peakcalling_dir"], '{file}_peaks.sorted.narrowPeak'))),
	conda: '../envs/macs2.yml'
	shell: 'sort -k8,8nr {input} > {output}'
	
#IDR as described in https://github.com/hbctraining/Intro-to-ChIPseq-flipped/blob/main/lessons/handling-replicates.md
rule idr_consitent_peaks:
	input:
		input_peaks=lambda wildcards: [ os.path.join(config["peakcalling_dir"], '{s}_peaks.sorted.narrowPeak') for s in condition2samle_ids[wildcards.condition] ],
	output:
		idr_res=os.path.join(config["peakcalling_dir"],'idr','{condition}.tsv'),
		idr_png=os.path.join(config["peakcalling_dir"],'idr','{condition}.tsv.png'),
		idr_log=os.path.join(config["peakcalling_dir"],'idr','{condition}.log'),
	conda: '../envs/macs2.yml'
	shell:
		'''
			idr --samples         {input_peaks}    \
			    --input-file-type narrowPeak       \
			    --rank            p.value          \
			    --output-file     {output.idr_res} \
			    --log-output-file {output.idr_log} 1
			    --plot
		'''
