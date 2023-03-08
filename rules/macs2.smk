from itertools import combinations

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
		pvalue = 1e-3, #Relaxed to get 'noisy' peaks for IDR
		outdir = lambda wildcards, output: os.path.dirname(output[0]),
		slocal = 2000, #Average fragment length ~1400
		prefix = lambda wildcards, output: os.path.basename(output[0]).replace('_peaks.narrowPeak', ''),
	conda: '../envs/macs2.yml'
	shell:
	 '''
		macs2 callpeak --treatment {input.ip_seq}    \
		               --control   {input.input_seq} \
		               --format    BAMPE             \
		               --gsize     hs                \
		               --keep-dup  auto              \
		               --outdir    {params.outdir}   \
		               --name      {params.prefix}   \
		               --pvalue    {params.pvalue}   \
		               --slocal    {params.slocal}   \
		               --bdg                         \
		               --SPMR                        \
		               --call-summits
	 '''

rule pooled_bams:
	input: [ os.path.join(config['align_dir'],f'{sample_id}.{config["aligner"]}_aligned.bam') for sample_id in ['{sample1_id}', '{sample2_id}'] ]
	output:
		bams  =                [ os.path.join(config["peakcalling_dir"],'idr','pooled','{sample1_id}_{sample2_id}_pooled'+f'{n}.bam') for n in [1,2] ],
		header=                  os.path.join(config["peakcalling_dir"],'idr','pooled','{sample1_id}_{sample2_id}_pooled.header'),
		splits=[ temporary(local(os.path.join(config["peakcalling_dir"],'idr','pooled','{sample1_id}_{sample2_id}_pooled_bamsplits'+f'{n}'))) for n in ['00','01'] ],
		merge =  temporary(local(os.path.join(config["peakcalling_dir"],'idr','pooled','{sample1_id}_{sample2_id}_pooled_merged.bam'))),
	params:
		split_prefix=lambda wildcards, output: output.splits[0].replace('00','')
	conda: '../envs/samtools.yml'
	shell:
	 '''
		samtools merge -u {output.merge} {input}

		samtools view -H {output.merge} > {output.header}

		linecount=$(samtools view {output.merge} | wc -l)
		linecount=$(( (linecount+1)/2 ))

		samtools view {output.merge} | shuf - | split -d -l ${{linecount}} - "{params.split_prefix}"
		cat {output.header} {output.splits[0]} | samtools view -bS - > {output.bams[0]}
		cat {output.header} {output.splits[1]} | samtools view -bS - > {output.bams[1]}
	 '''

rule pseudo_bams:
	input: os.path.join(config['align_dir'],'{sample_id}' + f'.{config["aligner"]}_aligned.bam'),
	output:
		bams  =                [ os.path.join(config["peakcalling_dir"],'idr','pseudo','{sample_id}_pseudo'+f'{n}.bam') for n in [1,2] ],
		header=                  os.path.join(config["peakcalling_dir"],'idr','pseudo','{sample_id}_pseudo.header'),
		splits=[ temporary(local(os.path.join(config["peakcalling_dir"],'idr','pseudo','{sample_id}_pseudo_bamsplits'+f'{n}'))) for n in ['00','01'] ],
	params:
		split_prefix=lambda wildcards, output: output.splits[0].replace('00','')
	conda: '../envs/samtools.yml'
	shell:
	 '''
		samtools view -H {input} > {output.header}

		linecount=$(samtools view {input} | wc -l)
		linecount=$(( (linecount+1)/2 ))

		samtools view {input} | shuf - | split -d -l ${{linecount}} - "{params.split_prefix}"
		cat {output.header} {output.splits[0]} | samtools view -bS - > {output.bams[0]}
		cat {output.header} {output.splits[1]} | samtools view -bS - > {output.bams[1]}
	 '''

use rule macs2_peakcalling as macs_peakcalling_pseudo with:
	input:
		ip_seq    =                   os.path.join(config["peakcalling_dir"],'idr','pseudo','{sample_id}_pseudo'+'{n}.bam'),
		input_seq = lambda wildcards: os.path.join(config["peakcalling_dir"],'idr','pseudo',ip_sample_id2input_sample_id[wildcards.sample_id]+'_pseudo{n}.bam'),
	output:                 multiext(os.path.join(config["peakcalling_dir"],'idr','{sample_id}_pseudo{n}'), '_peaks.narrowPeak', '_summits.bed', '_treat_pileup.bdg', '_control_lambda.bdg'),

use rule macs2_peakcalling as macs_peakcalling_pooled with:
	input:
		ip_seq    =                   os.path.join(config["peakcalling_dir"],'idr','pooled','{sample1_id}_{sample2_id}_pooled{n}.bam'),
		input_seq = lambda wildcards: os.path.join(config["peakcalling_dir"],'idr','pooled',ip_sample_id2input_sample_id[wildcards.sample1_id]+'_'+ip_sample_id2input_sample_id[wildcards.sample1_id]+'_pooled{n}.bam'),
	output:                 multiext(os.path.join(config["peakcalling_dir"],'idr','{sample1_id}_{sample2_id}_pooled{n}'), '_peaks.narrowPeak', '_summits.bed', '_treat_pileup.bdg', '_control_lambda.bdg'),

rule sort_narrow_peaks:
	input:             os.path.join(config["peakcalling_dir"], '{file}_peaks.narrowPeak'),
	output: temp(local(os.path.join(config["peakcalling_dir"], '{file}_peaks.sorted.narrowPeak'))),
	conda: '../envs/macs2.yml'
	shell: 'sort -k8,8nr {input} > {output}'

rule idr_peaks_true:
	input:
		input_peaks=lambda wildcards: [ os.path.join(config["peakcalling_dir"], f'{s}_peaks.sorted.narrowPeak') for s in [wildcards.sample1_id, wildcards.sample2_id] ],
	output:
		idr_res=os.path.join(config["peakcalling_dir"],'idr','{condition}_{sample1_id}_{sample2_id}_true.tsv'),
		idr_png=os.path.join(config["peakcalling_dir"],'idr','{condition}_{sample1_id}_{sample2_id}_true.png'),
		idr_log=os.path.join(config["peakcalling_dir"],'idr','{condition}_{sample1_id}_{sample2_id}_true.log'),
	conda: '../envs/macs2.yml'
	shell:
		'''
			idr --samples         {input_peaks}    \
			    --input-file-type narrowPeak       \
			    --rank            p.value          \
			    --output-file     {output.idr_res} \
			    --log-output-file {output.idr_log} 1
			    --plot
			mv {output.idr_res}.png {output.idr_png}
		'''

use rule idr_peaks_true as idr_peaks_pooled with:
	input:
		input_peaks=lambda wildcards: [ os.path.join(config["peakcalling_dir"], 'idr', '{sample1_id}_{sample2_id}_pooled'+f'{n}_peaks.sorted.narrowPeak') for n in [1,2] ],
	output:
		idr_res=os.path.join(config["peakcalling_dir"],'idr','{condition}_{sample1_id}_{sample2_id}_pooled.tsv'),
		idr_png=os.path.join(config["peakcalling_dir"],'idr','{condition}_{sample1_id}_{sample2_id}_pooled.png'),
		idr_log=os.path.join(config["peakcalling_dir"],'idr','{condition}_{sample1_id}_{sample2_id}_pooled.log'),

use rule idr_peaks_true as idr_peaks_pseudo with:
	input:
		input_peaks=lambda wildcards: [ os.path.join(config["peakcalling_dir"], 'idr', '{sample_id}_'+f'pseudo{n}_peaks.sorted.narrowPeak') for n in [1,2] ],
	output:
		idr_res=os.path.join(config["peakcalling_dir"],'idr','{condition}_{sample_id}_pseudo.tsv'),
		idr_png=os.path.join(config["peakcalling_dir"],'idr','{condition}_{sample_id}_pseudo.png'),
		idr_log=os.path.join(config["peakcalling_dir"],'idr','{condition}_{sample_id}_pseudo.log'),

rule idr_combination:
	input:
		idr_true  =os.path.join(config["peakcalling_dir"],'idr','{condition}_{sample1_id}_{sample2_id}_true.tsv'),
		idr_pool  =os.path.join(config["peakcalling_dir"],'idr','{condition}_{sample1_id}_{sample2_id}_pooled.tsv'),
		idr_pseudo1=os.path.join(config["peakcalling_dir"],'idr','{condition}_{sample1_id}_pseudo.tsv'),
		idr_pseudo2=os.path.join(config["peakcalling_dir"],'idr','{condition}_{sample2_id}_pseudo.tsv'),
	output:
		report=os.path.join(config["peakcalling_dir"],'idr','{condition}_{sample1_id}_{sample2_id}_report.txt')
	conda: '../envs/macs2.yml'
	shell: 'cat {input.idr_true} cat {input.idr_pool} cat {input.idr_pseudo1} cat {input.idr_pseudo2} > {output.report}'

#IDR as described in https://github.com/hbctraining/Intro-to-ChIPseq-flipped/blob/main/lessons/handling-replicates.md
rule idr_combine_condition:
	input:
		reports=lambda wildcards: [ os.path.join(config["peakcalling_dir"],'idr','{condition}_' + f'{sample1}_{sample2}_report.txt') for sample1, sample2 in combinations(condition2sample_ids[wildcards.condition], 2)],
	output:
		report=os.path.join(config["peakcalling_dir"],'idr','{condition}_report.txt'),
		
	conda: '../envs/macs2.yml'
	shell: 'cat {input} > {output}'
