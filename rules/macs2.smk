from itertools import combinations

macs2_type2pvalue = { 'normal':1e-5
                    , 'relaxed':1e-3 #Suitable to get 'noisy' peaks for IDR
                    , 'full':1 }
wildcard_constraints:
	peak_types = '|'.join(macs2_type2pvalue.keys())

#NOTE(kim) for atac-seq
# params.extsize=150
#		               --nomodel                         \
#		               --extsize   {params.extsize}      \
rule macs2_peakcalling:
	input:
		ip_seq    = lambda wildcards: sample_id2bam[wildcards.sample_id],
		input_seq = lambda wildcards: sample_id2bam[ip_sample_id2input_sample_id[wildcards.sample_id]],
	output: multiext(os.path.join(config["peakcalling_dir"], '{sample_id}_{peak_types}'), '_peaks.narrowPeak', '_peaks.xls', '_summits.bed', '_treat_pileup.bdg', '_control_lambda.bdg')
	params:
		pvalue  = lambda wildcards: macs2_type2pvalue[wildcards.peak_types],
		outdir  = lambda wildcards, output: os.path.dirname(output[0]),
		slocal  = 3000, #Average fragment length in ripseq fairly hight
		prefix  = lambda wildcards, output: os.path.basename(output[0]).replace('_peaks.narrowPeak', ''),
		fraglen = config['readlength']
	log: os.path.join(config["peakcalling_dir"], '{sample_id}_{peak_types}.log')
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
		               --call-summits > {log}
	 '''

rule pooled_bams:
	input: lambda wildcards: [ sample_id2bam[sample_id] for sample_id in [wildcards.sample1_id, wildcards.sample2_id] ]
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
	input: lambda wildcards: sample_id2bam[wildcards.sample_id],
	output:
		bams  =                [ os.path.join(config["peakcalling_dir"],'idr','pseudo','{sample_id}_pseudo'+f'{n}.bam') for n in [1,2] ],
		header=                  os.path.join(config["peakcalling_dir"],'idr','pseudo','{sample_id}_pseudo.header'),
		splits=[ temporary(local(os.path.join(config["peakcalling_dir"],'idr','pseudo','{sample_id}_pseudo_bamsplits'+f'{n}'))) for n in ['00','01'] ],
	params:
		split_prefix=lambda wildcards, output: output.splits[0].replace('00',''),
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
	output:                  multiext(os.path.join(config["peakcalling_dir"],'idr','{sample_id}_pseudo{n}_relaxed'), '_peaks.narrowPeak', '_summits.bed', '_treat_pileup.bdg', '_control_lambda.bdg'),
	log:                              os.path.join(config["peakcalling_dir"],'idr','{sample_id}_pseudo{n}_relaxed.log')

use rule macs2_peakcalling as macs_peakcalling_pooled with:
	input:
		ip_seq    =                   os.path.join(config["peakcalling_dir"],'idr','pooled','{sample1_id}_{sample2_id}_pooled{n}.bam'),
		input_seq = lambda wildcards: os.path.join(config["peakcalling_dir"],'idr','pooled',ip_sample_id2input_sample_id[wildcards.sample1_id]+'_'+ip_sample_id2input_sample_id[wildcards.sample2_id]+'_pooled{n}.bam'),
	output:                  multiext(os.path.join(config["peakcalling_dir"],'idr','{sample1_id}_{sample2_id}_pooled{n}_relaxed'), '_peaks.narrowPeak', '_summits.bed', '_treat_pileup.bdg', '_control_lambda.bdg'),
	log:                              os.path.join(config["peakcalling_dir"],'idr','{sample1_id}_{sample2_id}_pooled{n}_relaxed.log')

rule sort_narrow_peaks:
	input:             os.path.join(config["peakcalling_dir"], '{file}_peaks.narrowPeak'),
	output: temp(local(os.path.join(config["peakcalling_dir"], '{file}_peaks.sorted_head.narrowPeak'))),
	conda: '../envs/macs2.yml'
	shell: 'sort -k8,8nr {input} | head -n 100000 > {output}'

rule idr_peaks_true:
	input:
		peaks=lambda wildcards: [ os.path.join(config["peakcalling_dir"], f'{s}_relaxed_peaks.sorted_head.narrowPeak') for s in [wildcards.sample1_id, wildcards.sample2_id] ],
	output:
		idr_res=os.path.join(config["peakcalling_dir"],'idr','{condition}_{sample1_id}_{sample2_id}_true.tsv'),
		idr_png=os.path.join(config["peakcalling_dir"],'idr','{condition}_{sample1_id}_{sample2_id}_true.png'),
		idr_log=os.path.join(config["peakcalling_dir"],'idr','{condition}_{sample1_id}_{sample2_id}_true.log'),
	conda: '../envs/idr.yml'
	shell:
		'''
			idr --samples         {input.peaks}    \
			    --input-file-type narrowPeak       \
			    --rank            p.value          \
			    --output-file     {output.idr_res} \
			    --log-output-file {output.idr_log} \
			    --plot
			mv {output.idr_res}.png {output.idr_png}
		'''

use rule idr_peaks_true as idr_peaks_pooled with:
	input:
		peaks=lambda wildcards: [ os.path.join(config["peakcalling_dir"], 'idr', '{sample1_id}_{sample2_id}_pooled'+f'{n}_relaxed_peaks.sorted_head.narrowPeak') for n in [1,2] ],
	output:
		idr_res=os.path.join(config["peakcalling_dir"],'idr','{condition}_{sample1_id}_{sample2_id}_pooled.tsv'),
		idr_png=os.path.join(config["peakcalling_dir"],'idr','{condition}_{sample1_id}_{sample2_id}_pooled.png'),
		idr_log=os.path.join(config["peakcalling_dir"],'idr','{condition}_{sample1_id}_{sample2_id}_pooled.log'),

use rule idr_peaks_true as idr_peaks_pseudo with:
	input:
		peaks=lambda wildcards: [ os.path.join(config["peakcalling_dir"], 'idr', '{sample_id}_'+f'pseudo{n}_relaxed_peaks.sorted_head.narrowPeak') for n in [1,2] ],
	output:
		idr_res=os.path.join(config["peakcalling_dir"],'idr','{condition}_{sample_id}_pseudo.tsv'),
		idr_png=os.path.join(config["peakcalling_dir"],'idr','{condition}_{sample_id}_pseudo.png'),
		idr_log=os.path.join(config["peakcalling_dir"],'idr','{condition}_{sample_id}_pseudo.log'),

rule idr_combination:
	input:
		idr_true   =os.path.join(config["peakcalling_dir"],'idr','{condition}_{sample1_id}_{sample2_id}_true.tsv'),
		idr_pool   =os.path.join(config["peakcalling_dir"],'idr','{condition}_{sample1_id}_{sample2_id}_pooled.tsv'),
		idr_pseudo1=os.path.join(config["peakcalling_dir"],'idr','{condition}_{sample1_id}_pseudo.tsv'),
		idr_pseudo2=os.path.join(config["peakcalling_dir"],'idr','{condition}_{sample2_id}_pseudo.tsv'),
	output:
		report=os.path.join(config["peakcalling_dir"],'idr','{condition}_{sample1_id}_{sample2_id}_report.tsv')
	params:
		idr_cuttoff       =0.05,
		idr_cuttoff_pooled=0.01,
		sample1_id        =lambda wildcards: wildcards.sample1_id,
		sample2_id        =lambda wildcards: wildcards.sample2_id,
		condition         =lambda wildcards: wildcards.condition,
	conda: '../envs/environment.yml'
	script: '../scripts/idr2tsv.py'

#IDR as described in https://github.com/hbctraining/Intro-to-ChIPseq-flipped/blob/main/lessons/handling-replicates.md
rule idr_combine_condition:
	input:
		reports=lambda wildcards: [ os.path.join(config["peakcalling_dir"],'idr','{condition}_' + f'{sample1}_{sample2}_report.tsv') for sample1, sample2 in combinations(condition2sample_ids[wildcards.condition], 2)],
	output:
		report=os.path.join(config["peakcalling_dir"],'idr','{condition}_report.txt'),
	conda: '../envs/environment.yml'
	script: '../scripts/idr_report.py'

rule idr_combine:
	input:
		reports=[ os.path.join(config["peakcalling_dir"],'idr',f'{condition}_report.txt') for condition in conditions ]
	output:
		report=os.path.join(config["peakcalling_dir"],'idr','report.txt'),
	conda: '../envs/environment.yml'
	shell: 'cat {input} > {output}'

rule idr2tsv:
	input:
		idr_filenames = [ os.path.join(config['peakcalling_dir'],'idr',f'{condition}_{sample1_id}_{sample2_id}_true.tsv') for condition in conditions for sample1_id, sample2_id in combinations(condition2sample_ids[condition], 2) ],
	output:
		condition_peaks = os.path.join(config['peakcalling_dir'],'idr','merged_true_peaks.unnamed.tsv'),
	params:
		idr_conditions       = [ condition                          for condition in conditions for sample1_id, sample2_id in combinations(condition2sample_ids[condition], 2) ],
		idr_elements         = [ '_'.join([sample1_id, sample2_id]) for condition in conditions for sample1_id, sample2_id in combinations(condition2sample_ids[condition], 2) ],
		idr_cuttoff          = 0.01,
	run:
		idr_colnames = [ 'chrom', 'start', 'end', 'name', 'stat', 'strand', 'enrichment', 'p_value', 'q_value', 'peak', 'local_idr', 'global_idr', 's1_start', 's1_end', 's1_enrichment', 's1_summit', 's2_start', 's2_end', 's2_enrichment', 's2_summit' ]
		res = []
		condition2counts = Counter(params.idr_conditions)
		for idr_filename, idr_condition, idr_element in zip(input.idr_filenames, params.idr_conditions, params.idr_elements):
			df                = pd.read_csv(idr_filename, sep='\t', names=idr_colnames)
			df['strand']      = np.where( df['strand'] == '.', '*', df['strand'] )
			df['stat']         = 2**(df['stat']/-125)
			df['name']        = [ '_'.join(['idr', idr_condition, idr_element, str(i)]) for i,name in enumerate(df['name']) ]
			df['siblings']    = condition2counts[idr_condition]
			df['method']      = 'macs2-idr'
			df['condition']   = idr_condition
			df['significant'] = df['stat'] < params.idr_cuttoff
			df['stat_type']   = 'idr'
			res.append(df)
		res = pd.concat(res)
		res = res[['chrom', 'start', 'end', 'strand', 'name', 'method', 'condition', 'siblings', 'stat', 'stat_type', 'significant' ]]
		res.to_csv( output.condition_peaks, sep='\t', index=False )


