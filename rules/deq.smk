rule deq_macs2_beds:
	input:
		peaks = lambda wildcards: [ os.path.join(config["peakcalling_dir"], f'{sample_id}_normal_peaks.narrowPeak') for sample_id in chain( condition2sample_ids[wildcards.condition1], condition2sample_ids[wildcards.condition2] ) ],
	output:
		peaks = os.path.join(config['peakcalling_dir'], 'deq', '{condition1}_vs_{condition2}_combined_peaks.bed'),
	params:
		readlength = config['readlength']
	conda: '../envs/deq.yml'
	shell: 'cat {input.peaks} | cut -f 1-6 > {output.peaks}'

rule deq_diffbind:
	input:
		condition1_ips     = lambda wildcards: [ sample_id2bam[sample_id]   for sample_id in condition2sample_ids[wildcards.condition1] ],
		condition1_inputs  = lambda wildcards: [ sample_id2bam[sample_id]   for sample_id in condition2input_ids [wildcards.condition1] ],
		condition2_ips     = lambda wildcards: [ sample_id2bam[sample_id]   for sample_id in condition2sample_ids[wildcards.condition2] ],
		condition2_inputs  = lambda wildcards: [ unique_file_bam[sample_id] for sample_id in condition2input_ids [wildcards.condition2] ],
		peaks              = os.path.join(config['peakcalling_dir'], 'deq', '{condition1}_vs_{condition2}_combined_peaks.bed'),
		gtf                = os.path.join(config['reference_dir'], config['database'], config['assembly']+'.gtf'),
	output:
		counts = os.path.join(config['peakcalling_dir'], 'deq', '{condition1}_vs_{condition2}.counts.txt'),
		peaks  = os.path.join(config['peakcalling_dir'], 'deq', '{condition1}_vs_{condition2}.txt'),
	log: os.path.join(config['peakcalling_dir'], 'deq', '{condition1}_vs_{condition2}.log'),

	conda: '../envs/deq.yml'
	script: '../scripts/deq.R'

rule deq2tsv:
	input:
		deq_peaks  = [ os.path.join(config['peakcalling_dir'],'deq',f'{condition}_vs_{config["control_condition"]}.txt')        for condition in config['treatment_conditions'] ],
		deq_counts = [ os.path.join(config['peakcalling_dir'],'deq',f'{condition}_vs_{config["control_condition"]}.counts.txt') for condition in config['treatment_conditions'] ],
	output:
		diffbind_peaks = os.path.join(config['peakcalling_dir'], 'deq', 'merged_peaks.tsv')
	params:
		deq_conditions     = config['treatment_conditions'],
		deq_cuttoff        = 0.1,
		deq_lfc_cuttoff    = 1,
		deq_peak_min_count = 10,
	run:
		res = []
		for deq_filename, count_filename, deq_condition in zip(input.deq_peaks, input.deq_counts, params.deq_conditions):
			df                = pd.read_csv( deq_filename  , sep='\t' )
			df2               = pd.read_csv( count_filename, sep='\t' )
			print(df)
			df['chrom']       = df['seqnames']
			df['name']        = [ '_'.join(['deq', deq_condition, str(gn), str(s), str(e)]) for gn, s, e in zip(df['main.gene'],df['start'],df['end']) ]
			df['method']      = 'deq'
			df['condition']   = deq_condition

			# Filtered from https://www.biorxiv.org/content/10.1101/657130v1.full.pdf
			# One of the tools detect it and under that peak there is good lfc
			df['min_padj']    = np.nanmin(df[['deseq2.padj', 'edger.padj', 'qnb.padj']], axis=1) 
			df['tool_sig']    = df['min_padj'] < params.deq_cuttoff
			df['peak_sig']    = (np.abs(df['peak.l2fc']) > params.deq_lfc_cuttoff) # & (np.abs(df['gene.l2fc']) > params.deq_lfc_cuttoff)
			df['count_sig']   = np.sum(df2 / len(df2.columns) > params.deq_peak_min_count, axis=1)
			df['significant'] = df['tool_sig'] & df['peak_sig'] & df['count_sig']
			df['stat']        = df['min_padj']
			df['stat_type']   = 'min_adj_pvalue'
			res.append(df)
		res = pd.concat(res)
		res = res[['chrom', 'start', 'end', 'strand', 'name', 'method', 'condition', 'stat', 'stat_type', 'significant' ]]
		res = res.sort_values('stat')
		res.to_csv( output.diffbind_peaks, sep='\t', index=False )

