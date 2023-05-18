rule genrich_consitent_peaks:
	input:
		bam_ips    = lambda wildcards: [ sample_id2bam[sample_id].replace('.bam','.sorted_by_name.bam') for sample_id in condition2sample_ids[wildcards.condition] ],
		bam_inputs = lambda wildcards: [ sample_id2bam[sample_id].replace('.bam','.sorted_by_name.bam') for sample_id in condition2input_ids [wildcards.condition] ],
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

rule genrich2tsv:
	input:
		genrich_filenames = [ os.path.join(config['peakcalling_dir'],'genrich',f'{condition}_peaks.narrowPeak') for condition in conditions ],
	output:
		condition_peaks = os.path.join(config['peakcalling_dir'],'genrich','merged_peaks.unnamed.tsv'),
	params:
		genrich_conditions = conditions,
		genrich_cuttoff    = 1e-5,
	run:
		genrich_colnames = [ 'chrom', 'start', 'end', 'name', 'scaled_auc', 'strand', 'auc', 'stat', 'qvalue', 'peak' ]
		res = []
		for genrich_filename, genrich_condition in zip(input.genrich_filenames, params.genrich_conditions):
			df                = pd.read_csv(genrich_filename, sep='\t', names=genrich_colnames)
			df['strand']      = np.where( df['strand'] == '.', '*', df['strand'] )
			df['stat']        = 10**(-df['stat'])
			df['name']        = [ '_'.join(['genrich', genrich_condition, name]) for name in df['name'] ]
			df['siblings']    = 1
			df['method']      = 'genrich'
			df['condition']   = genrich_condition
			df['significant'] = df['stat'] < params.genrich_cuttoff
			df['stat_type']   = 'pvalue'
			res.append(df)
		res = pd.concat(res)
		res = res[['chrom', 'start', 'end', 'strand', 'name', 'method', 'condition', 'siblings', 'stat', 'stat_type', 'significant']]
		res.to_csv( output.condition_peaks, sep='\t', index=False )


