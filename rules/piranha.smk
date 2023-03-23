rule piranha_diffbind:
	input:
		ip_seq    = lambda wildcards: sample_id2bam[wildcards.sample_id]                              .replace('.bam','.bam2bed.weirdsort.bed'),
		input_seq = lambda wildcards: sample_id2bam[ip_sample_id2input_sample_id[wildcards.sample_id]].replace('.bam','.bam2bed.weirdsort.bed'),
	output: os.path.join(config["peakcalling_dir"], 'piranha', '{sample_id}.bed')
	conda: '../envs/piranha.yml'
	shell:
		'''
			Piranha -output        {output}         \
			        -bin_size_both 100              \
			        {input.ip_seq} {input.input_seq}
		'''

from collections import Counter
rule piranha2tsv:
	input:
		piranha_peaks = [ os.path.join(config["peakcalling_dir"], 'piranha', f'{sample_id}.bed') for condition in conditions for sample_id in condition2sample_ids[condition] ],
	output:
		diffbind_peaks =   os.path.join(config['peakcalling_dir'],'piranha','merged_peaks.tsv'),
	params:
		pirahna_conditions = [ condition for condition in conditions for sample_id in condition2sample_ids[condition] ],
	run:
		piranha_colnames = [ 'chrom', 'start', 'end', 'name', 'score', 'strand', 'stat', 'unknown' ]
		res = []
		condition2counts = Counter(params.idr_conditions)
		for piranha_filename, piranha_condition in zip(input.piranha_peaks, params.piranha_conditions):
			df                = pd.read_csv(piranha_filename, sep='\t', names=piranha_colnames)
			df['name']        = [ '_'.join(['piranha', piranha_condition, i]) for i,name in enumerate(df['name']) ]
			df['method']      = 'piranha'
			df['siblings']    = condition2counts[piranha_condition]
			df['condition']   = piranha_condition
			df['significant'] = True #df['stat'] < params.piranha_cuttoff #TODO(KIM): Need further cuttoff?
			df['stat_type']   = 'pvalue'
			res.append(df)
		res = pd.concat(res)
		res = res[['chrom', 'start', 'end', 'strand', 'name', 'method', 'condtion', 'siblings', 'stat', 'stat_type', 'significant' ]]
		res.to_csv( output.diffbind_peaks, sep='\t', index=False )
