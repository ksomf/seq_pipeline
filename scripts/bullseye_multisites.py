import pandas as pd

input_bedfile  = snakemake.input[0]
output_file    = snakemake.output[0]

df = pd.read_csv( input_bedfile, sep='\t' )
res = []
for g, d in df.groupby('gene_id'):
	if len(d) > 3:
		res.append({ 'chrom'                   : d['chrom'].to_list()[0]
		           , 'start'                   : d['start'].min()
		           , 'end'                     : d['end'].max()
		           , 'gene_id'                 : g
		           , 'gene_name'               : d['gene_name'].to_list()[0]
		           , 'mutation_sites'          : len(d)
		           , 'control_avg_coverage'    : d['control_coverage'].mean()
		           , 'control_avg_edit_freq'   : d['control_edit_freq'].mean()
		           , 'treatment_avg_coverage'  : d['treatment_coverage'].mean()
		           , 'treatment_avg_edit_freq' : d['treatment_edit_freq'].mean() })

res = pd.DataFrame(res)
res = res.sort_values( 'mutation_sites', ascending=False )
res.to_csv( output_file, sep='\t', index=False )
