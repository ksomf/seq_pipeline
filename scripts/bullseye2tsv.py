import pandas as pd

input_bedfile  = snakemake.input['bullseye']
input_genefile = snakemake.input['gene_names']
output_file    = snakemake.output[0]

bullseye_columns = [ 'chrom', 'start', 'end', 'meta', 'score', 'strand'
                   , 'control_edit_freq', 'control_coverage', 'treatment_edit_freq'
                   , 'treatment_coverage', 'edit_type', 'reps' ]
df = pd.read_csv( input_bedfile, sep='\t', names=bullseye_columns )
df_meta = pd.read_csv( input_genefile, sep='\t' )
df_split = df['meta'].str.split( '|', expand=True )
df['gene_id']          = df_split[0]
df['edit_location']    = df_split[1]
df['mutantion_events'] = df_split[3]
df['control_mutations'] = df['control_edit_freq'] * df['control_coverage']
df['treatment_mutations'] = df['treatment_edit_freq'] * df['treatment_coverage']
df = df.merge( df_meta, on='gene_id', how='left' )

df.to_csv( output_file, sep='\t', index=False )
