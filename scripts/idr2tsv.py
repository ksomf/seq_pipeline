import pandas as pd
import numpy  as np

filename_idr_true    = snakemake.input['idr_true']
filename_idr_pool    = snakemake.input['idr_pool']
filename_idr_pseudo1 = snakemake.input['idr_pseudo1']
filename_idr_pseudo2 = snakemake.input['idr_pseudo2']
filename_output      = snakemake.output['report']
idr_cuttoff          = snakemake.params['idr_cuttoff']
idr_cuttoff_pooled   = snakemake.params['idr_cuttoff_pooled']
sample1_id           = snakemake.params['sample1_id']
sample2_id           = snakemake.params['sample2_id']
condition            = snakemake.params['condition']

#https://github.com/nboley/idr
idr_column_names = ['chrom', 'start', 'end', 'name', 'scaled_idr', 'strand', 'enrichment', 'p_value', 'q_value', 'merged_peak_summit', 'local_idr', 'global_idr', 's1_start', 's1_end', 's1_enrichment', 's1_summit', 's2_start', 's2_end', 's2_enrichment', 's2_summit' ]
idr_true    = pd.read_csv(filename_idr_true   , sep='\t', names=idr_column_names )
idr_pool    = pd.read_csv(filename_idr_pool   , sep='\t', names=idr_column_names )
idr_pseudo1 = pd.read_csv(filename_idr_pseudo1, sep='\t', names=idr_column_names )
idr_pseudo2 = pd.read_csv(filename_idr_pseudo2, sep='\t', names=idr_column_names )

def idr2scaled_idr(idr):
	return min(int(-125*np.log2(idr)), 1000)

idr_true   ['significant'] = idr_true   ['scaled_idr'] > idr2scaled_idr(idr_cuttoff)
idr_pool   ['significant'] = idr_pool   ['scaled_idr'] > idr2scaled_idr(idr_cuttoff_pooled)
idr_pseudo1['significant'] = idr_pseudo1['scaled_idr'] > idr2scaled_idr(idr_cuttoff)
idr_pseudo2['significant'] = idr_pseudo2['scaled_idr'] > idr2scaled_idr(idr_cuttoff)

res = [ { 'condition':condition, 'sample1':sample1_id, 'sample2':sample2_id, 'type':'true_vs_pool'      , 'n_a_sig_peaks':sum(idr_true   ['significant']), 'n_a_peaks':len(idr_true   ), 'n_b_sig_peaks':sum(idr_pool   ['significant']), 'n_b_peaks':len(idr_pool   ) }
      , { 'condition':condition, 'sample1':sample1_id, 'sample2':sample2_id, 'type':'pseudo1_vs_pseudo2', 'n_a_sig_peaks':sum(idr_pseudo1['significant']), 'n_a_peaks':len(idr_pseudo1), 'n_b_sig_peaks':sum(idr_pseudo2['significant']), 'n_b_peaks':len(idr_pseudo2) } ]
pd.DataFrame(res).to_csv(filename_output, sep='\t', index=False)
