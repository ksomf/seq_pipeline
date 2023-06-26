import pandas as pd

from itertools import combinations

report_filenames = snakemake.input['reports']
output_filename  = snakemake.output['report']

df = pd.concat([ pd.read_csv(f, sep='\t') for f in report_filenames ])

condition = list(set(df['condition']))[0]

res = [ f'IDR reprot for {condition}'
      ,  '	true against pooled '    ]

for d in df[ df['type'] == 'true_vs_pool' ].to_dict(orient='records'):
	prop = d['n_b_sig_peaks'] / d['n_a_sig_peaks']
	res.append(f'		{d["sample1"]} {d["sample2"]}	{prop}		(true {d["n_a_sig_peaks"]}, pooled {d["n_b_sig_peaks"]})')


res.append('	pseudo replicates')
for d in df[ df['type'] == 'pseudo1_vs_pseudo2' ].to_dict(orient='records'):
	n_smaller = min(d['n_b_sig_peaks'], d['n_a_sig_peaks'])
	n_larger  = max(d['n_b_sig_peaks'], d['n_a_sig_peaks'])
	prop = n_larger / n_smaller
	res.append(f'		{d["sample1"]} {d["sample2"]}	{prop}		(p1 {d["n_a_sig_peaks"]}, p2 {d["n_b_sig_peaks"]})')

res.append('')
open(output_filename, 'w+').write('\n'.join(res))
