ule bullseye_parse_bam:
	input: 
		bam = lambda wildcards: sample_id2bam[wildcards.sample_id],
		bai = lambda wildcards: sample_id2bam[wildcards.sample_id] + '.bai',
	output:
		matrix = os.path.join(config["stamp_dir"], '{sample_id}.matrix.gz'),
		tbi    = os.path.join(config["stamp_dir"], '{sample_id}.matrix.gz.tbi'),
	params:
		min_coverage = 10,
		file_prefix = lambda wildcards, output: output['matrix'].replace('.gz','')
	conda: '../envs/bullseye.yml'
	threads: 8
	shell: 'parseBAM.pl --input {input.bam} --output {params.file_prefix} --cpu {threads} --minCoverage {params.min_coverage} --removeDuplicates'

rule bullseye_generate_genepresd:
	input:  os.path.join(config["reference_dir"], '{database}', '{assembly}.gtf'),
	output: os.path.join(config["reference_dir"], '{database}', '{assembly}.refFlat'),
	conda: '../envs/bullseye.yml'
	shell: 'perl $(which gtf2genepred.pl) --gtf {input} --out {output}'

bullseye_params = { 'normal'  : { 'minimum_rate':5, 'edit_threshold':1.5 }
                  , 'relaxed' : { 'minimum_rate':2, 'edit_threshold':1.2 } }

wildcard_constraints:
	bullseye_type     = '|'.join(bullseye_params.keys()),

rule bullseye_find_edit_sites:
	input:
		matrix_cond = os.path.join(config["stamp_dir"], '{sample1_id}.matrix.gz'),
		matrix_ctrl = os.path.join(config["stamp_dir"], '{sample2_id}.matrix.gz'),
		ref_flat    = os.path.join(config["reference_dir"], f'{config["database"]}', f'{config["assembly"]}.refFlat'),
		genome      = os.path.join(config["reference_dir"], f'{config["database"]}', f'{config["assembly"]}.fasta'),
	output: os.path.join(config["stamp_dir"], 'sample_{bullseye_type}_{sample1_id}_vs_{sample2_id}.bed')
	params:
		minimum_editing_rate   = lambda wildcards: bullseye_params[wildcards.bullseye_type]['minimum_rate'],
		maximum_editing_rate   = 90,
		minimum_edit_threshold = lambda wildcards: bullseye_params[wildcards.bullseye_type]['edit_threshold'],
		minimum_edit_sites     = 3,
	conda: '../envs/bullseye.yml'
	threads: 8
	shell:
		'''
			Find_edit_site.pl --annotationFile    {input.ref_flat}                \
			                  --EditedMatrix      {input.matrix_cond}             \
 			                  --controlMatrix     {input.matrix_ctrl}             \
 			                  --minEdit           {params.minimum_editing_rate}   \
 			                  --maxEdit           {params.maximum_editing_rate}   \
 			                  --editFoldThreshold {params.minimum_edit_threshold} \
 			                  --MinEditSites      {params.minimum_edit_sites}     \
 			                  --cpu               {threads}                       \
 			                  --outfile           {output}                        \
 			                  --fallback          {input.genome}                  \
 			                  --verbose
		'''

rule bullseye_merge_edit_sites_complex:
	input: lambda wildcards: [ os.path.join(config["stamp_dir"], 'sample_{bullseye_type}_{sample_id}_vs_'+f'{sample_id2matching_ids[wildcards.sample_id][condition]}.bed') for condition in config["complex_comparisons"][wildcards.named_comparison][1] ]
	output: os.path.join(config["stamp_dir"], 'complex_{bullseye_type}_condition_{named_comparison}_{sample_id}.bed')
	params:
		minimum_reps      = lambda wildcards: len(config["complex_comparisons"][wildcards.named_comparison][1]),
		minimum_mutations = 3,
	conda: '../envs/bullseye.yml'
	shell:
		'''
			perl $(which summarize_sites.pl) --MinRep {params.minimum_reps}      \
			                                 --mut    {params.minimum_mutations} \
			                                 --repOnly                           \
			                                 {input} > {output}
		'''

rule bullseye_filter_edit_sites_simple:
	input: lambda wildcards: [os.path.join(config["stamp_dir"], 'sample_{bullseye_type}_'+f'{sample_id}_vs_{sample_id2matching_ids[sample_id][wildcards.condition2]}.bed') for sample_id in condition2sample_ids[wildcards.condition1]]
	output:                   os.path.join(config["stamp_dir"], 'simple_{bullseye_type}_condition_{condition1}_vs_{condition2}.bed')
	params:
		minimum_reps      = lambda wildcards: len(condition2sample_ids[wildcards.condition1]),
		minimum_mutations = 3,
	conda: '../envs/bullseye.yml'
	shell:
		'''
			perl $(which summarize_sites.pl) --MinRep {params.minimum_reps}      \
			                                 --mut    {params.minimum_mutations} \
			                                 {input} > {output}
		'''

use rule bullseye_filter_edit_sites_simple as rule bullseye_filter_edit_sites_complex with:
	input: lambda wildcards: [ os.path.join(config["stamp_dir"], 'complex_{bullseye_type}_condition_{named_comparison}_'+f'{sample_id}.bed') for sample_id in condition2sample_ids[config["complex_comparisons"][wildcards.named_comparison][0]] ]
	output: os.path.join(config["stamp_dir"], 'complex_{bullseye_type}_condition_{named_comparison}.bed')
	params:
		minimum_reps      = lambda wildcards: len(condition2sample_ids[config["complex_comparisons"][wildcards.named_comparison][0]]),
		minimum_mutations = 3,

rule bullseye2tsv_simple:
	input:
		gene_names = os.path.join(config["reference_dir"], config["database"], f'{config["assembly"]}_gene_id2gene_name.tsv'),
		bullseye   = os.path.join(config["stamp_dir"], 'simple_{bullseye_type}_condition_{condition1}_vs_{condition2}.bed'),
	output:          os.path.join(config["stamp_dir"], 'simple_{bullseye_type}_condition_{condition1}_vs_{condition2}.tsv'),
	conda: '../envs/biocondapy.yml'
	script: '../scripts/bullseye2tsv.py'

use rule bullseye2tsv_simple as rule bullseye2tsv_complex with:
	input:
		gene_names = os.path.join(config["reference_dir"], config["database"], f'{config["assembly"]}_gene_id2gene_name.tsv'),
		bullseye   = os.path.join(config["stamp_dir"], 'complex_{bullseye_type}_condition_{named_comparison}.bed'),
	output:          os.path.join(config["stamp_dir"], 'complex_{bullseye_type}_condition_{named_comparison}.tsv'),

rule bullseye_find_edited_genes_simple:
	input:  os.path.join(config["stamp_dir"], 'simple_relaxed_condition_{condition1}_vs_{condition2}.tsv'),
	output: os.path.join(config["stamp_dir"], 'simple_relaxed_condition_{condition1}_vs_{condition2}_edited_genes.tsv'),
	conda: '../envs/biocondapy.yml'
	script: '../scripts/bullseye_multisites.py'

use rule bullseye_find_edited_genes_simple as rule bullseye_find_edited_genes_complex with:
	input:  os.path.join(config["stamp_dir"], 'complex_relaxed_condition_{named_comparison}.tsv'),
	output: os.path.join(config["stamp_dir"], 'complex_relaxed_condition_{named_comparison}_edited_genes.tsv'),

rule bullseye_pileups:
	input:
		simple_bullseye             = [os.path.join(config["stamp_dir"], f'simple_normal_condition_{condition1}_vs_{condition2}.tsv') for condition1, condition2 in config["simple_comparisons"]],
		complex_gbullseye           = [os.path.join(config["stamp_dir"], f'complex_normal_condition_{name}.tsv') for name in config["complex_comparisons"]],
		multisite_simple_bullseye   = [os.path.join(config["stamp_dir"], f'simple_relaxed_condition_{condition1}_vs_{condition2}_edited_genes.tsv') for condition1, condition2 in config["simple_comparisons"]],
		multisite_complex_gbullseye = [os.path.join(config["stamp_dir"], f'complex_relaxed_condition_{name}_edited_genes.tsv') for name in config["complex_comparisons"]],
		counts                      = [os.path.join(config["stamp_dir"], f'{s}.matrix.gz') for condition in condition2sample_ids for s in conditions2sample_ids[condition] ],
		gtf                         = os.path.join(config["reference_dir"], f'{config["database"]}', f'{config["assembly"]}.gtf'),
		genome                      = os.path.join(config["reference_dir"], f'{config["database"]}', f'{config["assembly"]}.fasta'),
		genome_fai                  = os.path.join(config["reference_dir"], f'{config["database"]}', f'{config["assembly"]}.fasta.fai'),
	output:
		plot_dir = directory(os.path.join(config["stamp_dir"], 'plots')),
		flag     =     touch(os.path.join(config["stamp_dir"], 'plots.flag'))
	params:
		conditions = [ condition for condition in condition2sample_ids for s in conditions2sample_ids[condition] ],
		condition_order = config["display_order"]
		gtf_database = config["database"]
	conda: '../envs/pileups.yml'
	script: '../scripts/stamp_plots.R'

