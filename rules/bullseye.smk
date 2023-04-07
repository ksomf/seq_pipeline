rule bullseye_parse_bam:
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

rule bullseye_find_edit_sites:
	input:
		matrix_cond = os.path.join(config["stamp_dir"], '{sample1_id}.matrix.gz'),
		matrix_ctrl = os.path.join(config["stamp_dir"], '{sample2_id}.matrix.gz'),
		ref_flat    = os.path.join(config["reference_dir"], f'{config["database"]}', f'{config["assembly"]}.refFlat'),
		genome      = os.path.join(config["reference_dir"], f'{config["database"]}', f'{config["assembly"]}.fasta'),
	output: os.path.join(config["stamp_dir"], 'sample_{sample1_id}_vs_{sample2_id}.bed')
	params:
		minimum_editing_rate   = 5,
		maximum_editing_rate   = 90,
		minimum_edit_threshold = 1.5,
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
	input: lambda wildcards: [ os.path.join(config["stamp_dir"], 'sample_{sample_id}_vs_'+f'{sample_id2matching_ids[wildcards.sample_id][condition]}.bed') for condition in config["complex_comparisons"][wildcards.named_comparison][1] ]
	output: os.path.join(config["stamp_dir"], 'complex_condition_{named_comparison}_{sample_id}.bed')
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
	input: lambda wildcards: [os.path.join(config["stamp_dir"], f'sample_{sample_id}_vs_{sample_id2matching_ids[sample_id][wildcards.condition2]}.bed') for sample_id in condition2sample_ids[wildcards.condition1]]
	output: os.path.join(config["stamp_dir"], 'simple_condition_{condition1}_vs_{condition2}.bed')
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
	input: lambda wildcards: [ os.path.join(config["stamp_dir"], 'complex_condition_{named_comparison}_'+f'{sample_id}.bed') for sample_id in condition2sample_ids[config["complex_comparisons"][wildcards.named_comparison][0]] ]
	output: os.path.join(config["stamp_dir"], 'complex_condition_{named_comparison}.bed')
	params:
		minimum_reps      = lambda wildcards: len(condition2sample_ids[config["complex_comparisons"][wildcards.named_comparison][0]]),
		minimum_mutations = 3,
